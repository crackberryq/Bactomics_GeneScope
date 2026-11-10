import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from micpscope.retention import publish, retained_attempt
from micpscope.reference import search
from test_reference import fixture


class RetentionTests(unittest.TestCase):
    def test_directory_publication_preserves_native_extras_without_carrying_stale_results(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            current = root / "current"; current.mkdir()
            (current / "native.extra").write_text("old native artifact")
            generated = root / "generated"; generated.mkdir()
            (generated / "new.result").write_text("new result")
            publish(generated, current)
            self.assertFalse((current / "native.extra").exists())
            self.assertEqual((current / "new.result").read_text(), "new result")
            previous = list((root / "retained/publication").glob("*/previous/native.extra"))
            self.assertEqual([p.read_text() for p in previous], ["old native artifact"])
            self.assertEqual((generated / "new.result").read_text(), "new result")

    def test_publication_preserves_current_file_and_durable_source(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            old = root / 'report.html'; old.write_text('old')
            generated = root / 'generated'; generated.write_text('new')
            publish(generated, old)
            self.assertEqual(old.read_text(), 'new')
            self.assertEqual(generated.read_text(), 'new')
            previous = list((root / 'retained/publication').glob('*/previous'))
            self.assertEqual([p.read_text() for p in previous], ['old'])

    def test_publication_copy_failure_never_removes_original(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            old = root / 'result'; old.write_text('old')
            generated = root / 'generated'; generated.write_text('new')
            with patch('micpscope.retention._copy', side_effect=OSError('disk full')):
                with self.assertRaises(OSError):
                    publish(generated, old)
            self.assertEqual(old.read_text(), 'old')
            self.assertTrue(list((root / 'retained/publication').glob('*/failed.json')))

    def test_failed_attempt_leaves_generated_data_and_failure_receipt(self):
        with tempfile.TemporaryDirectory() as td:
            with self.assertRaisesRegex(ValueError, 'failure'):
                with retained_attempt(td, 'test') as attempt:
                    (attempt / 'partial.tsv').write_text('partial')
                    raise ValueError('failure')
            self.assertEqual((attempt / 'partial.tsv').read_text(), 'partial')
            self.assertEqual(json.loads((attempt / 'failed.json').read_text())['status'], 'failed')
            self.assertFalse((attempt / 'completed.json').exists())

    def test_reference_failure_and_retry_keep_native_query_hits_and_logs(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td); db = fixture(root)
            proteins = root / 'input.faa'; proteins.write_text('>p\nACDEFGHIKLMNPQRSTVWY\n')
            output = root / 'search.json'; output.write_text('old result')
            def failed(command, **kwargs):
                Path(command[command.index('--out') + 1]).write_text('partial native hits')
                self.assertIn('--no-unlink', command)
                (Path(command[command.index('--tmpdir') + 1]) / 'native.tmp').write_text('native intermediate')
                kwargs['stderr'].write(b'failure detail')
                raise subprocess.CalledProcessError(1, command)
            with patch('micpscope.reference.subprocess.check_output', return_value='test version'), patch('micpscope.reference.subprocess.run', side_effect=failed):
                with self.assertRaises(subprocess.CalledProcessError):
                    search(proteins, db, output)
            self.assertEqual(output.read_text(), 'old result')
            attempts = list((root / 'retained/reference-search').iterdir())
            self.assertEqual(len(attempts), 1)
            self.assertEqual((attempts[0] / 'hits.tsv').read_text(), 'partial native hits')
            self.assertEqual((attempts[0] / 'stderr.log').read_text(), 'failure detail')
            self.assertTrue((attempts[0] / 'queries.faa').is_file())
            self.assertEqual((attempts[0] / 'native-tmp/native.tmp').read_text(), 'native intermediate')
            def success(command, **kwargs):
                Path(command[command.index('--out') + 1]).write_text('')
            with patch('micpscope.reference.subprocess.check_output', return_value='test version'), patch('micpscope.reference.subprocess.run', side_effect=success):
                search(proteins, db, output)
            self.assertEqual(len(list((root / 'retained/reference-search').iterdir())), 2)
            self.assertEqual((attempts[0] / 'hits.tsv').read_text(), 'partial native hits')
            self.assertIn('old result', [p.read_text() for p in (root / 'retained/publication').glob('*/previous')])


if __name__ == '__main__':
    unittest.main()
