import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from micpscope.reference import search, verify_database, attach
from micpscope.cli import main


def fixture(root):
    db = root / 'db'; db.mkdir()
    spec = dict(schema_version='1.0', name='synthetic', version='test', proteins=[
        dict(id='p1', gene='one', protein_name='<candidate>', references=['paper']),
        dict(id='p2', gene='two', protein_name='alternative', references=[])],
        annotations={'literature_records':[dict(id='paper', fields=[dict(heading='Function',text='<script>')])]},
        references=[dict(id='paper',citation='Synthetic reference')])
    for name, content in {'proteins.faa':'>p1\nACDEFGHIKLMNPQRSTVWY\n', 'proteins.dmnd':'mock index',
                          'diamond-version.txt':'mock', 'reference-manifest.json':json.dumps(dict(schema_version='1.0',specification=spec))}.items():
        (db/name).write_text(content)
    (db/'checksums.json').write_text(json.dumps({p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in db.iterdir()}))
    (db/'completed.json').write_text(json.dumps(dict(status='completed',checksums_sha256=hashlib.sha256((db/'checksums.json').read_bytes()).hexdigest())))
    return db


class ReferenceTests(unittest.TestCase):
    def test_ambiguity_coverage_skips_provenance_and_report(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);db=fixture(root)
            proteins=root/'input.faa';proteins.write_text('>sample|id\nACDEFGHIKLMNPQRSTVWY*\n>bad\nAX\n>none\nACDEFGHIKLMNPQRSTVWY\n')
            def run(command, check):
                Path(command[command.index('--out')+1]).write_text(
                    'q0\tp1\t90\t1\t20\t1\t20\t20\t20\t1e-20\t100\n'
                    'q0\tp2\t89\t1\t20\t1\t20\t20\t20\t1e-19\t97\n'
                    'q2\tp1\t100\t1\t5\t1\t5\t20\t20\t1e-20\t200\n')
            output=root/'search.json'
            with patch('micpscope.reference.subprocess.check_output',return_value='mock diamond'), patch('micpscope.reference.subprocess.run',side_effect=run):
                result=search(proteins,db,output)
            self.assertEqual(result['queries']['sample|id']['status'],'ambiguous_candidates')
            self.assertEqual(result['queries']['bad']['status'],'skipped_unsupported_residues')
            self.assertEqual(result['queries']['none']['status'],'no_qualifying_hit')
            self.assertEqual(len(result['queries']['sample|id']['hits']),2)
            main(['--proteins',str(proteins),'--sample','test','--reference-results',str(output),'--out',str(root/'report')])
            html=(root/'report/report.html').read_text()
            self.assertIn('Reference search · 1 proteins with candidates',html)
            self.assertIn('&lt;candidate&gt;',html)
            self.assertIn('&lt;script&gt;',html)
            self.assertIn('reference_search_status',(root/'report/proteins.tsv').read_text())
            with self.assertRaises(ValueError): attach([],output,'wrong')
            (db/'proteins.dmnd').write_text('tampered')
            with self.assertRaises(ValueError):verify_database(db)

    def test_pinned_identity_and_invalid_threshold(self):
        with tempfile.TemporaryDirectory() as td:
            db=fixture(Path(td))
            with self.assertRaises(ValueError):verify_database(db,'wrong')
            with self.assertRaises(ValueError):search('missing',db,'output',identity=float('nan'))

if __name__=='__main__':unittest.main()
