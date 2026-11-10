"""Controller identity must describe executing code, including Conda builds."""
import importlib.metadata
import sys
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from micpscope.artifact_retention import runtime_version


class ControllerIdentityTests(unittest.TestCase):
    def check_version(self, runtime, metadata):
        with patch.dict(sys.modules, {'snakemake': SimpleNamespace(__version__=runtime)}):
            with patch('micpscope.artifact_retention.importlib.metadata.version', return_value=metadata):
                return runtime_version()

    def test_conda_placeholder_accepts_only_supported_runtime(self):
        for version in ('9.6.2', '9.27.1'):
            self.assertEqual(self.check_version(version, '0.0.0'), version)
            self.assertEqual(self.check_version(version, version), version)

    def test_missing_distribution_metadata_uses_known_runtime(self):
        with patch.dict(sys.modules, {'snakemake': SimpleNamespace(__version__='9.6.2')}):
            with patch('micpscope.artifact_retention.importlib.metadata.version', side_effect=importlib.metadata.PackageNotFoundError('snakemake')):
                self.assertEqual(runtime_version(), '9.6.2')

    def test_unknown_runtime_never_falls_back_to_supported_distribution(self):
        for version in (None, '', '0.0.0', '9.6.3', 'unknown'):
            with self.subTest(version=version), self.assertRaisesRegex(RuntimeError, 'runtime'):
                self.check_version(version, '9.6.2')

    def test_concrete_conflicting_distribution_fails_closed(self):
        for metadata in ('9.27.1', '9.6.3', 'unknown', ''):
            with self.subTest(metadata=metadata), self.assertRaisesRegex(RuntimeError, 'Conflicting'):
                self.check_version('9.6.2', metadata)


if __name__ == '__main__':
    unittest.main()
