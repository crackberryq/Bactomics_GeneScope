import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from micpscope.application import application, archive_name, catalog
from micpscope.cli import main
import test_database_builder as db
builder = db.builder
import yaml


class ApplicationTests(unittest.TestCase):
    def test_custom_application_build_and_report_without_micp_defaults(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td); path=db.DatabaseBuilderTests().fixture(root)
            spec=yaml.safe_load(path.read_text());spec.pop('name');spec.pop('version')
            spec.update(schema_version='2.0',application=dict(id='custom',name='Custom <engineering>',version='2.1',description='Independent application fixture'))
            path.write_text(yaml.safe_dump(spec));builder.validate(path,root/'plan.json');builder.retrieve(root/'plan.json',root/'seq')
            (root/'index.dmnd').write_text('fixture');(root/'index.version.txt').write_text('fixture')
            builder.package(root/'plan.json',root/'seq',root/'index.dmnd',root/'db');builder.verify(root/'db',root/'db/completed.json')
            self.assertTrue((root/'db/custom-2.1.zip').is_file())
            receipt=json.loads((root/'db/completed.json').read_text());self.assertEqual(receipt['application'],spec['application'])
            self.assertEqual(receipt['archive_name'],'custom-2.1.zip')
            proteins=root/'query.faa';proteins.write_text('>p gene=ureC supplied label\nACDE\n')
            result=dict(schema_version='1.0',input_sha256=hashlib.sha256(proteins.read_bytes()).hexdigest(),application=spec['application'],
                knowledge_catalog=catalog(spec),queries={'p':{'status':'no_qualifying_hit','hits':[]}},database=spec['application'],
                references={},interpretation='Synthetic fixture',tool_version='fixture',thresholds={},bibliography=[])
            (root/'search.json').write_text(json.dumps(result))
            main(['--proteins',str(proteins),'--sample','fixture','--reference-results',str(root/'search.json'),'--out',str(root/'out')])
            data=json.loads((root/'out/results.json').read_text());self.assertEqual(data['application'],spec['application'])
            self.assertEqual(data['proteins'][0]['category'],'unassigned')
            html=(root/'out/report.html').read_text();self.assertIn('GeneScope — Custom &lt;engineering&gt;',html)
            self.assertNotIn('MICP',html);self.assertNotIn('Catalytic component of the multimeric',html)

    def test_legacy_identity_and_invalid_schema(self):
        old={'schema_version':'1.0','name':'legacy','version':'1'}
        self.assertEqual(application(old)['name'],'legacy')
        self.assertEqual(archive_name(old),'legacy-1.zip')
        for spec in ({'schema_version':'3.0'}, {'schema_version':'2.0'}, {'schema_version':'2.0','application':{'id':'../x','name':'Name','version':'1','description':'Description'}}):
            with self.assertRaises(ValueError):application(spec)

    def test_unsafe_catalog_links_rejected(self):
        spec={'schema_version':'2.0','application':{'id':'test','name':'Test','version':'1','description':'Test'},
              'annotations':{'knowledge_catalog':{'categories':[{'id':'unassigned','name':'Unknown','description':''}],
                'entries':[],'sources':[{'id':'bad','title':'Bad','type':'test','url':'javascript:alert(1)'}]}}}
        with self.assertRaises(ValueError):catalog(spec)


class ReportCoverageTests(unittest.TestCase):
    def test_descriptor_only_is_not_a_negative_application_result(self):
        from micpscope.report import assessment_summary
        html = assessment_summary({'proteins': [], 'evidence': []})
        self.assertIn('Descriptor-only run', html)
        self.assertIn('Unassigned proteins are not negative findings', html)
        self.assertIn('existing frozen run', html)

    def test_zero_matches_differs_from_search_not_requested(self):
        from micpscope.report import assessment_summary
        html = assessment_summary({'proteins': [], 'application': {'name': '<custom>'},
            'reference_search': {'queries': {'q': {'hits': []}}}})
        self.assertIn('0 proteins have qualifying reference candidates', html)
        self.assertNotIn('Descriptor-only run', html)
        self.assertIn('&lt;custom&gt;', html)
