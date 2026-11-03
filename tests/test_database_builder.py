import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
import yaml

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'database_builder'))
import build as builder

class DatabaseBuilderTests(unittest.TestCase):
    def test_bundled_starter_builds_offline_with_complete_evidence(self):
        from unittest.mock import patch
        spec=ROOT/'database_builder/references/micp-starter.yaml'
        with tempfile.TemporaryDirectory() as td:
            r=Path(td)
            builder.validate(spec,r/'plan.json')
            with patch.object(builder,'download',side_effect=AssertionError('Bundled build must be offline')):
                builder.retrieve(r/'plan.json',r/'seq')
            plan=json.loads((r/'plan.json').read_text())['spec']
            self.assertEqual(len(builder.fasta(r/'seq/proteins.faa')),164)
            self.assertEqual(len(plan['annotations']['literature_records']),207)
            self.assertEqual({x['source_id'] for x in plan['annotations']['literature_records']},{'processes','regulation'})
            ids={x['id'] for x in plan['annotations']['literature_records']}
            for protein in plan['proteins']:
                self.assertTrue(set(protein['references'][1:])<=ids)
            self.assertEqual((r/'seq/regulatory.fna').read_text(),'')

    def fixture(self,root):
        (root/'protein.faa').write_text('>protein\nACDEFGHIKLMNPQRSTVWY\n')
        (root/'dna.fna').write_text('>contig\nACGTACGT\n')
        spec={'schema_version':'1.0','name':'synthetic','version':'test','proteins':[{'id':'p1','path':str(root/'protein.faa'),'record':'protein','references':['synthetic'],'evidence':'synthetic','license':'test'}], 'regulatory':[{'id':'r1','path':str(root/'dna.fna'),'record':'contig','references':['synthetic'],'evidence':'synthetic','license':'test','feature_type':'upstream_region','accession_version':'synthetic.1','start':2,'end':5,'strand':'-','sequence_sha256':hashlib.sha256(b'TACG').hexdigest()}]}
        p=root/'spec.yaml';p.write_text(yaml.safe_dump(spec)); return p
    def test_local_sequences_and_package_tampering(self):
        with tempfile.TemporaryDirectory() as td:
            r=Path(td);p=self.fixture(r);plan=r/'plan.json'
            builder.validate(p,plan);builder.retrieve(plan,r/'seq')
            self.assertEqual((r/'seq/regulatory.fna').read_text(),'>r1\nTACG\n')
            (r/'proteins.dmnd').write_bytes(b'synthetic-index-only')
            (r/'proteins.version.txt').write_text('synthetic fixture')
            builder.package(plan,r/'seq',r/'proteins.dmnd',r/'db')
            builder.verify(r/'db',r/'receipt.json')
            (r/'db/proteins.faa').write_text('tampered')
            with self.assertRaises(ValueError): builder.verify(r/'db',r/'receipt2.json')
            self.assertFalse((r/'receipt2.json').exists())
    def test_changes_after_validation_and_invalid_coordinates(self):
        with tempfile.TemporaryDirectory() as td:
            r=Path(td);p=self.fixture(r);plan=r/'plan.json';builder.validate(p,plan)
            (r/'protein.faa').write_text('>protein\nAAAA\n')
            with self.assertRaises(ValueError): builder.retrieve(plan,r/'seq')
            spec=yaml.safe_load(p.read_text());spec['regulatory'][0]['end']=99;p.write_text(yaml.safe_dump(spec))
            with self.assertRaises(ValueError): builder.validate(p,plan)
    def test_duplicate_ids(self):
        with tempfile.TemporaryDirectory() as td:
            r=Path(td);p=self.fixture(r);s=yaml.safe_load(p.read_text());s['regulatory'][0]['id']='p1';p.write_text(yaml.safe_dump(s))
            with self.assertRaises(ValueError): builder.validate(p,r/'plan')

if __name__=='__main__': unittest.main()
