"""Synthetic metadata/state tests, not biological performance benchmarks."""
import csv
import hashlib
import json
from pathlib import Path
import re
import tempfile
import unittest

import yaml

from micpscope.application import catalog
from micpscope.analysis import characterize
from micpscope.cli import main
from micpscope.knowledge import enrich
from micpscope.reference import attach, apply_category_evidence
from micpscope.report import render_report, assessment_summary

ROOT = Path(__file__).resolve().parents[1]


def knowledge():
    return dict(categories=[dict(id=c, name=c, description='Synthetic grouping')
                            for c in ('one', 'two', 'unassigned')], sources=[], entries=[])


def evidence(proteins, status='reference_candidate', categories=None, gene=''):
    return dict(schema_version='1.0', input_sha256='fixture', knowledge_catalog=knowledge(),
                database={'version': 'fixture'}, references={'r': dict(id='r', gene=gene,
                protein_name='Synthetic candidate', categories=categories or ['two'])},
                queries={p['protein_id']: dict(status=status, hits=[dict(reference_id='r',
                identity=99, query_coverage=100, reference_coverage=100, evalue=1e-10, bitscore=100)]) for p in proteins},
                interpretation='Synthetic fixture', tool_version='fixture', thresholds={})


class EvidenceIntegrityTests(unittest.TestCase):
    def test_top_summary_counts_unresolved_states_without_closed_details(self):
        states = ['annotation_reference_conflict', 'ambiguous_candidates', 'ambiguous_candidates',
                  'unresolved_reference_category', 'reference_candidate', 'no_qualifying_hit']
        result = dict(proteins=[dict(protein_id=str(i), status='calculated', evidence_status=s)
                                for i, s in enumerate(states)], reference_search={'queries': {}})
        html = assessment_summary(result)
        visible = html.split('<details>', 1)[0]
        self.assertIn('Unresolved reference evidence: 4 proteins.', visible)
        self.assertIn('1 annotation/reference conflicts; 2 ambiguous candidate sets; 1 unknown reference categories.', visible)
        self.assertIn('Unresolved evidence is not a negative finding', visible)

    def test_every_starter_record_resolves_with_neutral_query_id(self):
        spec = yaml.safe_load((ROOT / 'database_builder/references/micp-starter.yaml').read_text())
        k = catalog(spec)
        self.assertEqual(len(spec['proteins']), 164)
        for ref in spec['proteins']:
            with self.subTest(reference=ref['id']):
                p = enrich([characterize('neutral', 'neutral', 'AVIL')], k)[0]
                result = evidence([p], categories=ref['categories'])
                result['knowledge_catalog'] = k
                result['references']['r'] = dict(ref, id='r')
                apply_category_evidence(p, result)
                self.assertEqual(p['category'], ref['categories'][0])
                self.assertEqual(p['reference_categories'], ref['categories'])
                self.assertEqual(p['evidence_status'], 'reference_candidate')
                self.assertIn('first category', p['category_basis'])

    def test_unknown_category_rejected_before_build_or_search(self):
        spec = dict(schema_version='1.0', name='fixture', version='1', proteins=[{'id': 'r', 'categories': ['typo']}],
                    annotations={'knowledge_catalog': knowledge()})
        with self.assertRaisesRegex(ValueError, 'reference categories'):
            catalog(spec)
        with tempfile.TemporaryDirectory() as td:
            p = enrich([characterize('p', 'p', 'AVIL')], knowledge())[0]
            result = evidence([p], categories=['typo'])
            path = Path(td) / 'search.json'; path.write_text(json.dumps(result))
            with self.assertRaisesRegex(ValueError, 'reference categories'):
                attach([p], path, 'fixture')

    def test_conflict_retains_both_sources_and_excludes_summary_group(self):
        p = enrich([characterize('p', 'p gene=first', 'AVIL')], knowledge())[0]
        p['supplied_annotation']['category'] = 'one'
        result = evidence([p], gene='second')
        apply_category_evidence(p, result)
        self.assertEqual(p['category'], 'unassigned')
        self.assertEqual(p['annotation_category'], 'one')
        self.assertEqual(p['reference_categories'], ['two'])
        self.assertEqual(p['gene_name'], 'first')
        self.assertEqual(p['reference_gene'], 'second')
        self.assertEqual(p['evidence_status'], 'annotation_reference_conflict')
        self.assertEqual(len(p['evidence_conflicts']), 2)

    def test_ambiguous_and_skipped_states_never_count_as_resolved(self):
        for status in ('ambiguous_candidates', 'skipped_unsupported_residues', 'no_qualifying_hit'):
            p = enrich([characterize('p', 'p', 'AVIL')], knowledge())[0]
            result = evidence([p], status=status)
            if status != 'ambiguous_candidates':
                result['queries']['p']['hits'] = []
            apply_category_evidence(p, result)
            self.assertEqual(p['category'], 'unassigned')
            self.assertEqual(p['evidence_status'], status)

    def test_competing_categories_are_unresolved_even_with_same_gene(self):
        p = enrich([characterize('p', 'p', 'AVIL')], knowledge())[0]
        result = evidence([p])
        result['references']['r2'] = dict(result['references']['r'], id='r2', categories=['one'])
        result['queries']['p']['hits'].append(dict(result['queries']['p']['hits'][0], reference_id='r2', bitscore=98))
        apply_category_evidence(p, result)
        self.assertEqual(p['evidence_status'], 'ambiguous_candidates')
        self.assertEqual(p['category'], 'unassigned')
        self.assertEqual(set(p['reference_categories']), {'one', 'two'})

    def test_json_tsv_html_and_visual_payload_agree(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            fasta = root / 'p.faa'; fasta.write_text('>p gene=first\nAVIL\n')
            p = enrich([characterize('p', 'p gene=first', 'AVIL')], knowledge())[0]
            result = evidence([p], gene='second')
            result['input_sha256'] = hashlib.sha256(fasta.read_bytes()).hexdigest()
            search = root / 'search.json'; search.write_text(json.dumps(result))
            main(['--proteins', str(fasta), '--sample', 'fixture', '--reference-results', str(search), '--out', str(root / 'out')])
            saved = json.loads((root / 'out/results.json').read_text())
            self.assertEqual(saved['proteins'][0]['evidence_status'], 'annotation_reference_conflict')
            with (root / 'out/proteins.tsv').open() as handle:
                row = next(csv.DictReader(handle, delimiter='\t'))
            self.assertEqual(row['category'], 'unassigned')
            self.assertEqual(row['evidence_status'], 'annotation_reference_conflict')
            html = (root / 'out/report.html').read_text()
            self.assertIn('Annotation / reference conflict · unresolved', html)
            payload = json.loads(re.search(r'<script id="gx-data" type="application/json">(.*?)</script>', html, re.S)[1])
            self.assertEqual(payload['proteins'][0]['evidence_status'], row['evidence_status'])
            self.assertEqual(payload['proteins'][0]['category'], 'unassigned')
            # Refresh of old saved results must also expose disagreement, without mutating data.
            saved['proteins'][0]['category'] = 'one'
            before = json.dumps(saved, sort_keys=True)
            self.assertIn('Annotation / reference conflict · unresolved', render_report(saved))
            self.assertEqual(json.dumps(saved, sort_keys=True), before)
