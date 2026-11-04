import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from micpscope.analysis import characterize
from micpscope.cli import legacy_main as main
from micpscope.knowledge import annotation_context, enrich, interpret_descriptors, load_catalog


class KnowledgeTests(unittest.TestCase):
    def setUp(self):
        self.catalog, self.digest = load_catalog()

    def test_gene_product_and_uniprot_headers(self):
        for header in ('p [gene=ureC] urease subunit alpha',
                       'p Urease subunit alpha', 'p Urease OS=Example GN=ureC PE=3'):
            row = annotation_context('p', header, self.catalog)
            self.assertEqual(row['gene_name'], 'ureC')
            self.assertEqual(row['category'], 'urease')
            self.assertIn('not sequence-confirmed', row['annotation_basis'])

    def test_unknown_and_ambiguous_never_claim_detection(self):
        for identifier, header in [('contig_ureC_1', 'contig_ureC_1 hypothetical protein'),
                                   ('p', 'p ureC-like unrelated annotation'), ('p', 'p')]:
            self.assertEqual(annotation_context(identifier, header, self.catalog)['category'], 'unassigned')
        row = annotation_context('p', 'p ureC and ureD', self.catalog)
        self.assertEqual(row['annotation_status'], 'ambiguous annotation')
        self.assertEqual(row['category'], 'unassigned')

    def test_extreme_descriptors_do_not_produce_optima(self):
        rows = enrich([characterize('p','p [gene=ureC]', 'AVIL')], self.catalog)
        self.assertGreater(rows[0]['aliphatic_index'], 200)
        self.assertTrue(rows[0]['temperature_indication'].startswith('Undetermined'))
        self.assertTrue(rows[0]['ph_indication'].startswith('Undetermined'))
        for n, label in [(39.99,'below 40'), (40,'40: boundary'), (40.01,'above 40')]:
            p = dict(status='calculated', instability_index=n, gravy=-1)
            self.assertTrue(interpret_descriptors(p)['instability_class'].startswith(label))
        self.assertEqual(interpret_descriptors(dict(status='not_calculated'))['instability_class'], 'not calculated')

    def test_catalog_references_and_categories_resolve(self):
        categories={c['id'] for c in self.catalog['categories']}
        refs={r['id'] for r in self.catalog['sources']}
        self.assertEqual(len({e['id'] for e in self.catalog['entries']}), len(self.catalog['entries']))
        for entry in self.catalog['entries']:
            self.assertIn(entry['category'], categories)
            self.assertTrue(set(entry['references']) <= refs)
            self.assertTrue(entry['scope'])
            self.assertTrue(entry['evidence_level'])

    def test_report_json_tsv_agree_and_escape_annotation(self):
        with tempfile.TemporaryDirectory() as directory:
            base=Path(directory)
            fasta=base/'proteins.faa'
            fasta.write_text('>p [gene=ureC] <img src=x onerror=alert(1)>\nAVIL\n>q hypothetical protein\nMX\n')
            main(['--proteins',str(fasta),'--sample','demo','--out',str(base/'out')])
            result=json.loads((base/'out/results.json').read_text())
            self.assertEqual(result['proteins'][0]['gene_name'], 'ureC')
            self.assertEqual(result['proteins'][1]['category'], 'unassigned')
            self.assertEqual(result['provenance']['catalog_sha256'], self.digest)
            report=(base/'out/report.html').read_text()
            self.assertIn('Urease: catalysis', report)
            self.assertIn('Best operating temperature: undetermined', report)
            self.assertIn('urease2020', report)
            self.assertNotIn('<img src=x',report)
            self.assertIn('&lt;img src=x',report)
            tsv=(base/'out/proteins.tsv').read_text()
            self.assertIn('gene_name\tprotein_name\tfunction',tsv)
            self.assertIn('ureC',tsv)


if __name__ == '__main__':
    unittest.main()

