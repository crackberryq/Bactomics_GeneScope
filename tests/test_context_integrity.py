from pathlib import Path
import json
import re
import tempfile
import unittest

from micpscope.analysis import characterize
from micpscope.context import add_context, neighborhood
from micpscope.explorer import explorer


class ContextIntegrityTests(unittest.TestCase):
    def parse(self, text):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / 'p.gff'; path.write_text(text)
            p = characterize('p', 'p', 'AVIL')
            return p, add_context([p], path)

    def test_impossible_bounds_and_conflicting_declarations_reject(self):
        for text in (
            '##sequence-region c 1 10\nc\tf\tCDS\t20\t31\t.\t+\t0\tID=p\n',
            '##sequence-region c 1 100\n##sequence-region c 1 110\n',
        ):
            with self.subTest(text=text), self.assertRaises(ValueError):
                self.parse(text)

    def test_bounds_checked_even_when_directive_follows_feature(self):
        with self.assertRaisesRegex(ValueError, 'outside declared'):
            self.parse('c\tf\tCDS\t20\t31\t.\t+\t0\tID=p\n##sequence-region c 1 10\n')

    def test_circular_origin_and_multipart_are_retained(self):
        p, context = self.parse('##sequence-region c 1 100\nc\tf\tregion\t1\t100\t.\t+\t.\tID=c;Is_circular=true\n'
                               'c\tf\tCDS\t95\t106\t.\t+\t0\tID=p\n')
        self.assertEqual(p['coordinates'][0]['display_segments'], [[95, 100], [1, 6]])
        self.assertEqual(neighborhood(p, context), [])
        self.assertIn('sequence pairing unverified', p['coordinate_linkage'])
        self.assertTrue(context['warnings'])
        html = explorer({'context': context, 'proteins': [p]})
        payload = json.loads(re.search(r'<script id="gx-data" type="application/json">(.*?)</script>', html, re.S)[1])
        self.assertEqual(len(payload['features']), 2)
        self.assertEqual(len({f['source_feature_key'] for f in payload['features']}), 1)
        p, context = self.parse('##sequence-region c 1 100\nc\tf\tCDS\t90\t100\t.\t+\t0\tID=p\n'
                               'c\tf\tCDS\t1\t10\t.\t+\t0\tID=p\n')
        self.assertEqual(len(p['coordinates']), 2)
        self.assertIn('multiple CDS', p['coordinate_status'])

    def test_missing_bounds_and_exact_id_do_not_claim_pair_verification(self):
        p, context = self.parse('c\tf\tCDS\t1\t12\t.\t+\t0\tID=p\n')
        self.assertIn('identifier-only', p['coordinate_linkage'])
        self.assertTrue(any('bounds are unverified' in w for w in context['warnings']))

    def test_partial_region_declarations_allow_unbounded_contig_with_warning(self):
        p, context = self.parse('##sequence-region c 1 100\nother\tf\tCDS\t1\t12\t.\t+\t0\tID=p\n')
        self.assertEqual(p['coordinates'][0]['contig'], 'other')
        self.assertTrue(any('contig other' in w and 'bounds are unverified' in w for w in context['warnings']))
