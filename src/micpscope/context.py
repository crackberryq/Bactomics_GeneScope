"""Read recorded genomic coordinates and protein signatures without inferring function."""
import csv
import hashlib
from collections import defaultdict
from urllib.parse import unquote


def add_context(proteins, gff=None, interpro=None):
    by_id = {p['protein_id']: p for p in proteins}
    context = {'coordinates': 'not supplied', 'domains': 'not supplied', 'warnings': [],
               'input_sha256': {}, 'features': [], 'operons': [], 'sequence_regions': {}}
    copies = defaultdict(list)
    for p in proteins:
        p['domains'] = []
        p['coordinates'] = []
        p['domain_status'] = 'not assessed'
        p['coordinate_status'] = 'not assessed'
        copies[p['sequence_sha256']].append(p['protein_id'])
    for p in proteins:
        p['identical_sequence_records'] = copies[p['sequence_sha256']]
    if gff:
        context['input_sha256']['gff'] = hashlib.sha256(gff.read_bytes()).hexdigest()
        features = []
        for number, line in enumerate(gff.read_text(encoding='utf-8-sig').splitlines(), 1):
            if line.startswith('##FASTA'):
                break
            if line.startswith('##sequence-region '):
                parts = line.split()
                if len(parts) != 4 or not 1 <= int(parts[2]) <= int(parts[3]):
                    raise ValueError(f'GFF3 row {number}: invalid sequence-region')
                context['sequence_regions'][parts[1]] = [int(parts[2]), int(parts[3])]
            if not line or line.startswith('#'):
                continue
            cols = line.split('\t')
            if len(cols) != 9:
                raise ValueError(f'GFF3 row {number}: expected nine fields')
            if cols[2] not in {'CDS', 'operon'}:
                continue
            start, end = int(cols[3]), int(cols[4])
            if start < 1 or end < start or cols[6] not in {'+', '-', '.', '?'}:
                raise ValueError(f'GFF3 row {number}: invalid coordinates or strand')
            attrs = {}
            for item in cols[8].split(';'):
                if '=' in item:
                    key, value = item.split('=', 1)
                    attrs[key] = [unquote(v) for v in value.split(',')]
            candidates = set(v for k in ('ID', 'protein_id', 'locus_tag') for v in attrs.get(k, [])) & by_id.keys()
            feature = dict(contig=cols[0], start=start, end=end, strand=cols[6],
                           feature_id=next(iter(attrs.get('ID', [])), f'CDS-row-{number}'),
                           product=next(iter(attrs.get('product', [])), 'CDS'),
                           gene_name=next(iter(attrs.get('gene', [])), ''),
                           protein_id=next(iter(candidates)) if len(candidates) == 1 else None)
            if cols[2] == 'operon':
                feature['source'] = cols[1]
                feature['product'] = next(iter(attrs.get('Name', [])), 'Supplied operon annotation')
                feature['protein_id'] = None
                context['operons'].append(feature)
                continue
            features.append(feature)
            if len(candidates) == 1:
                by_id[feature['protein_id']]['coordinates'].append(feature)
            elif len(candidates) > 1:
                context['warnings'].append(f'GFF3 row {number}: multiple protein identifiers; left unresolved')
        context['features'] = sorted(features, key=lambda f: (f['contig'], f['start'], f['end']))
        context['coordinates'] = 'imported GFF3 CDS coordinates; no operon prediction'
        for p in proteins:
            p['coordinate_status'] = ('recorded single CDS' if len(p['coordinates']) == 1 else
                                      'multiple CDS segments; neighborhood unresolved' if p['coordinates'] else
                                      'no exact identifier match')
    if interpro:
        context['input_sha256']['interpro'] = hashlib.sha256(interpro.read_bytes()).hexdigest()
        unknown = set()
        seen = set()
        with interpro.open(encoding='utf-8-sig', newline='') as handle:
            for number, cols in enumerate(csv.reader(handle, delimiter='\t'), 1):
                if not cols or cols[0].startswith('#'):
                    continue
                if not 11 <= len(cols) <= 15:
                    raise ValueError(f'InterProScan TSV row {number}: expected 11–15 fields')
                if cols[0] not in by_id:
                    unknown.add(cols[0]); continue
                p = by_id[cols[0]]
                if int(cols[2]) != p['length'] or cols[1].lower() != p['sequence_md5']:
                    raise ValueError(f'InterProScan TSV row {number}: sequence identity mismatch')
                start, end = int(cols[6]), int(cols[7])
                if not 1 <= start <= end <= p['length']:
                    raise ValueError(f'InterProScan TSV row {number}: signature outside protein')
                key = tuple(cols)
                if key in seen:
                    continue
                seen.add(key)
                p['domains'].append(dict(analysis=cols[3], accession=cols[4], description=cols[5],
                                         start=start, end=end, score=cols[8], status=cols[9], date=cols[10],
                                         interpro=cols[11] if len(cols) > 11 else '-'))
        if unknown:
            context['warnings'].append(f'{len(unknown)} InterProScan identifiers did not match this protein input and were ignored')
        context['domains'] = 'imported InterProScan TSV signatures; no new prediction performed'
        for p in proteins:
            p['domain_status'] = 'recorded signature matches' if p['domains'] else 'no row supplied; assessment coverage unknown'
    return context


def neighborhood(protein, context):
    if len(protein.get('coordinates', [])) != 1:
        return []
    focal = protein['coordinates'][0]
    features = [f for f in context.get('features', []) if f['contig'] == focal['contig']]
    index = features.index(focal)
    return features[max(0, index - 2):index + 3]

