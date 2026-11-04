"""Traceable protein similarity search; reference associations are not activity calls."""
import hashlib
import json
import math
from pathlib import Path
import subprocess
import tempfile


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def verify_database(folder, expected=None):
    folder = Path(folder).resolve()
    receipt = json.loads((folder / 'completed.json').read_text())
    identity = digest(folder / 'checksums.json')
    if receipt.get('status') != 'completed' or receipt.get('checksums_sha256') != identity:
        raise ValueError('Reference database completion receipt does not match checksums')
    if expected and expected != identity:
        raise ValueError('Reference database changed since this run was created')
    checksums = json.loads((folder / 'checksums.json').read_text())
    required = {'proteins.faa', 'proteins.dmnd', 'reference-manifest.json', 'diamond-version.txt'}
    if not isinstance(checksums, dict) or not required <= checksums.keys():
        raise ValueError('Reference database is missing required checksums')
    for name, checksum in checksums.items():
        if Path(name).name != name or name in {'.', '..'}:
            raise ValueError('Unsafe reference database filename')
        if digest(folder / name) != checksum:
            raise ValueError('Reference database checksum mismatch: ' + name)
    manifest = json.loads((folder / 'reference-manifest.json').read_text())
    if manifest.get('schema_version') != '1.0':
        raise ValueError('Unsupported reference database schema')
    return manifest, identity


def search(proteins, database, output, threads=1, expected=None,
           identity=30.0, coverage=70.0, evalue=1e-5):
    from .analysis import read_fasta
    for value, lower, upper in ((identity, 0, 100), (coverage, 0, 100), (evalue, 0, 1)):
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or not lower < value <= upper:
            raise ValueError('Invalid reference search threshold')
    if type(threads) is not int or threads < 1:
        raise ValueError('Search threads must be a positive integer')
    manifest, db_hash = verify_database(database, expected)
    spec = manifest['specification']
    refs = {r['id']: r for r in spec['proteins']}
    if len(refs) != len(spec['proteins']):
        raise ValueError('Duplicate reference identifiers')
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    raw = Path(proteins).read_bytes()
    input_hash = hashlib.sha256(raw).hexdigest()
    records = read_fasta(raw.decode('utf-8-sig'))
    queries, valid = {}, {}
    for number, (pid, description, sequence) in enumerate(records):
        sequence = sequence.removesuffix('*')
        supported = bool(sequence) and not set(sequence) - set('ACDEFGHIKLMNPQRSTVWY')
        queries[pid] = {'status': 'no_qualifying_hit' if supported else 'skipped_unsupported_residues', 'hits': []}
        if supported:
            valid['q' + str(number)] = (pid, sequence)
    version = subprocess.check_output(['diamond', 'version'], text=True).strip()
    with tempfile.TemporaryDirectory(prefix='micpscope-search-') as temporary:
        query = Path(temporary) / 'queries.faa'
        hits = Path(temporary) / 'hits.tsv'
        query.write_text(''.join('>' + key + '\n' + value[1] + '\n' for key, value in valid.items()))
        command = ['diamond', 'blastp', '--db', str(Path(database).resolve() / 'proteins.dmnd'),
                   '--query', str(query), '--out', str(hits), '--threads', str(threads),
                   '--very-sensitive', '--max-target-seqs', '0', '--max-hsps', '1',
                   '--evalue', str(evalue), '--outfmt', '6', 'qseqid', 'sseqid', 'pident',
                   'qstart', 'qend', 'sstart', 'send', 'qlen', 'slen', 'evalue', 'bitscore']
        if valid:
            subprocess.run(command, check=True)
            for line in hits.read_text().splitlines():
                q, subject, ident, qa, qb, sa, sb, qlen, slen, ev, bits = line.split('\t')
                if q not in valid or subject not in refs:
                    raise ValueError('Search returned an unknown query/reference identifier')
                qcov = 100 * (abs(int(qb) - int(qa)) + 1) / int(qlen)
                scov = 100 * (abs(int(sb) - int(sa)) + 1) / int(slen)
                ident, ev, bits = float(ident), float(ev), float(bits)
                if not all(math.isfinite(x) for x in (ident, ev, bits, qcov, scov)):
                    raise ValueError('Non-finite search score')
                if ident >= identity and min(qcov, scov) >= coverage and ev <= evalue:
                    queries[valid[q][0]]['hits'].append(dict(reference_id=subject, identity=ident,
                        query_coverage=qcov, reference_coverage=scov, evalue=ev, bitscore=bits,
                        query_start=int(qa), query_end=int(qb), reference_start=int(sa), reference_end=int(sb)))
        raw_hits = hits.read_text() if valid else ''
    for row in queries.values():
        row['hits'].sort(key=lambda h: (-h['bitscore'], h['evalue'], h['reference_id']))
        if row['hits']:
            best = row['hits'][0]['bitscore']
            near = [refs[h['reference_id']] for h in row['hits'] if h['bitscore'] >= best * .95]
            labels = {(r.get('gene', ''), r.get('protein_name', r['id'])) for r in near}
            row['status'] = 'ambiguous_candidates' if len(labels) > 1 else 'reference_candidate'
    verify_database(database, db_hash)
    if digest(proteins) != input_hash:
        raise ValueError('Protein input changed during reference search')
    annotations = spec.get('annotations', {})
    literature = annotations.get('literature_records', []) if isinstance(annotations, dict) else []
    result = dict(schema_version='1.0', input_sha256=input_hash, queries=queries,
        database=dict(name=spec['name'], version=spec['version'], checksums_sha256=db_hash,
                      manifest_sha256=digest(Path(database) / 'reference-manifest.json')),
        tool_version=version, command=command, thresholds=dict(identity=identity, query_coverage=coverage,
        reference_coverage=coverage, evalue=evalue, ambiguity_bitscore_fraction=.95),
        references=refs, literature_records=literature,
        bibliography=spec.get('references', []),
        interpretation='Similarity candidates only. Search filters are not validated functional cutoffs. No qualifying hit does not establish biological absence.')
    output.with_suffix('.hits.tsv').write_text(raw_hits)
    output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    return result


def attach(proteins, path, input_sha256):
    result = json.loads(Path(path).read_text())
    if result.get('schema_version') != '1.0' or result.get('input_sha256') != input_sha256:
        raise ValueError('Reference results do not match this protein input')
    if set(result['queries']) != {p['protein_id'] for p in proteins}:
        raise ValueError('Reference query identifiers do not match this input')
    for protein in proteins:
        row = result['queries'][protein['protein_id']]
        protein['reference_search_status'] = row['status']
        protein['reference_hits'] = row['hits']
        if row['hits']:
            hit = row['hits'][0]
            protein['best_reference'] = hit['reference_id']
            for key in ('identity', 'query_coverage', 'reference_coverage', 'evalue', 'bitscore'):
                protein['reference_' + key] = hit[key]
        protein['reference_database_version'] = result['database']['version']
    return result
