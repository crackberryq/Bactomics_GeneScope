"""Traceable protein similarity search; reference associations are not activity calls."""
from .application import application, catalog, validate_catalog, validate_reference_categories
import hashlib
import json
import math
from pathlib import Path
import subprocess
from .retention import retained_attempt, publish, write_receipt


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
    application(manifest['specification'])
    catalog(manifest['specification'])
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
    with retained_attempt(output.parent, 'reference-search') as temporary:
        with (temporary / 'version.stderr.log').open('wb') as stderr:
            try:
                version_output = subprocess.check_output(['diamond', 'version'], text=True, stderr=stderr)
            except subprocess.CalledProcessError as error:
                (temporary / 'diamond-version.txt').write_text(error.output or '')
                raise
        (temporary / 'diamond-version.txt').write_text(version_output)
        version = version_output.strip()
        query = Path(temporary) / 'queries.faa'
        hits = Path(temporary) / 'hits.tsv'
        native_tmp = temporary / 'native-tmp'
        native_tmp.mkdir()
        query.write_text(''.join('>' + key + '\n' + value[1] + '\n' for key, value in valid.items()))
        command = ['diamond', 'blastp', '--db', str(Path(database).resolve() / 'proteins.dmnd'),
                   '--query', str(query), '--out', str(hits), '--threads', str(threads),
                   '--tmpdir', str(native_tmp), '--no-unlink',
                   '--very-sensitive', '--max-target-seqs', '0', '--max-hsps', '1',
                   '--evalue', str(evalue), '--outfmt', '6', 'qseqid', 'sseqid', 'pident',
                   'qstart', 'qend', 'sstart', 'send', 'qlen', 'slen', 'evalue', 'bitscore']
        write_receipt(temporary / 'command.json', {'argv': command, 'input_sha256': input_hash, 'database_sha256': db_hash})
        if valid:
            with (temporary / 'stdout.log').open('wb') as stdout, (temporary / 'stderr.log').open('wb') as stderr:
                subprocess.run(command, check=True, stdout=stdout, stderr=stderr, cwd=temporary)
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
        if not valid:
            hits.write_text('')
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
            application=application(spec), knowledge_catalog=catalog(spec),
            database=dict(**application(spec), checksums_sha256=db_hash,
                          manifest_sha256=digest(Path(database) / 'reference-manifest.json')),
            tool_version=version, command=command, thresholds=dict(identity=identity, query_coverage=coverage,
            reference_coverage=coverage, evalue=evalue, ambiguity_bitscore_fraction=.95),
            references=refs, literature_records=literature,
            bibliography=spec.get('references', []),
            interpretation='Similarity candidates only. Search filters are not validated functional cutoffs. No qualifying hit does not establish biological absence.')
        generated = temporary / 'reference-search.json'
        generated.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
        publish(hits, output.with_suffix('.hits.tsv'))
        publish(generated, output)
        return result


def attach(proteins, path, input_sha256):
    result = json.loads(Path(path).read_text())
    if result.get('schema_version') != '1.0' or result.get('input_sha256') != input_sha256:
        raise ValueError('Reference results do not match this protein input')
    if set(result['queries']) != {p['protein_id'] for p in proteins}:
        raise ValueError('Reference query identifiers do not match this input')
    knowledge = validate_catalog(result.get('knowledge_catalog', {}))
    validate_reference_categories(result['references'].values(), knowledge)
    for protein in proteins:
        row = result['queries'][protein['protein_id']]
        protein['reference_search_status'] = row['status']
        protein['reference_hits'] = row['hits']
        if row['hits']:
            hit = row['hits'][0]
            protein['best_reference'] = hit['reference_id']
            for key in ('identity', 'query_coverage', 'reference_coverage', 'evalue', 'bitscore'):
                protein['reference_' + key] = hit[key]
        apply_category_evidence(protein, result)
        protein['reference_database_version'] = result['database']['version']
    return result


def apply_category_evidence(protein, result):
    """Keep supplied labels and reference grouping separate, including saved reports."""
    row = result.get('queries', {}).get(protein['protein_id'], {})
    hits = row.get('hits', protein.get('reference_hits', []))
    status = row.get('status', protein.get('reference_search_status', 'unassessed'))
    supplied = protein.get('supplied_annotation')
    if supplied is None:
        supplied = {k: protein.get(k) for k in ('gene_name', 'protein_name', 'function', 'category', 'annotation_status', 'annotation_basis')}
        if protein.get('category_basis', '').startswith('Reference similarity'):
            supplied['category'] = 'unassigned'
    protein['supplied_annotation'] = supplied
    protein['annotation_category'] = supplied.get('category', 'unassigned')
    protein['reference_search_status'] = status
    protein['reference_hits'] = hits
    protein['reference_categories'] = []
    protein['reference_gene'] = ''
    protein['reference_protein_name'] = ''
    protein['category'] = 'unassigned'
    protein['evidence_status'] = status
    protein['evidence_conflicts'] = []
    protein['category_basis'] = 'Unresolved or unassessed reference grouping; supplied annotation retained separately'
    if not hits:
        return
    ref = result['references'][hits[0]['reference_id']]
    protein['reference_gene'] = ref.get('gene', '')
    protein['reference_protein_name'] = ref.get('protein_name', '')
    best = hits[0]['bitscore']
    near = [result['references'][h['reference_id']] for h in hits if h['bitscore'] >= best * .95]
    categories = list(dict.fromkeys(c for r in near for c in r.get('categories', [])))
    protein['reference_categories'] = categories
    known = {c['id'] for c in result.get('knowledge_catalog', {}).get('categories', [])}
    if set(categories) - known:
        protein['evidence_status'] = 'unresolved_reference_category'
        return
    if status == 'ambiguous_candidates' or len({tuple(sorted(r.get('categories', []))) for r in near}) > 1:
        protein['evidence_status'] = 'ambiguous_candidates'
        return
    if status != 'reference_candidate':
        return
    annotation = protein['annotation_category']
    if annotation not in (None, '', 'unassigned') and annotation not in categories:
        protein['evidence_conflicts'].append('Supplied annotation category differs from the reference categories')
    supplied_genes = {s.strip().casefold() for s in str(supplied.get('gene_name', '')).split(',')} - {'', 'unassigned', 'family only'}
    for entry in result.get('knowledge_catalog', {}).get('entries', []):
        if entry['id'] in protein.get('catalog_ids', []):
            supplied_genes.update(g.casefold() for g in entry.get('gene_symbols', []))
    reference_genes = {r.get('gene', '').casefold() for r in near} - {''}
    if supplied_genes and reference_genes and not supplied_genes & reference_genes:
        protein['evidence_conflicts'].append('Supplied and reference gene labels differ; nomenclature or annotation requires review')
    if protein['evidence_conflicts']:
        protein['evidence_status'] = 'annotation_reference_conflict'
        return
    if categories:
        protein['category'] = categories[0]
        protein['category_basis'] = ('Reference similarity candidate; primary grouping uses the first category declared by the reference. '
                                     'All reference categories are retained; function is not experimentally established')


def evidence_label(protein, searched):
    if not searched:
        return 'Supplied annotation · search not requested'
    return {
        'reference_candidate': 'Reference candidate · function unconfirmed',
        'ambiguous_candidates': 'Ambiguous reference candidates · unresolved',
        'annotation_reference_conflict': 'Annotation / reference conflict · unresolved',
        'unresolved_reference_category': 'Unknown reference category · unresolved',
        'skipped_unsupported_residues': 'Reference search skipped · unsupported residues',
        'no_qualifying_hit': 'No qualifying reference hit · absence unproven',
        'not_requested': 'Reference search not requested',
    }.get(protein.get('evidence_status', protein.get('reference_search_status')), 'Reference assessment unverified')
