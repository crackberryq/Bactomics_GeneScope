"""Versioned reference packager. Python 3.11+, standard library only."""
import argparse
import hashlib
import json
import re
import subprocess
import tempfile
import time
import urllib.error
import urllib.request
import zipfile
from pathlib import Path

VERSION = '0.1.0'

def digest(data):
    return hashlib.sha256(data).hexdigest()

def encoded(obj):
    return (json.dumps(obj, indent=2, sort_keys=True, ensure_ascii=False) + '\n').encode()

def validate_record(row):
    for key in ('id', 'accession', 'taxon_id', 'sequence_version', 'sequence_sha256', 'references', 'license', 'curation_status'):
        if not row.get(key):
            raise ValueError(f'Missing required field: {key}')
    if not re.fullmatch(r'[A-Za-z0-9_.-]+', row['id']):
        raise ValueError('Unsafe record identifier')
    if not re.fullmatch(r'[A-Z0-9]{6,10}', row['accession']):
        raise ValueError('Expected canonical UniProt accession')
    if not re.fullmatch(r'[a-f0-9]{64}', row['sequence_sha256']):
        raise ValueError('Expected lowercase sequence SHA256')
    if row['curation_status'] != 'approved':
        raise ValueError('Sequence records must be explicitly approved')
    if not isinstance(row['references'], list) or not all(isinstance(x, str) and x.strip() for x in row['references']):
        raise ValueError('References must be a nonempty list of source identifiers')

def download(accession):
    url = f'https://rest.uniprot.org/uniprotkb/{accession}.json'
    for attempt in range(4):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'MICP-reference-builder/' + VERSION})
            with urllib.request.urlopen(req, timeout=45) as response:
                data = response.read(10_000_001)
            if len(data) > 10_000_000:
                raise ValueError('Oversized reference response')
            return data
        except urllib.error.HTTPError as exc:
            if exc.code not in (429, 500, 502, 503, 504) or attempt == 3:
                raise
        except (urllib.error.URLError, TimeoutError):
            if attempt == 3:
                raise
        time.sleep(2 ** attempt)

def check_sequence(row, data):
    obj = json.loads(data)
    if obj['primaryAccession'] != row['accession']:
        raise ValueError('Accession changed or was merged; curate again')
    if obj['organism']['taxonId'] != row['taxon_id']:
        raise ValueError('Taxon mismatch')
    if obj['entryAudit']['sequenceVersion'] != row['sequence_version']:
        raise ValueError('Sequence version mismatch')
    seq = obj['sequence']['value']
    if not seq or not re.fullmatch(r'[ACDEFGHIKLMNPQRSTVWYUBZOX]+', seq):
        raise ValueError('Invalid protein sequence')
    if obj['sequence']['length'] != len(seq) or digest(seq.encode()) != row['sequence_sha256']:
        raise ValueError('Sequence length/hash mismatch')
    return seq

def build(manifest_path, output, cache, offline=False, diamond=False):
    manifest_path, output, cache = map(Path, (manifest_path, output, cache))
    manifest = json.loads(manifest_path.read_text())
    if manifest.get('schema_version') != '1.0':
        raise ValueError('Unsupported manifest schema')
    if not re.fullmatch(r'[A-Za-z0-9_.-]+', manifest.get('release', '')):
        raise ValueError('Invalid release name')
    if output.exists():
        raise ValueError('Output already exists; use a new release filename')
    rows = manifest['proteins']
    if not rows:
        raise ValueError('No approved reference sequences; this is a draft catalog')
    ids = set()
    for row in rows:
        validate_record(row)
        if row['id'].casefold() in ids:
            raise ValueError('Duplicate reference ID')
        ids.add(row['id'].casefold())
    cache.mkdir(parents=True, exist_ok=True)
    files, fasta, provenance = {}, [], []
    for row in sorted(rows, key=lambda r: r['id']):
        key = row['accession'] + '-' + row['sequence_sha256'] + '.json'
        saved = cache / key
        if saved.exists():
            data = saved.read_bytes()
        elif offline:
            raise ValueError('Missing cached accession: ' + row['accession'])
        else:
            data = download(row['accession'])
        seq = check_sequence(row, data)
        if not saved.exists():
            saved.write_bytes(data)
        fasta.append('>' + row['id'] + '\n' + '\n'.join(seq[i:i+80] for i in range(0, len(seq), 80)) + '\n')
        files['sources/' + key] = data
        provenance.append({'id': row['id'], 'source_url': f'https://rest.uniprot.org/uniprotkb/{row["accession"]}.json', 'response_sha256': digest(data)})
    files['proteins.faa'] = ''.join(fasta).encode()
    files['manifest.json'] = encoded(manifest)
    files['provenance.json'] = encoded({'builder_version': VERSION, 'records': provenance})
    # Evidence stays verbatim and explicitly draft; never infer accessions from names.
    for name in ('evidence_catalog.json',):
        source = manifest_path.parent / name
        if source.exists():
            json.loads(source.read_text())
            files[name] = source.read_bytes()
    if diamond:
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            (root / 'proteins.faa').write_bytes(files['proteins.faa'])
            version = subprocess.check_output(['diamond', 'version'], text=True).strip()
            subprocess.run(['diamond', 'makedb', '--in', str(root / 'proteins.faa'), '--db', str(root / 'proteins')], check=True)
            files['proteins.dmnd'] = (root / 'proteins.dmnd').read_bytes()
            files['index.json'] = encoded({'tool': version, 'input_sha256': digest(files['proteins.faa'])})
    files['SHA256SUMS'] = ''.join(f'{digest(data)}  {name}\n' for name, data in sorted(files.items())).encode()
    output.parent.mkdir(parents=True, exist_ok=True)
    # Exclusive creation prevents silently replacing an existing database release.
    with output.open('xb') as handle:
        with zipfile.ZipFile(handle, 'w', compression=zipfile.ZIP_DEFLATED) as archive:
            for name, data in sorted(files.items()):
                info = zipfile.ZipInfo(manifest['release'] + '/' + name, date_time=(1980, 1, 1, 0, 0, 0))
                info.compress_type = zipfile.ZIP_DEFLATED
                archive.writestr(info, data)
    return len(rows)

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--cache', type=Path, default=Path('reference-cache'))
    p.add_argument('--offline', action='store_true')
    p.add_argument('--diamond', action='store_true')
    a = p.parse_args()
    print(f'Packaged {build(a.manifest, a.output, a.cache, a.offline, a.diamond)} references')
