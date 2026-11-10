"""Validated, staged reference builds. Does not infer biological claims."""
import sys
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1] / "src"))
from micpscope.application import application, archive_name, catalog
from micpscope.retention import retained_attempt, publish
import argparse
import hashlib
import json
import re
import shutil
from pathlib import Path
import zipfile
import yaml
from downloads import download, check_sequence, validate_record


def sha(data):
    return hashlib.sha256(data).hexdigest()


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')


def fasta(path):
    result = {}; key = None
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            key = line[1:].split()[0]
            if not re.fullmatch(r'[A-Za-z0-9_.:|+-]+', key) or key in result:
                raise ValueError('Invalid or duplicate FASTA identifier')
            result[key] = ''
        elif line.strip():
            if key is None:
                raise ValueError('Sequence before FASTA header')
            result[key] += line.strip().upper()
    if not result or any(not seq for seq in result.values()):
        raise ValueError('Empty FASTA sequence/file')
    return result


def _validate(spec_path, output):
    spec_path = Path(spec_path).resolve()
    s = yaml.safe_load(spec_path.read_text())
    if not isinstance(s, dict):
        raise ValueError('Reference specification must be a mapping')
    application(s)
    catalog(s)
    if set(s) - {'schema_version','application','name','version','proteins','regulatory','annotations','modules','references'}:
        raise ValueError('Unknown specification fields')
    rows = s.get('proteins', [])
    if not isinstance(rows, list) or not rows:
        raise ValueError('At least one protein source is required')
    ids = set(); assets = {}
    for row in rows + s.get('regulatory', []):
        if not isinstance(row,dict) or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,79}',row.get('id','')) or row['id'].casefold() in ids:
            raise ValueError('Unique safe record IDs required')
        ids.add(row['id'].casefold())
        if not row.get('references') or not row.get('evidence') or not row.get('license'):
            raise ValueError('Each sequence requires references, evidence and license')
        if 'path' in row:
            p = Path(row['path']).expanduser()
            p = (spec_path.parent / p).resolve() if not p.is_absolute() else p.resolve()
            if not p.is_file():
                raise ValueError('Missing input file: '+str(p))
            row['path'] = str(p); assets[str(p)] = sha(p.read_bytes())
            seqs = fasta(p)
            if row.get('record') not in seqs:
                raise ValueError('Select an exact FASTA record for '+row['id'])
        elif row in rows:
            validate_record(row)
        else:
            raise ValueError('Regulatory records require a local reference FASTA')
    for row in s.get('regulatory', []):
        for field in ('accession_version','feature_type','start','end','strand','sequence_sha256'):
            if not row.get(field): raise ValueError('Regulatory record missing '+field)
        if row['feature_type'] not in {'promoter','binding_site','upstream_region','regulatory_rna'}:
            raise ValueError('Unknown regulatory feature type')
        seq=fasta(row['path'])[row['record']]
        a,b=row['start'],row['end']
        if type(a) is not int or type(b) is not int or not 1<=a<=b<=len(seq) or row['strand'] not in {'+','-'}:
            raise ValueError('Expected 1-based inclusive interval and explicit strand')
        if set(seq)-set('ACGTN'): raise ValueError('Reference DNA must contain ACGTN only')
    for name in ('annotations','modules','references'):
        if name in s and not isinstance(s[name],(dict,list)):
            raise ValueError(name+' must be a mapping or list')
    write(output, {'spec':s,'input_sha256':assets,'spec_sha256':sha(spec_path.read_bytes()),'coordinate_system':'1-based-inclusive'})


def _retrieve(plan_path, outdir):
    plan=json.loads(Path(plan_path).read_text()); s=plan['spec']; out=Path(outdir); out.mkdir(parents=True,exist_ok=True)
    for path,digest in plan['input_sha256'].items():
        if sha(Path(path).read_bytes()) != digest: raise ValueError('Input changed after validation: '+path)
    proteins=[]; dna=[]; provenance=[]
    for row in s['proteins']:
        if 'path' in row:
            seq=fasta(row['path'])[row['record']]
            if not re.fullmatch('[ACDEFGHIKLMNPQRSTVWYUBZOX]+',seq): raise ValueError('Invalid protein alphabet')
            if row.get('sequence_sha256') and sha(seq.encode())!=row['sequence_sha256']: raise ValueError('Protein checksum mismatch')
        else:
            raw=download(row['accession'])
            (out/(row['id']+'.source.json')).write_bytes(raw)
            seq=check_sequence(row,raw)
        proteins.append('>'+row['id']+'\n'+seq+'\n')
        provenance.append({'id':row['id'],'sequence_sha256':sha(seq.encode()),'length':len(seq)})
    for row in s.get('regulatory',[]):
        seq=fasta(row['path'])[row['record']][row['start']-1:row['end']]
        if row['strand']=='-': seq=seq.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
        if sha(seq.encode()) != row['sequence_sha256']: raise ValueError('Regulatory sequence checksum mismatch')
        dna.append('>'+row['id']+'\n'+seq+'\n')
        provenance.append({'id':row['id'],'sequence_sha256':sha(seq.encode()),'length':len(seq)})
    (out/'proteins.faa').write_text(''.join(proteins)); (out/'regulatory.fna').write_text(''.join(dna))
    write(out/'sequence-provenance.json',provenance)


def _package(plan_path, seqdir, index, output):
    plan=json.loads(Path(plan_path).read_text()); out=Path(output); out.mkdir(parents=True,exist_ok=True)
    for p in Path(seqdir).iterdir():
        if p.is_file(): shutil.copyfile(p,out/p.name)
    shutil.copyfile(index,out/'proteins.dmnd')
    version=Path(index).with_suffix('.version.txt')
    shutil.copyfile(version,out/'diamond-version.txt')
    write(out/'reference-manifest.json',{'schema_version':'1.0','builder_version':'0.1.0','specification':plan['spec'],'input_sha256':plan['input_sha256'],'spec_sha256':plan['spec_sha256'],'coordinate_system':plan['coordinate_system'],'scientific_validation':'User-supplied evidence; build completion is not biological validation'})
    archive=archive_name(plan['spec'])
    hashes={p.name:sha(p.read_bytes()) for p in out.iterdir() if p.is_file() and p.name not in {'checksums.json','completed.json','database.zip',archive}}
    write(out/'checksums.json',hashes)
    with zipfile.ZipFile(out/archive,'w',zipfile.ZIP_DEFLATED) as z:
        for name in sorted([*hashes,'checksums.json']):
            z.write(out/name,name)
    # Compatibility alias for older Bactologix installations.
    shutil.copyfile(out/archive,out/'database.zip')


def _verify(root, receipt):
    root=Path(root); hashes=json.loads((root/'checksums.json').read_text())
    for name,digest in hashes.items():
        if Path(name).name!=name or sha((root/name).read_bytes())!=digest: raise ValueError('Database checksum mismatch')
    for name in ('proteins.faa','proteins.dmnd','reference-manifest.json','diamond-version.txt'):
        if name not in hashes or (root/name).stat().st_size==0: raise ValueError('Incomplete database')
    spec=json.loads((root/'reference-manifest.json').read_text())['specification']
    application(spec)
    archive=root/archive_name(spec)
    if not archive.is_file() and spec.get('schema_version')=='1.0': archive=root/'database.zip'
    with zipfile.ZipFile(archive) as z:
        if set(z.namelist()) != set(hashes)|{'checksums.json'}: raise ValueError('Package file set mismatch')
        for name,digest in hashes.items():
            if sha(z.read(name))!=digest: raise ValueError('Package checksum mismatch')
        if z.read('checksums.json')!=(root/'checksums.json').read_bytes(): raise ValueError('Package manifest mismatch')
    write(receipt,{'status':'completed','schema_version':'1.0','checksums_sha256':sha((root/'checksums.json').read_bytes()),'archive_sha256':sha(archive.read_bytes()),'archive_name':archive.name,'application':application(spec)})


def validate(spec_path, output):
    output = Path(output)
    with retained_attempt(output.parent, "validate") as attempt:
        (attempt / "product").mkdir()
        generated = attempt / "product" / output.name
        _validate(spec_path, generated)
        publish(generated, output)


def retrieve(plan_path, outdir):
    outdir = Path(outdir)
    with retained_attempt(outdir.parent, "retrieve") as attempt:
        generated = attempt / "sequences"
        _retrieve(plan_path, generated)
        publish(generated, outdir)


def package(plan_path, seqdir, index, output):
    output = Path(output)
    with retained_attempt(output.parent, "package") as attempt:
        generated = attempt / "database"
        _package(plan_path, seqdir, index, generated)
        publish(generated, output)


def verify(root, receipt):
    receipt = Path(receipt)
    with retained_attempt(receipt.parent, "verify") as attempt:
        (attempt / "product").mkdir()
        generated = attempt / "product" / receipt.name
        _verify(root, generated)
        publish(generated, receipt)


if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('stage',choices=['validate','retrieve','package','verify']); p.add_argument('args',nargs='+'); a=p.parse_args()
    globals()[a.stage](*a.args)
