"""CI integration check against the bundled database, with real DIAMOND."""
import json
from pathlib import Path
import subprocess
import sys

root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'src'))
from micpscope.analysis import read_fasta
from micpscope.reference import verify_database

database=root/'starter-output/database'
manifest,_=verify_database(database)
first=read_fasta((database/'proteins.faa').read_text())[0]
work=root/'reference-smoke';work.mkdir(exist_ok=True)
query=work/'sample.faa';query.write_text('>query_with_an_unrelated_annotation\n'+first[2]+'\n')
sheet=work/'samples.tsv';sheet.write_text('sample_id\tproteins\tevidence\nsmoke\t'+str(query)+'\t\n')
cmd=[sys.executable,'-m','snakemake','-s',str(root/'MICPscope.smk'),'--cores','1','--config',
     'reference_database='+str(database),'samples='+str(sheet),'output_dir='+str(work/'results')]
subprocess.run(cmd,check=True,cwd=root)
data=json.loads((work/'results/smoke/micpscope/results.json').read_text())
hits=data['proteins'][0]['reference_hits']
assert any(h['reference_id']==first[0] and h['identity']==100 for h in hits),hits
assert data['reference_search']['database']['version']==manifest['specification']['version']
assert 'Reference search · 1 proteins with candidates' in (work/'results/smoke/micpscope/report.html').read_text()
subprocess.run(cmd+['--dry-run'],check=True,cwd=root)
