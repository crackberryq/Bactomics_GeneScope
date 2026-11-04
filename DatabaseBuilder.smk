"""GeneScope Database Builder: sample-free, user-defined reference builds."""
from pathlib import Path
import yaml
import sys
ROOT = Path(workflow.basedir)
SPEC = str(Path(config.get('reference_spec', ROOT / 'database_builder/references/micp-starter.yaml')).resolve())
OUT = str(Path(config['build_dir']).resolve())
sys.path.insert(0, str(ROOT / 'src'))
from micpscope.application import archive_name
SCRIPT = str(ROOT / 'database_builder/build.py')
CODE = [SCRIPT, str(ROOT / 'database_builder/downloads.py'), str(ROOT / 'src/micpscope/application.py')]
DOCUMENT = yaml.safe_load(Path(SPEC).read_text())
ASSETS = sorted({str((Path(SPEC).parent / row['path']).resolve()) for group in ('proteins', 'regulatory') for row in DOCUMENT.get(group, []) if 'path' in row})

rule all:
    input: OUT + '/database/completed.json'

rule validate_references:
    input: spec=SPEC, code=CODE, assets=ASSETS
    output: OUT + '/validated.json'
    log: OUT + '/logs/validate.log'
    benchmark: OUT + '/benchmarks/validate.tsv'
    conda: 'envs/database-builder.yml'
    shell: 'python {SCRIPT:q} validate {input.spec:q} {output:q} > {log:q} 2>&1'

rule retrieve_references:
    input: plan=rules.validate_references.output, code=CODE
    output: directory(OUT + '/sequences')
    log: OUT + '/logs/retrieve.log'
    benchmark: OUT + '/benchmarks/retrieve.tsv'
    conda: 'envs/database-builder.yml'
    shell: 'python {SCRIPT:q} retrieve {input.plan:q} {output:q} > {log:q} 2>&1'

rule index_proteins:
    input: rules.retrieve_references.output
    output: index=OUT + '/index/proteins.dmnd', version=OUT + '/index/proteins.version.txt'
    log: OUT + '/logs/index.log'
    benchmark: OUT + '/benchmarks/index.tsv'
    threads: config.get('threads', 2)
    resources: mem_mb=config.get('memory_mb', 2048)
    conda: 'envs/database-builder.yml'
    shell: 'diamond version > {output.version:q}; diamond makedb --in {input:q}/proteins.faa --db {output.index:q} --threads {threads} > {log:q} 2>&1'

rule package_database:
    input: plan=rules.validate_references.output, sequences=rules.retrieve_references.output, index=rules.index_proteins.output.index, version=rules.index_proteins.output.version, code=CODE
    output: checksums=OUT + '/database/checksums.json', archive=OUT + '/database/' + archive_name(DOCUMENT), legacy=OUT + '/database/database.zip'
    log: OUT + '/logs/package.log'
    benchmark: OUT + '/benchmarks/package.tsv'
    conda: 'envs/database-builder.yml'
    params: destination=OUT + '/database'
    shell: 'python {SCRIPT:q} package {input.plan:q} {input.sequences:q} {input.index:q} {params.destination:q} > {log:q} 2>&1'

rule verify_database:
    input: rules.package_database.output, code=CODE
    output: OUT + '/database/completed.json'
    log: OUT + '/logs/verify.log'
    benchmark: OUT + '/benchmarks/verify.tsv'
    conda: 'envs/database-builder.yml'
    params: destination=OUT + '/database'
    shell: 'python {SCRIPT:q} verify {params.destination:q} {output:q} > {log:q} 2>&1'
