"""GeneScope: protein descriptors followed by evidence reporting."""
import csv
from pathlib import Path
import re
import sys

configfile: "config.yaml"

ROOT = Path(workflow.basedir)
SHEET = Path(config["samples"]).resolve()
OUT = str(Path(config["output_dir"]).resolve())
sys.path.insert(0, str(ROOT / "src"))
from micpscope.artifact_retention import install
install(workflow, OUT)

DATABASE = config.get("reference_database", "")
MEMORY = config.get("memory_mb", 2048)
if isinstance(MEMORY, bool) or not isinstance(MEMORY, int) or MEMORY < 512:
    raise ValueError("memory_mb must be an integer of at least 512")
SAMPLES = {}
with SHEET.open(newline="", encoding="utf-8-sig") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    if reader.fieldnames != ["sample_id", "proteins", "evidence"]:
        raise ValueError("Sample sheet columns must be sample_id, proteins, evidence")
    for row in reader:
        if None in row or any(value is None for value in row.values()):
            raise ValueError("Malformed sample sheet row")
        sid = row["sample_id"]
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,99}", sid):
            raise ValueError("Invalid sample ID: " + sid)
        if sid.casefold() in {s.casefold() for s in SAMPLES}:
            raise ValueError("Duplicate sample ID: " + sid)
        if not row["proteins"].strip():
            raise ValueError("Protein FASTA is required: " + sid)
        SAMPLES[sid] = row
if not SAMPLES:
    raise ValueError("Sample sheet contains no samples")

rule all:
    input:
        expand(OUT + "/{sample}/micpscope/report.html", sample=SAMPLES)

rule protein_properties:
    input:
        proteins=lambda w: SAMPLES[w.sample]["proteins"],
        evidence=lambda w: [SAMPLES[w.sample]["evidence"]] if SAMPLES[w.sample]["evidence"] else [],
        sidecars=lambda w: [str(p) for suffix in (".gff", ".interpro.tsv", ".fna") if (p := Path(SAMPLES[w.sample]["proteins"]).with_suffix(suffix)).is_file()],
        reference=lambda w: [OUT + "/" + w.sample + "/micpscope/reference-search.json"] if DATABASE else [],
        sheet=str(SHEET),
        code=[str(ROOT / "src/micpscope" / name) for name in ["application.py", "analysis.py", "retention.py", "artifact_retention.py", "reference.py", "cli.py", "knowledge.py", "context.py", "data/micp_catalog.json", "__init__.py"]] + [str(ROOT / "scripts/properties.py")]
    output:
        data=OUT + "/{sample}/micpscope/results.json",
        table=OUT + "/{sample}/micpscope/proteins.tsv"
    log:
        OUT + "/{sample}/micpscope/logs/properties.log"
    benchmark:
        OUT + "/{sample}/micpscope/benchmarks/properties.tsv"
    threads: 1
    resources:
        mem_mb=MEMORY
    conda: "envs/properties.yml"
    params:
        script=str(ROOT / "scripts/properties.py"), root=str(ROOT),
        reference=lambda w: "--reference-results " + __import__("shlex").quote(OUT + "/" + w.sample + "/micpscope/reference-search.json") if DATABASE else ""
    shell:
        "python {params.script:q} {params.root:q} {wildcards.sample:q} {input.proteins:q} {output.data:q} {output.table:q} {input.evidence:q} {params.reference} > {log:q} 2>&1"

rule report:
    input:
        data=rules.protein_properties.output.data,
        table=rules.protein_properties.output.table,
        code=[str(ROOT / "src/micpscope" / name) for name in ["cli.py", "retention.py", "artifact_retention.py", "report.py", "context.py", "reference.py", "application.py", "explorer.py", "report_assets.py"]] + [str(ROOT / "scripts/report.py")]
    output:
        OUT + "/{sample}/micpscope/report.html"
    log:
        OUT + "/{sample}/micpscope/logs/report.log"
    benchmark:
        OUT + "/{sample}/micpscope/benchmarks/report.tsv"
    threads: 1
    resources:
        mem_mb=MEMORY
    conda: "envs/properties.yml"
    params:
        script=str(ROOT / "scripts/report.py"), root=str(ROOT)
    shell:
        "python {params.script:q} {params.root:q} {input.data:q} {output:q} > {log:q} 2>&1"


if DATABASE:
    rule reference_search:
        input:
            proteins=lambda w: SAMPLES[w.sample]["proteins"],
            database=[str(Path(DATABASE) / name) for name in ("proteins.dmnd", "proteins.faa", "reference-manifest.json", "checksums.json", "completed.json", "diamond-version.txt")],
            code=[str(ROOT / "src/micpscope/retention.py"), str(ROOT / "src/micpscope/artifact_retention.py"), str(ROOT / "src/micpscope/application.py"), str(ROOT / "src/micpscope/reference.py"), str(ROOT / "src/micpscope/analysis.py"), str(ROOT / "scripts/reference_search.py")]
        output:
            data=OUT + "/{sample}/micpscope/reference-search.json",
            hits=OUT + "/{sample}/micpscope/reference-search.hits.tsv"
        log: OUT + "/{sample}/micpscope/logs/reference-search.log"
        benchmark: OUT + "/{sample}/micpscope/benchmarks/reference-search.tsv"
        threads: config.get("threads", 1)
        resources: mem_mb=MEMORY
        conda: "envs/reference-search.yml"
        params:
            script=str(ROOT / "scripts/reference_search.py"), root=str(ROOT), database=str(Path(DATABASE).resolve()),
            expected=("--expected " + __import__("shlex").quote(str(config["reference_database_sha256"]))) if config.get("reference_database_sha256") else "",
            identity=config.get("reference_min_identity", 30), coverage=config.get("reference_min_coverage", 70),
            evalue=config.get("reference_max_evalue", 1e-5)
        shell:
            "python {params.script:q} {params.root:q} {input.proteins:q} {params.database:q} {output.data:q} --threads {threads} {params.expected} --identity {params.identity:q} --coverage {params.coverage:q} --evalue {params.evalue:q} > {log:q} 2>&1"
