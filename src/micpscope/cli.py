import argparse
import csv
import hashlib
import io
import json
import platform
from datetime import datetime, timezone
from pathlib import Path

import Bio
from . import __version__
from .analysis import characterize, read_fasta
from .knowledge import load_catalog, enrich
from .context import add_context
from .reference import attach
from .report import render_report, table

LIMITS = ("Sequence descriptors do not establish temperature or pH optima, shelf life, "
          "mineralization activity, material performance, or deployment safety. "
          "Protein completeness is unassessed. Header-derived functional associations are not sequence-confirmed. "
          "Evidence entries are user-supplied and are not independently verified.")
EVIDENCE_FIELDS = ["sample_id", "subject_id", "evidence_type", "property", "value", "unit", "source", "notes"]


def read_evidence(path, sample):
    if path is None:
        return []
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != EVIDENCE_FIELDS:
            raise ValueError("Evidence CSV columns must be: " + ",".join(EVIDENCE_FIELDS))
        rows = list(reader)
    for number, row in enumerate(rows, 2):
        if None in row or any(value is None for value in row.values()):
            raise ValueError(f"Malformed evidence row {number}")
        if row["sample_id"] != sample:
            raise ValueError(f"Evidence row {number} belongs to a different sample")
        if row["evidence_type"] not in {"measured", "predicted", "literature-derived", "calculated"}:
            raise ValueError(f"Invalid evidence type at row {number}")
        if any(not row[key].strip() for key in ("subject_id", "property", "value", "source")):
            raise ValueError(f"Missing subject, property, value or source at row {number}")
    return rows


def render(result):
    return render_report(result)


def main(argv=None):
    parser = argparse.ArgumentParser(description="MICP Scope protein characterization and evidence report")
    parser.add_argument("--version", action="version", version=__version__)
    parser.add_argument("--proteins", type=Path, required=True, help="Annotated protein FASTA (.faa), not genome DNA")
    parser.add_argument("--data-only", action="store_true", help="Write analysis data for the Snakemake report stage")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--evidence", type=Path, help="Optional evidence CSV")
    parser.add_argument('--gff', type=Path, help='Optional recorded GFF3 CDS coordinates')
    parser.add_argument('--interpro', type=Path, help='Optional InterProScan TSV signatures')
    parser.add_argument("--reference-results", type=Path, help="Verified DIAMOND search JSON for this exact protein input")
    parser.add_argument("--out", type=Path, required=True, help="New output directory; existing paths are never overwritten")
    args = parser.parse_args(argv)
    try:
        raw = args.proteins.read_bytes()
        proteins = [characterize(*record) for record in read_fasta(raw.decode("utf-8-sig"))]
        catalog, catalog_sha256 = load_catalog()
        enrich(proteins, catalog)
        gff = args.gff or (args.proteins.with_suffix('.gff') if args.proteins.with_suffix('.gff').is_file() else None)
        interpro = args.interpro or (args.proteins.with_suffix('.interpro.tsv') if args.proteins.with_suffix('.interpro.tsv').is_file() else None)
        context = add_context(proteins, gff, interpro)
        evidence = read_evidence(args.evidence, args.sample)
        result = dict(schema_version="1.3", software_version=__version__, sample_id=args.sample,
                      context=context,
                      proteins=proteins, evidence=evidence, operating_conditions="undetermined",
                      knowledge_catalog=catalog,
                      limitations=LIMITS, provenance=dict(
                          input_sha256=hashlib.sha256(raw).hexdigest(),
                          catalog_sha256=catalog_sha256,
                          evidence_sha256=hashlib.sha256(args.evidence.read_bytes()).hexdigest() if args.evidence else None,
                          created_utc=datetime.now(timezone.utc).isoformat(),
                          python_version=platform.python_version(), biopython_version=Bio.__version__))
        if args.reference_results:
            result["reference_search"] = attach(proteins, args.reference_results, result["provenance"]["input_sha256"])
        assembly = args.proteins.with_suffix(".fna")
        if assembly.is_file():
            result["provenance"]["assembly_sha256"] = hashlib.sha256(assembly.read_bytes()).hexdigest()
        report = None if args.data_only else render(result)
        columns = ["protein_id", "gene_name", "protein_name", "function", "category", "annotation_status", "annotation_basis", "literature_scope", "description", "length", "status", "aliphatic_index", "theoretical_pi", "gravy", "instability_index", "molecular_weight_da", "aromaticity", "sequence_sha256", "evidence_type", "completeness", "instability_class", "hydropathy_class", "temperature_indication", "ph_indication", "references", "notes"]
        columns += ["reference_search_status", "best_reference", "reference_identity", "reference_query_coverage", "reference_reference_coverage", "reference_evalue", "reference_bitscore", "reference_database_version"]
        buffer = io.StringIO(newline="")
        writer = csv.DictWriter(buffer, fieldnames=columns, extrasaction="ignore", delimiter="\t")
        writer.writeheader()
        # Prevent spreadsheet formula interpretation in exported text fields.
        writer.writerows({key: ("'" + value if isinstance(value, str) and value.startswith(("=", "+", "-", "@")) else value)
                         for key, value in row.items()} for row in proteins)
        args.out.mkdir(parents=True, exist_ok=False)
        (args.out / "results.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        (args.out / "proteins.tsv").write_text(buffer.getvalue(), encoding="utf-8")
        if report is not None:
            (args.out / "report.html").write_text(report, encoding="utf-8")
    except (OSError, ValueError) as error:
        parser.exit(2, f"MICP Scope: {error}\n")
    print(f"Outputs: {args.out}")


if __name__ == "__main__":
    main()

