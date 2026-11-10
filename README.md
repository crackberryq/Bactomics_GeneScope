# GeneScope

Canonical repository: [crackberryq/Bactomics_GeneScope](https://github.com/crackberryq/Bactomics_GeneScope).

**0.1.0a17 — GeneScope alpha.** Application-driven protein characterization and traceable genomic evidence, managed through Bactologix alongside HybAs AWARE and BaSIS AWARE.

## Available now

- Versioned literature catalog and categorized gene/function report with searchable, expandable protein details.
- Protein FASTA input, explicit validation and per-record ambiguity flags.
- Aliphatic index, theoretical pI, GRAVY, instability index, molecular weight, aromaticity and amino-acid counts.
- Optional sample-linked evidence CSV import with source attribution.
- Self-contained HTML report, machine-readable JSON and spreadsheet-readable protein TSV.
- Input checksums, software versions and UTC timestamps.
- Optional DIAMOND similarity candidates with separate supplied annotations, explicit ambiguity/conflict states and visible unresolved-evidence counts.
- MICP starter 1.0.1 with a consistent category vocabulary; existing 1.0.0 database builds need rebuilding before new searches.
- The direct descriptor CLI refuses existing output directories. Workflow retries and refreshes preserve prior products before publishing replacements.

## Run through Bactologix

Select GeneScope in **Install & update**, install the named `genescope-v0.1.0a17` compatible snapshot. First installation prepares the runner and rule environment. Register a sample folder containing `annotation/<sample_id>.faa`; optionally include `annotation/<sample_id>.evidence.csv`. Create a GeneScope draft, validate, prepare missing environments and run. Reports appear under Results & downloads.

## Snakemake structure

- `GeneScope.smk`: native DAG, separate properties and report rules.
- `environment.yml`: Python 3.11 / Snakemake 9.6.2 runner.
- `envs/properties.yml`: rule-local Python / Biopython environment and post-deploy check.
- `config.yaml` and `examples/samples.tsv`: configuration and batch inputs.
- `scripts/`: workflow entry points; `src/micpscope/`: reusable analysis code.
- Outputs: `<output_dir>/<sample_id>/micpscope/`, including JSON, TSV, HTML, logs and benchmarks.

For standalone development, create the Conda runner from `environment.yml`, then use:

```bash
snakemake --snakefile GeneScope.smk --cores 1 --use-conda --dry-run
snakemake --snakefile GeneScope.smk --cores 1 --use-conda --conda-create-envs-only
snakemake --snakefile GeneScope.smk --cores 1 --use-conda
```

The default configuration runs synthetic fixtures only. Sample-sheet paths resolve relative to the workflow working directory; absolute paths are recommended. Supply annotated protein FASTA, not genome DNA. The previous `micpscope` Python CLI remains available for developer use. Snakemake manages resumability and reruns; the standalone CLI still refuses existing output directories.

## Evidence

Copy `examples/evidence_template.csv` and retain its column order. Add one observation per row. `sample_id` must match `--sample`. `subject_id` identifies the protein, preparation, formulation or material tested. Allowed evidence types: `measured`, `predicted`, `literature-derived`, `calculated`. Property, value, source and subject are required; use a unit where applicable and notes for measurement context. Supply the file with `--evidence your_evidence.csv`.

Entries are displayed as supplied, not independently verified or converted into recommended conditions. No numeric response models are fitted in this release.

## Interpretation limits

Aliphatic index is not an optimum temperature. Theoretical pI is not an optimum pH. Instability index is not shelf life. Sequence properties do not demonstrate mineralization, material compatibility or safe deployment. Completeness remains unassessed. Header labels are linked conservatively to the catalog and explicitly labeled as unverified annotation; this is not sequence-based function detection. Operating conditions remain undetermined.

One terminal stop is removed and documented; sequences with other unsupported residues are retained as flagged records without property calculations. Duplicate IDs, empty records and malformed FASTA fail validation. Original sequence hashes and input-file hashes provide traceability. TSV text fields beginning with formula characters are prefixed with an apostrophe for spreadsheet safety; JSON preserves original text.

## Planned stages

1. Broader primary-literature catalog coverage and sequence-confirmed functional annotation.
2. HybAs/BaSIS result adapters with explicit schema compatibility checks.
3. Measured material-performance and formulation comparison, with uncertainty and units.
4. Expand application-specific evidence views.

These are planned features, not implemented capabilities. This alpha does not identify pathways, rank strains, infer formulation conditions, design genetic changes or certify deployment.

## Methods and development

Biopython 1.86 `ProteinAnalysis` provides standard descriptors. Aliphatic index uses `100 × [A + 2.9V + 3.9(I + L)] / length` with residue counts.

- [ExPASy ProtParam methodology](https://web.expasy.org/protparam/protparam-doc.html)
- [Biopython 1.86 API](https://biopython.org/docs/1.86/api/Bio.SeqUtils.ProtParam.html)

```bash
python -m unittest discover -s tests -v
```

No license has yet been selected by the repository owner.

## Literature catalog and annotation

The bundled MICP application specification carries its annotation catalog and literature evidence. GeneScope only loads application-specific interpretation from the selected database; an analysis without a database provides generic protein descriptors. See [application packages](docs/application-packages.md).

Annotation uses supplied FASTA headers: explicit `gene=` / UniProt `GN=` symbols, standalone known symbols, or specific product phrases. For example, `>locus_1 [gene=ureC] urease subunit alpha` can link to the MICP starter ureolysis category. This does not validate function from sequence. Conflicting matches remain unresolved; missing labels are not biological absence. Plain locus-only FASTA cannot provide gene names. Generic product descriptions remain visible even without a catalog match.

The report separates direct urease context, candidate associations and indirect matrix/sporulation context. Instability index uses the documented below/above-40 heuristic, with exactly 40 left unclassified. No universal aliphatic-index temperature bins or pI activity-pH bins are supported by the reviewed methodology. Numeric operating optima remain undetermined; sample-linked measured or literature observations can be supplied through the existing evidence CSV with identity and assay context.

For later table curation, each catalog entry has `id`, `gene_symbols`, `protein_name`, `product_aliases`, `category`, `function`, `references`, `scope`, and `evidence_level`. Update cited sources and catalog version together. Avoid broad aliases that conflate unrelated families.

## Compact report and optional recorded context

Reports open with a compact overview. Unresolved-evidence counts remain visible at the top. Categories, unresolved records, quality flags, protein details, evidence and methods are collapsed. Empty categories live under a separate disclosure. Search opens matching categories and records; category filters, Collapse all, filtered TSV export and Print / Save PDF work offline.

New optional inputs are discovered beside the protein FASTA: `annotation/<sample_id>.gff` for GFF3 CDS coordinates and `annotation/<sample_id>.interpro.tsv` for classic InterProScan TSV (11–15 columns). Bactologix stages existing files automatically; no manual installation of report assets is needed. Protein-only runs remain supported. These are sample analysis inputs, not downloadable reference databases.

GFF3 joins use exact CDS `ID`, `protein_id` or `locus_tag`. Multiple matching CDS segments remain unresolved for neighborhood display. Coordinates use the supplied annotation, without rechecking the genomic sequence. Neighborhoods show two adjacent CDS records on either side, and do not establish an operon or co-transcription.

Imported InterProScan matches require the same protein identifier, sequence length and MD5 as the normalized input protein. Duplicate identical rows are merged. Signature coordinates must lie within the sequence; mismatches fail explicitly. Unmatched identifiers are counted, and a missing TSV row means assessment coverage is unknown. Signature tracks include families/sites as well as domains. No external analysis is run implicitly.

Protein cards also show amino-acid composition, residue flags and identical-sequence record groups. They do not infer frameshifts, truncation, genomic copy number or biological completeness. See [analysis-tools.md](docs/analysis-tools.md) for implemented tools and remaining analysis stages.


## Database Builder

Build user-defined reference packages through Bactologix setup using [DatabaseBuilder.smk](DatabaseBuilder.smk). See [input format and limitations](docs/database-builder.md).

## Pre-beta integrity release

Version 0.1.0a16 carries the reviewed evidence-integrity fixes in both the current release and the development channel. See [release notes](docs/release-notes.md) and [integrity details](docs/prebeta-integrity.md). It remains an alpha for controlled testing; these fixes and software checks do not establish analytical accuracy or deployment safety. Bactologix and Bactologix Dev select their own pinned sources and workspaces.

## Generated data are retained

Retention is a core feature of the current workflow. GeneScope does not clean up
generated research products after successful or failed execution. Before a
Snakemake job replaces its declared outputs, logs or benchmark files, their
previous versions are moved into `output_dir/retained/<attempt>/files/`, with a
manifest identifying their original paths. Database Builder uses the same policy
under `build_dir/retained/`. This includes retries within one invocation. Failure
to preserve a prior product stops execution; no copy/delete fallback is used.

Each descriptor, report, reference-search and database-build stage also keeps a
durable attempt directory. For sample analyses these are under
`<sample>/micpscope/retained/`. Reference-search attempts contain the exact query
FASTA, unfiltered native hit table, command, tool version, standard output/error,
and DIAMOND temporary files. DIAMOND runs with `--no-unlink` and a persistent
temporary directory. Generated JSON/TSV/HTML products remain in their attempt
directories after publication. Failed attempts keep their partial products and
an explicit failure receipt; they are not presented as successful results.

The retention adapter is checked against Snakemake 9.6.2 and 9.27.1 and fails
closed for unsupported versions or unsafe output ownership. It keeps incomplete
outputs and disables automatic temporary-output deletion. It does not preserve
unwritten in-memory state, recover files removed by older versions, or turn a
retained failed artifact into a validated result. Export the complete output tree,
including every `retained/` directory, to carry the history with the isolate.
Disk-space cleanup must be a separate, explicit user action after reviewing the
files; the workflows never infer permission to delete from a successful run.
