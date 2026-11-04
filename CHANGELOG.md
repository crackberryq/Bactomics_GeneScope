## GeneScope 0.1.0a14

- Align the offline report typography, navy sidebar, cards, tables and spacing with BaSIS.
- Link protein records to their genome coordinates; search the explorer by gene name or identifier and choose among multiple coordinate matches.
- Add labeled gene arrows, selected-gene highlighting, pan and zoom controls, and a collapsible feature table.
- Place proteins without qualifying database hits in a collapsed section at the end. Preserve all source annotations and evidence.
- Presentation and navigation changes only; analysis thresholds and biological interpretation are unchanged.

# Changelog

## 0.1.0a5 — compact report and recorded protein/genome context

Fully collapsible categories and ancillary sections, category filtering, search reveal, filtered TSV export, print support, residue quality flags, composition and identical-sequence groups. Add optional exact-ID GFF3 neighborhood import and sequence-verified InterProScan signature import with offline diagrams. Keep existing inputs and output paths compatible.

## 0.1.0a4 — MICP Scope

Named `micpscope-v0.1.0a4` release, versioned 20-entry literature catalog, conservative header annotation mapping, categorized offline report with search and expandable protein details, source-linked functions and descriptor interpretations. No descriptor-to-temperature/pH optimum conversion is claimed. Catalog hash and full catalog are embedded in JSON.

## 0.1.0a3

Use **MICP Scope** as the display name in reports, CLI messages and documentation. Repository, workflow filenames and command names remain compatible.

## 0.1.0a2

Native Snakemake DAG, pinned controller and rule environments, sample-sheet batch input, separate report stage, logs and benchmarks. Preserve the existing characterization API.


## 0.1.0a1

Initial protein descriptor CLI, provenance-aware JSON/TSV exports, self-contained HTML report, evidence CSV import, synthetic demo and validation tests. This is a foundation release, not a complete MICP genome pipeline.

