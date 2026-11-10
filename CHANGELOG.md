## GeneScope 0.1.0a17

Generated-data preservation is a core current-branch feature, applied to both
main and development. Search algorithms, thresholds and interpretation are unchanged.

- Preserve prior declared products, logs and benchmark records before Snakemake
  removes them on a rerun, retry or report refresh; retain incomplete outputs.
- Retain descriptor/report attempts and native DIAMOND query files, raw hit tables,
  commands, versions, standard output/error and temporary files. Use the supported
  `--no-unlink` flag and persistent temporary directories in DIAMOND 2.1.11.
- Preserve native Database Builder outputs and failed build attempts. Publishing
  replacements keeps the previous file or directory, including native sidecars.
- Refuse replacement when preservation fails or output ownership is unsafe.
  The adapter supports pinned Snakemake 9.6.2 and 9.27.1; it refuses other versions.

Current outputs remain at their existing paths. Previous workflow products are
under `output_dir/retained/` (or `build_dir/retained/` for databases); durable
sample-stage attempts are under `<sample>/micpscope/retained/`. Export the complete
output tree to include history. No existing lost files can be recovered by this
update. Cleanup requires a separate explicit user decision; completion is never
permission to delete generated evidence.

DIAMOND retention flags were checked against the pinned upstream source:
https://github.com/bbuchfink/diamond/blob/v2.1.11/src/basic/config.cpp

## GeneScope 0.1.0a16

This alpha release applies the pre-beta evidence-integrity fixes to the current
installable channel. The same corrections are retained in development for
parallel testing. It is not a declaration of analytical validation or a beta
promotion.

- Align all 164 MICP starter reference records with the annotation category vocabulary and reject unresolved category identifiers during database build/load. The corrected starter is **MICP 1.0.1**, catalog `2026.11-category.2`; reference sequences and search thresholds are unchanged.
- Preserve supplied labels separately from reference candidates. Conflicting labels/categories, ambiguous candidates and unknown categories remain unresolved, with visible counts at the top of the report and consistent JSON/TSV/HTML states. All reference categories are retained; the first declared category is the documented primary display group for an unambiguous, nonconflicting candidate.
- Check GFF features against declared region bounds and reject contradictory declarations. Missing bounds and identifier-only linkage are explicit warnings. Circular-origin features display as two intervals while counting as one source feature; multipart neighborhoods remain unresolved.
- Track renderer modules and entry scripts so report-only code changes regenerate HTML without rerunning protein analysis. Preserve explicit unassessed operon/promoter and engineering-performance states.
- Run regression, browser, real DIAMOND smoke and Snakemake dependency checks on both main and development. Browser CI uses the Playwright container's own Python to retain matching SSL libraries.

**Updating existing work:** install the new compatible source through Bactologix.
Rebuild an old MICP 1.0.0 database from the corrected bundled starter before new
searches, select the resulting 1.0.1 package for a new draft, and retain the old
run and database. Saved reports can be refreshed to show clearer unresolved
states, but refreshing does not rerun reference searches or replace the saved
reference database. No existing results should be deleted or overwritten.

See [integrity details](docs/prebeta-integrity.md) for the exact interpretation
boundaries and regression scope. Software checks do not establish sensitivity,
specificity, functional activity, material performance or deployment safety.

## GeneScope 0.1.0a15

- Replace whole-proteome category meters with clickable application-category bars and named protein comparisons for aliphatic index, pI, GRAVY and instability index.
- Move the whole-proteome scatter plot into collapsed background content at the end of the report.
- Add a complete chromosome density overview and a default 20 kb gene view. Dense intervals show all features as density rather than silently truncating the first 500 genes.
- Search saved gene names, protein descriptions and identifiers with suggestions, and navigate from category comparisons to gene locations.
- Prefer gene names, then protein names, on records and gene arrows; preserve identifiers in details. Preserve supplied GFF gene names in new result bundles.
- Hide empty operon and unknown-strand tracks. Supplied operon annotations remain distinct from gene proximity.
- Existing compatible saved results can be regenerated through Bactologix; missing labels use available protein annotations. No new reference search is required for these presentation changes.

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

