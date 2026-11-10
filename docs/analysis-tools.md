# Analysis tools and evidence coverage

| Component | Implementation in 0.1.0a18 | Inputs / interpretation |
|---|---|---|
| Snakemake 9.6.2 | Existing native workflow, retained | Tracks optional context inputs and code changes |
| Biopython 1.86 ProteinAnalysis | Runs locally | Descriptors and residue composition; no operating optima |
| Protein quality checks | Runs locally | Residue flags and exact sequence groups; no inferred frameshifts or truncation |
| GFF3 context reader | Runs when a matching sample GFF exists | Declared-bound checks, identifier-only linkage warnings and circular display; no sequence pairing verification or operon prediction |
| DIAMOND reference search | Runs only with a selected completed database | Similarity candidates with recorded coverage filters, competing matches and explicit annotation conflicts; no experimentally confirmed function |
| InterProScan TSV reader | Runs when a matching sample TSV exists | Imports classic 11–15-column signatures, verifies sequence MD5/length and bounds; does not execute InterProScan |
| Report renderer | Runs locally, offline | Native HTML disclosures, SVG tracks and JavaScript search/export; no CDN dependency |

## Additional analysis stages still to implement

- Operon prediction needs genomic sequence/coordinates and an explicitly validated prediction method. Gene inventory, adjacency and predicted co-transcription must remain separate outputs.
- Fresh domain annotation needs a pinned analysis engine and reference release, with Bactologix resource installation, licensing review, checksums, resource requirements and reproducibility tests. The current classic TSV reader must not be assumed compatible with every future InterProScan format.
- Localization, secretion signals, transmembrane topology and disorder need separate validated predictors and versioned models. Their outputs must retain prediction confidence and method scope.
- Structure and experimental proteomics need explicit source/identity mapping, confidence or peptide coverage, and provenance. No generic overall performance score is implemented.

These stages are not silently enabled or represented as completed analyses. Existing sample annotation files are reused only as documented above. Report assets and the literature catalog are installed through the normal pipeline update.

## Input specifications

- GFF3: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
- Classic InterProScan TSV: https://interproscan-docs.readthedocs.io/en/v5/OutputFormats.html
- Protein descriptors: https://biopython.org/docs/1.86/api/Bio.SeqUtils.ProtParam.html
