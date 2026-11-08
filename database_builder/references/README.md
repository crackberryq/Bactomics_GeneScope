# MICP literature starter 1.0.1

`micp-starter.yaml` is the default Database Builder input. `starter-proteins.faa`
is shipped beside it, so sequence import works offline after the environment
has been installed. Bactologix loads the installed YAML automatically and
resolves the bundled FASTA path before snapshotting a build. Edits affect the
build specification, not this shipped file.

The starter contains 164 reviewed UniProt protein references: 157 from
Bacillus subtilis 168 and seven urease structural/accessory references from
Sporosarcina pasteurii. Accession, sequence version, entry version, taxonomy,
sequence checksum, source URL, gene aliases and the retrieval date are recorded.
Sequences are redistributed under UniProt CC BY 4.0; attribution and license
URL are in the specification. References were retrieved on 3 November 2026.

All 207 table entries from the two user-supplied literature reviews are
preserved with their original field headings, evidence boundaries and source
document hashes. These include biological candidates, definitions and methods;
207 is not a gene count. Bibliography IDs are namespaced by source document.
Protein selections point to the relevant review table entries. Functional
modules group references descriptively; they do not assert complete operons.

The sequence set is a starter, not exhaustive coverage of all organisms,
pathways or genes discussed in the reviews. B. subtilis 168 sequences represent
annotation references; experiments in other strains are not assigned to 168.
Unresolved selections are listed explicitly, including ambiguous or unavailable
reviewed matches. Historical gene aliases are not sufficient by themselves to
identify a protein; CapA/PgsA, CheA and DacB ambiguities were resolved by protein
identity. Literature claims were preserved, not independently revalidated.

The regulatory review includes relative positions, binding-site observations,
consensus motifs and inferred relationships. These are retained in the
literature evidence layer. No absolute genome intervals were established from
these documents alone, so the regulatory sequence list is empty. A build must
not present these observations as mapped sample promoters or measured activity.
GeneScope can search a completed build when it is explicitly selected for a run.
Version 1.0.1 aligns the annotation catalog and reference category vocabulary;
reference sequences are unchanged. Rebuild older packages from the corrected
specification to obtain the new version and checksums. See
[the integrity notes](../../docs/prebeta-integrity.md) for grouping and review states.
