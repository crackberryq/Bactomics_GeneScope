"""Compatibility entry point for existing MICP Scope integrations."""
include: "GeneScope.smk"

rule legacy_entry:
    input: rules.all.input
    default_target: True
