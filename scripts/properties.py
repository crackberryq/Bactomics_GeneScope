"""Run the existing analysis API in an isolated staging directory."""
import os
from pathlib import Path
import sys
import tempfile

root, sample, proteins, data, table, *extra = sys.argv[1:]
evidence = extra[:1] if extra and extra[0] != "--reference-results" else []
reference_args = extra[len(evidence):]
sys.path.insert(0, str(Path(root) / "src"))
from micpscope.cli import main

parent = Path(data).parent
parent.mkdir(parents=True, exist_ok=True)
with tempfile.TemporaryDirectory(prefix=".properties-", dir=parent) as staging:
    output = Path(staging) / "analysis"
    args = ["--proteins", proteins, "--sample", sample, "--out", str(output), "--data-only"]
    if evidence:
        args += ["--evidence", evidence[0]]
    main(args + reference_args)
    os.replace(output / "results.json", data)
    os.replace(output / "proteins.tsv", table)

