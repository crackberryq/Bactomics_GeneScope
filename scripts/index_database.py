"""Build a reference index in a persistent attempt, preserving failure products."""
from pathlib import Path
import subprocess
import sys

root, source, destination, version_destination, threads = sys.argv[1:]
sys.path.insert(0, str(Path(root) / "src"))
from micpscope.retention import retained_attempt, publish, write_receipt

with retained_attempt(Path(destination).parent, "index") as attempt:
    version = attempt / "diamond-version.txt"
    with version.open("wb") as stdout, (attempt / "version.stderr.log").open("wb") as stderr:
        subprocess.run(["diamond", "version"], check=True, stdout=stdout, stderr=stderr, cwd=attempt)
    index = attempt / "proteins.dmnd"
    command = ["diamond", "makedb", "--in", str(Path(source).resolve()), "--db", str(index.resolve()), "--threads", threads,
               "--tmpdir", str(attempt), "--no-unlink"]
    write_receipt(attempt / "command.json", {"argv": command})
    with (attempt / "stdout.log").open("wb") as stdout, (attempt / "stderr.log").open("wb") as stderr:
        subprocess.run(command, check=True, stdout=stdout, stderr=stderr, cwd=attempt)
    publish(index, destination)
    publish(version, version_destination)
