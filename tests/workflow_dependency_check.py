"""Execute the small descriptor fixture and verify isolated rendering code reruns."""
from pathlib import Path
import os
import re
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
with tempfile.TemporaryDirectory(prefix='genescope-render-dependencies-') as td:
    root = Path(td) / 'source'
    root.mkdir()
    for name in ('src', 'scripts', 'examples', 'envs'):
        shutil.copytree(ROOT / name, root / name)
    for name in ('GeneScope.smk', 'config.yaml'):
        shutil.copyfile(ROOT / name, root / name)
    if os.environ.get('GENESCOPE_TEST_NO_BENCHMARK') == '1':
        workflow = root / 'GeneScope.smk'
        workflow.write_text(re.sub(r'^    benchmark:\n        .*\n', '', workflow.read_text(), flags=re.M))
    command = [sys.executable, '-m', 'snakemake', '-s', 'GeneScope.smk', '--cores', '1']
    subprocess.run(command, cwd=root, check=True)
    data = root / 'results/demo/micpscope/results.json'
    table = data.with_name('proteins.tsv')
    original = [(p.read_bytes(), p.stat().st_mtime_ns) for p in (data, table)]
    for name in ('src/micpscope/explorer.py', 'src/micpscope/report_assets.py', 'scripts/report.py'):
        path = root / name
        path.write_text(path.read_text() + '\n# Rendering dependency regression fixture.\n')
        # Make the edit newer than all saved outputs even on coarse filesystem clocks.
        newest = max(p.stat().st_mtime for p in data.parent.iterdir() if p.is_file())
        os.utime(path, (newest + 1, newest + 1))
        dry = subprocess.run(command + ['--dry-run'], cwd=root, check=True, text=True, capture_output=True)
        output = dry.stdout + dry.stderr
        assert 'rule report:' in output, output
        assert 'rule protein_properties:' not in output, output
        subprocess.run(command, cwd=root, check=True)
        assert [(p.read_bytes(), p.stat().st_mtime_ns) for p in (data, table)] == original
    print('Renderer dependency edits regenerate HTML without recomputing protein data.')
