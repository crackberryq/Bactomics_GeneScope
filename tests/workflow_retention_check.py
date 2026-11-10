"""Exercise retry, native sidecars, refresh and no-op retention in the real engine.

The DIAMOND executable is a filesystem stub; this tests orchestration/data
preservation only, not sequence-search correctness.
"""
from pathlib import Path
import os
import re
import shutil
import subprocess
import sys
import tempfile

from test_database_builder import DatabaseBuilderTests

ROOT = Path(__file__).resolve().parents[1]
with tempfile.TemporaryDirectory(prefix='genescope-retention-check-') as td:
    root = Path(td)
    for name in ('src', 'scripts', 'database_builder', 'envs'):
        shutil.copytree(ROOT / name, root / name)
    shutil.copyfile(ROOT / 'DatabaseBuilder.smk', root / 'DatabaseBuilder.smk')
    if os.environ.get('GENESCOPE_TEST_NO_BENCHMARK') == '1':
        workflow = root / 'DatabaseBuilder.smk'
        workflow.write_text(re.sub(r'^    benchmark:.*\n', '', workflow.read_text(), flags=re.M))
    fixture = root / 'fixture'; fixture.mkdir()
    spec = DatabaseBuilderTests().fixture(fixture)
    bin_dir = root / 'bin'; bin_dir.mkdir()
    tool = bin_dir / 'diamond'
    tool.write_text(f'#!{sys.executable}\n' + '''import os, sys
from pathlib import Path
if sys.argv[1] == 'version':
    print('DIAMOND filesystem stub; no biological validation')
    raise SystemExit(0)
assert sys.argv[1] == 'makedb'
assert '--no-unlink' in sys.argv
target = Path(sys.argv[sys.argv.index('--db') + 1])
target.write_bytes(b'synthetic-index')
temp = Path(sys.argv[sys.argv.index('--tmpdir') + 1])
(temp / 'native.extra').write_text('native tool intermediate')
marker = Path(os.environ['FAIL_DIAMOND_ONCE'])
if not marker.exists():
    marker.write_text('failed once')
    print('first-attempt failure detail', file=sys.stderr)
    raise SystemExit(1)
print('second-attempt success')
''')
    tool.chmod(0o755)
    env = dict(os.environ, PATH=str(bin_dir) + os.pathsep + str(Path(sys.executable).parent) + os.pathsep + os.environ['PATH'], FAIL_DIAMOND_ONCE=str(root / 'failed-once'))
    output = root / 'output'
    command = [sys.executable, '-m', 'snakemake', '-s', 'DatabaseBuilder.smk', '--cores', '1', '--retries', '1',
               '--config', 'reference_spec=' + str(spec), 'build_dir=' + str(output)]
    subprocess.run(command, cwd=root, env=env, check=True)
    index_attempts = list((output / 'index/retained/index').iterdir())
    assert len(index_attempts) == 2
    failed = [p for p in index_attempts if (p / 'failed.json').exists()]
    assert len(failed) == 1
    assert (failed[0] / 'proteins.dmnd').read_bytes() == b'synthetic-index'
    assert (failed[0] / 'native.extra').read_text() == 'native tool intermediate'
    assert 'first-attempt failure detail' in (failed[0] / 'stderr.log').read_text()
    (output / 'sequences/user-inspectable.native').write_text('native prior product')
    subprocess.run(command + ['--forceall'], cwd=root, env=env, check=True)
    archived = list((output / 'retained').glob('*/files/sequences/user-inspectable.native'))
    assert len(archived) == 1 and archived[0].read_text() == 'native prior product'
    assert list((output / 'retained').glob('*/files/logs/index.log'))
    assert (failed[0] / 'native.extra').read_text() == 'native tool intermediate'
    before = sorted(p.relative_to(output).as_posix() for p in output.rglob('*'))
    subprocess.run(command, cwd=root, env=env, check=True)
    subprocess.run(command + ['--forceall', '--dry-run'], cwd=root, env=env, check=True)
    after = sorted(p.relative_to(output).as_posix() for p in output.rglob('*'))
    assert before == after, 'No-op/dry-run changed retained outputs'
    print('Engine retries, native sidecars, previous logs, and no-op/dry-run preservation passed.')
