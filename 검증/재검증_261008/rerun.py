"""Reproduce tracked project scripts in an isolated git archive; never overwrite originals."""
from pathlib import Path
import hashlib, io, json, os, re, subprocess, sys, tarfile, tempfile, time

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
WORK = Path(tempfile.mkdtemp(prefix='microbiome_reaudit_261008_'))
# 2026-10-09 14–23 이 code/ 로 옮겨져 경로만 고침(옛 목록: + 온라인멘토링/14_피지_점검.py · 최종발표_준비/추가분석 · 최종발표_준비/staphylococcus 종 수준 확인)
archive = subprocess.check_output(['git', 'archive', 'HEAD', 'code', 'data', 'results'], cwd=ROOT)
with tarfile.open(fileobj=io.BytesIO(archive)) as tf:
    tf.extractall(WORK)
scripts = sorted((p for p in WORK.rglob('*.py') if re.match(r'^\d+_', p.name)), key=lambda p: int(p.name.split('_')[0]))
record = {'commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(), 'python': sys.version, 'workdir': str(WORK), 'scripts': []}
logdir = OUT / '실행로그'
logdir.mkdir(exist_ok=True)
for p in scripts:
    number = int(p.name.split('_')[0])
    args = [sys.executable, str(p)] + (['--summarize-only'] if number == 16 else [])
    start = time.monotonic()
    with (logdir / f'{number:02d}.txt').open('w') as f:
        result = subprocess.run(args, cwd=p.parent, stdout=f, stderr=subprocess.STDOUT, env={**os.environ, 'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1'})
    row = {'number': number, 'script': str(p.relative_to(WORK)), 'exit_code': result.returncode, 'seconds': round(time.monotonic()-start, 2)}
    record['scripts'].append(row)
    (OUT / '재실행_결과.json').write_text(json.dumps(record, ensure_ascii=False, indent=2))
    print(json.dumps(row, ensure_ascii=False), flush=True)
record['comparisons'] = []
for p in sorted((ROOT / 'results').iterdir()):
    other = WORK / 'results' / p.name
    if p.is_file() and other.exists():
        a, b = p.read_bytes(), other.read_bytes()
        record['comparisons'].append({'file': p.name, 'equal': a == b, 'original_sha256': hashlib.sha256(a).hexdigest(), 'rerun_sha256': hashlib.sha256(b).hexdigest()})
(OUT / '재실행_결과.json').write_text(json.dumps(record, ensure_ascii=False, indent=2))
print('COMPLETE', len(record['scripts']), sum(r['equal'] for r in record['comparisons']), '/', len(record['comparisons']), flush=True)
