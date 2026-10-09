"""팀 분석 스크립트 1–23 재실행 · 결과 비교 (클로드 검증 v4, 2026-10-08)

- 원본 저장소는 건드리지 않는다: `git archive <ref>`로 만든 사본(--work)에서만 실행한다.
- 실행 전 사본의 results/ 를 results_committed/ 로 복사해 두고, 재실행 결과와 비교한다.
- 스크립트마다 종료 코드 · 소요시간(Hour) · peak 메모리(GB, 자식 프로세스 peak working set)를 기록한다.
- 출력은 matrix 형태 TSV 두 개(행 = 스크립트 / 결과 파일)와 실행 로그.

실행: conda activate HH_python3 후
  python rerun_all_20261008.py --repo <저장소> --ref main --work <임시폴더> --out <출력폴더>
"""
import argparse, glob, json, math, os, shutil, subprocess, sys, time, zipfile
import psutil

T0 = time.time()

# 실행 순서 — (사본 루트 기준 경로, 추가 인자). 16번은 NCBI 재제출 없이 저장된 결과만 요약.
# 2026-10-09 14–23 이 code/ 로 옮겨져 경로만 고침(옛 위치: 온라인멘토링/ · 최종발표_준비/staphylococcus 종 수준 확인/ ·
# 최종발표_준비/추가분석/). 그 전 커밋(--ref)을 다시 돌리려면 옛 경로가 필요하다. 10/8 기록 파일의 경로는 옛 위치 그대로.
SCRIPTS = (
    [(f"code/{n}", []) for n in [
        "1_분석표_생성.py", "2_게이트_판정.py", "3_설문라벨_검증.py", "4_조성효과_점검.py",
        "5_후보선정규칙.py", "6_오염점검.py", "7_종수준_한계.py", "8_진균_판별력.py",
        "9_나이보정_판별력.py", "10_판정표_이식성.py", "11_진균_종수준.py", "12_시간안정성.py",
        "13_다양성_분리.py"]]
    + [("code/14_피지_점검.py", [])]
    + [(f"code/{n}", a) for n, a in [
        ("15_포도상구균_종수준.py", []), ("16_포도상구균_BLAST.py", ["--summarize-only"]),
        ("17_속별_종라벨비율_261007.py", [])]]
    + [(f"code/{n}", []) for n in [
        "18_검정력_표본크기_261007.py", "19_심도_군비교_261007.py", "20_로지스틱_표준대조_261007.py",
        "21_나이보정_로그비_261007.py", "22_랜덤포레스트_비교_261007.py", "23_정규성_점검_261007.py"]]
)


def find_git():
    """GitHub Desktop 에 들어 있는 git.exe 를 찾는다 (이 PC에는 PATH 에 git 이 없음)."""
    hits = sorted(glob.glob(os.path.join(os.environ.get("LOCALAPPDATA", ""), "GitHubDesktop",
                                         "app-*", "resources", "app", "git", "cmd", "git.exe")))
    if not hits:
        sys.exit("git.exe 를 찾지 못함")
    return hits[-1]


def export_snapshot(repo, ref, work):
    """git archive 로 커밋 그대로의 사본을 만든다 (작업 트리의 미커밋 변경은 포함하지 않음)."""
    if os.path.exists(work):
        sys.exit(f"작업 폴더가 이미 있음 — 덮어쓰지 않음: {work}")
    os.makedirs(work)
    zpath = os.path.join(work, "_snapshot.zip")
    subprocess.run([find_git(), "-C", repo, "archive", "--format=zip", "-o", zpath, ref], check=True)
    src = os.path.join(work, "src")
    with zipfile.ZipFile(zpath) as z:
        z.extractall(src)
    os.remove(zpath)
    shutil.copytree(os.path.join(src, "results"), os.path.join(work, "results_committed"))
    return src


def run_one(src, rel, args, logdir, idx):
    path = os.path.join(src, rel)
    # 검수 I3 반영: BLAS 스레드를 1개로 고정해 부동소수점 순서 차이를 줄인다(팀 rerun.py 와 같은 설정)
    env = dict(os.environ, PYTHONUTF8="1", PYTHONIOENCODING="utf-8",
               OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    log = os.path.join(logdir, f"{idx:02d}.txt")
    t = time.time()
    with open(log, "w", encoding="utf-8") as fh:
        p = subprocess.Popen([sys.executable, path, *args], cwd=os.path.dirname(path),
                             stdout=fh, stderr=subprocess.STDOUT, env=env)
        # 검수 I2 반영: 즉시 끝나는 스크립트에서 NoSuchProcess 로 러너가 죽지 않게 하고,
        # joblib 자식 프로세스까지 포함한 동시 RSS 합계와 부모 peak_wset 중 큰 값을 0.1초 간격으로 기록
        try:
            proc = psutil.Process(p.pid)
        except psutil.NoSuchProcess:
            proc = None
        peak = 0
        while p.poll() is None and proc is not None:
            try:
                tree = [proc] + proc.children(recursive=True)
                peak = max(peak, proc.memory_info().peak_wset, sum(c.memory_info().rss for c in tree))
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                pass
            time.sleep(0.1)
        p.wait()
    return p.returncode, (time.time() - t) / 3600, peak / 1024**3, log


def flatten(x, prefix=""):
    """JSON 을 (경로, 값) 목록으로 펼친다."""
    if isinstance(x, dict):
        for k, v in x.items():
            yield from flatten(v, f"{prefix}/{k}")
    elif isinstance(x, list):
        for i, v in enumerate(x):
            yield from flatten(v, f"{prefix}[{i}]")
    else:
        yield prefix, x


def as_num(v):
    if isinstance(v, bool):
        return None
    if isinstance(v, (int, float)):
        return float(v)
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def compare(a, b):
    """두 결과 파일 비교 → (바이트 동일, 수치 차이 수, 최대 절대차, 최대 상대차, 비수치 차이 수)."""
    ba, bb = open(a, "rb").read(), open(b, "rb").read()
    if ba == bb:
        return 1, 0, 0.0, 0.0, 0
    if a.endswith(".json"):
        try:
            fa, fb = dict(flatten(json.loads(ba))), dict(flatten(json.loads(bb)))
        except ValueError:
            return 0, "", "", "", "parse_error"
    elif a.endswith((".tsv", ".csv")):
        sep = "\t" if a.endswith(".tsv") else ","
        rows = lambda raw: [r.split(sep) for r in raw.decode("utf-8", "replace").replace("\r\n", "\n").split("\n")]
        ra, rb = rows(ba), rows(bb)
        fa = {(i, j): c for i, r in enumerate(ra) for j, c in enumerate(r)}
        fb = {(i, j): c for i, r in enumerate(rb) for j, c in enumerate(r)}
    else:
        return 0, "", "", "", "binary_or_text_differs"
    n_num = n_other = 0
    max_abs = max_rel = 0.0
    for k in set(fa) | set(fb):
        va, vb = fa.get(k), fb.get(k)
        if va == vb:
            continue
        na, nb = as_num(va), as_num(vb)
        if na is not None and nb is not None:
            # 검수 I1 반영: 둘 다 NaN 이면 같음, 한쪽만 NaN 이면 비수치 차이로 센다
            if math.isnan(na) and math.isnan(nb):
                continue
            if math.isnan(na) or math.isnan(nb):
                n_other += 1
                continue
            d = abs(na - nb)
            if d == 0:
                continue
            n_num += 1
            max_abs = max(max_abs, d)
            max_rel = max(max_rel, d / max(abs(na), abs(nb), 1e-300))
        else:
            n_other += 1
    return 0, n_num, max_abs, max_rel, n_other


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True)
    ap.add_argument("--ref", default="main")
    ap.add_argument("--work", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    logdir = os.path.join(a.out, "rerun_logs_20261008")
    os.makedirs(logdir, exist_ok=True)

    src = export_snapshot(a.repo, a.ref, a.work)
    commit = subprocess.check_output([find_git(), "-C", a.repo, "rev-parse", a.ref], text=True).strip()
    # 검수 C1 반영: 추출 직후 results/ 파일의 수정 시각을 기록 → 재실행이 실제로 다시 쓴 파일만 「재생성」으로 센다
    res_dir = os.path.join(src, "results")
    mt0 = {n: os.stat(os.path.join(res_dir, n)).st_mtime_ns for n in os.listdir(res_dir)}
    with open(os.path.join(a.out, "rerun_status_20261008.tsv"), "w", encoding="utf-8", newline="\n") as f:
        f.write(f"# commit {commit} · python {sys.version.split()[0]} · peak_mem_gb = max(부모 peak_wset, 자식 포함 동시 RSS 합계), 0.1초 표본\n")
        f.write("script\treturncode\telapsed_hour\tpeak_mem_gb\tlog\n")
        for i, (rel, args) in enumerate(SCRIPTS, 1):
            rc, h, gb, log = run_one(src, rel, args, logdir, i)
            f.write(f"{rel}\t{rc}\t{h:.4f}\t{gb:.3f}\t{os.path.basename(log)}\n")
            f.flush()
            print(f"[{i:02d}] rc={rc} {h:.4f} h {gb:.3f} GB  {rel}", flush=True)

    old, new = os.path.join(a.work, "results_committed"), os.path.join(src, "results")
    names = sorted(n for n in set(os.listdir(old)) | set(os.listdir(new))
                   if not os.path.isdir(os.path.join(old, n)) and not os.path.isdir(os.path.join(new, n)))
    with open(os.path.join(a.out, "rerun_diff_20261008.tsv"), "w", encoding="utf-8", newline="\n") as f:
        # regenerated=0 이면 이번 재실행이 그 파일을 쓰지 않았다는 뜻 → byte_identical=1 이어도 「재현 확인」이 아님
        # numeric_ok = 비수치 차이 없음 + 수치 상대차 ≤ 1e-9 (검수 I3)
        f.write("result_file\tin_committed\tin_rerun\tregenerated\tbyte_identical\tnumeric_ok\t"
                "n_numeric_diff\tmax_abs_diff\tmax_rel_diff\tn_other_diff\n")
        for n in names:
            po, pn = os.path.join(old, n), os.path.join(new, n)
            regen = int(os.path.exists(pn) and os.stat(pn).st_mtime_ns != mt0.get(n))
            if os.path.exists(po) and os.path.exists(pn):
                row = compare(po, pn)
                ok = int(row[4] == 0 and (row[1] == 0 or row[3] <= 1e-9)) if isinstance(row[1], int) else 0
            else:
                row, ok = ("", "", "", "", "missing"), 0
            f.write("\t".join([n, str(int(os.path.exists(po))), str(int(os.path.exists(pn))), str(regen),
                               str(row[0]), str(ok)] + [str(x) for x in row[1:]]) + "\n")

    me = psutil.Process()
    print(f"total {(time.time() - T0) / 3600:.4f} Hour · runner peak {me.memory_info().peak_wset / 1024**3:.3f} GB "
          f"· python {sys.version.split()[0]} · {sys.executable}")


if __name__ == "__main__":
    main()
