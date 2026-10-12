import os, sys, runpy, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.figure, matplotlib.pyplot as plt
matplotlib.figure.Figure.savefig = lambda *a, **k: None
plt.savefig = lambda *a, **k: None
def run(path, args=()):
    d = os.path.dirname(path); cwd = os.getcwd(); argv = sys.argv
    os.chdir(d); sys.argv = [path, *args]; sys.path.insert(0, d)
    import io, contextlib, glob, shutil, tempfile
    keep = tempfile.mkdtemp()                      # figure PDFs in the script's folder, restored if the run leaves them empty
    saved = {p: shutil.copy2(p, keep) for p in glob.glob(os.path.join(d, "*.pdf")) if os.path.getsize(p) > 2000}
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            ns = runpy.run_path(path, run_name="__main__")
    finally:
        for p, c in saved.items():
            if not os.path.exists(p) or os.path.getsize(p) < 2000:
                shutil.copy2(c, p); print(f"_ns_probe: restored {os.path.basename(p)} (left empty by the run)", file=sys.stderr)
        shutil.rmtree(keep, ignore_errors=True)
        os.chdir(cwd); sys.argv = argv; sys.path.remove(d)
    plt.close("all")
    return ns
if __name__ == "__main__":
    ns = run(sys.argv[1], sys.argv[2:])
    for k, v in ns.items():
        if k.startswith("_"): continue
        if isinstance(v, pd.DataFrame): print(f"DF {k}: {v.shape} {list(v.columns)[:12]}")
        elif isinstance(v, np.ndarray) and v.size > 1: print(f"ARR {k}: {v.shape}")
        elif isinstance(v, dict) and v and all(isinstance(x, (np.ndarray, pd.DataFrame, list, tuple)) for x in v.values()):
            print(f"DICT {k}: keys {list(v)[:10]}")
        elif isinstance(v, list) and len(v) > 1 and not callable(v[0]): print(f"LIST {k}: len {len(v)} type {type(v[0]).__name__}")
