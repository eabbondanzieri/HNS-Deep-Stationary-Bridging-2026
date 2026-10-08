"""Rotation test of every loop caller: how many loops each one keeps, loses and gains when the same
circular map starts at a different genome position (rotation_test/build_rotated_inputs.py,
run_rotated.sh). Calls on the map rotated by k are mapped back to genome bins (+k mod n) and
compared with the calls on the unrotated map (k = 0), matching loops within 1 bin (circularly).

A loop on a circular chromosome has a short arc (d <= n/2 bins). On a map rotated so that the
short arc crosses the map edge, a linear caller sees the loop at linear distance n - d (near the
far corner of the map) instead of d. Lost loops are split by whether this happened.

Writes rotation_test/rotation_summary.csv (per tool and dataset) and rotation_detail.csv (per rotation).
"""
import glob
import os
import numpy as np
import pandas as pd

R = os.path.dirname(os.path.abspath(__file__)) + "/runs/"
OUT = os.path.dirname(os.path.abspath(__file__)) + "/"
N = {"wt96": 461, "bsub": 422}
SHIFTS = {"wt96": [1, 2, 3, 57, 118, 230, 333, 400], "bsub": [1, 2, 3, 97, 150, 211, 300]}
TOOLS = ["Bactoloops", "cooltools", "Chromosight", "HiCCUPS", "Mustache", "SIP", "Peakachu"]
BIN = 10000


def read(path, **kw):
    try:
        return pd.read_csv(path, sep="\t", **kw)
    except (FileNotFoundError, pd.errors.EmptyDataError):
        return None


def frame_pairs(tool, name):
    """loops in the map's own bins (rotated frame), or None if the run produced no output"""
    if tool == "Bactoloops":     # the macro already reports genome bins (recenter bin = k + 1)
        f = R + f"bactoloops/loop_calls_{name}.csv"
        if not os.path.exists(f):
            return None
        L = pd.read_csv(f)
        return ("genome", np.c_[L.genome_bin1.round().astype(int), L.genome_bin2.round().astype(int)])
    cols = {"cooltools": (f"cooltools/{name}.tsv", "start1", "start2", {}),
            "Chromosight": (f"chromosight/{name}.tsv", "start1", "start2", {}),
            "Mustache": (f"mustache/{name}.tsv", "BIN1_START", "BIN2_START", {}),
            "Peakachu": (f"peakachu/{name}.tsv", 1, 4, {"header": None}),
            "HiCCUPS": (f"hiccups/{name}/merged_loops.bedpe", 1, 4, {"header": None, "comment": "#"})}
    if tool == "SIP":
        f = glob.glob(R + f"sip/{name}/*Loops*.txt")
        T = read(f[0]) if f else None
        c1, c2 = "x1", "y1"
    else:
        f, c1, c2, kw = cols[tool]
        T = read(R + f, **kw)
    if T is None:
        return None
    if len(T) == 0:
        return ("frame", np.zeros((0, 2), int))
    a = np.round(T[c1].astype(float).values / BIN).astype(int)
    b = np.round(T[c2].astype(float).values / BIN).astype(int)
    return ("frame", np.c_[a, b])


def genome_set(tool, d, k):
    r = frame_pairs(tool, f"{d}_r{k}")
    if r is None:
        return None
    kind, P = r
    n = N[d]
    if kind == "frame":
        P = (P + k) % n
    return {tuple(sorted(p)) for p in P.tolist()}


def near(p, S, n):
    return any(max(min(abs(p[i] - q[i]), n - abs(p[i] - q[i])) for i in (0, 1)) <= 1 for q in S)


def crosses_edge(p, k, n):
    """does the loop's short arc cross the edge of the map rotated by k?"""
    a, b = (p[0] - k) % n, (p[1] - k) % n
    lin = abs(a - b)
    return lin > n - lin


rows = []
for d in N:
    n = N[d]
    for tool in TOOLS:
        S0 = genome_set(tool, d, 0)
        if S0 is None:
            continue
        for k in SHIFTS[d]:
            Sk = genome_set(tool, d, k)
            if Sk is None:
                rows.append(dict(dataset=d, tool=tool, k=k, missing=True))
                continue
            lost = [p for p in S0 if not near(p, Sk, n)]
            gained = [p for p in Sk if not near(p, S0, n)]
            rows.append(dict(dataset=d, tool=tool, k=k, missing=False, n0=len(S0), nk=len(Sk),
                             exact=len(S0 & Sk), kept=len(S0) - len(lost), lost=len(lost), gained=len(gained),
                             lost_edge=sum(crosses_edge(p, k, n) for p in lost),
                             edge_loops0=sum(crosses_edge(p, k, n) for p in S0)))
D = pd.DataFrame(rows)
D.to_csv(OUT + "rotation_detail.csv", index=False)
ok = D[~D.missing]
S = ok.groupby(["dataset", "tool"]).agg(n0=("n0", "first"), rotations=("k", "count"),
        kept_mean=("kept", "mean"), kept_min=("kept", "min"), lost_mean=("lost", "mean"), lost_max=("lost", "max"),
        gained_mean=("gained", "mean"), gained_max=("gained", "max"), lost_edge_total=("lost_edge", "sum"),
        lost_total=("lost", "sum"), edge_loops_total=("edge_loops0", "sum")).reset_index()
S["pct_kept_mean"] = (100 * S.kept_mean / S.n0).round(1)
S["pct_kept_min"] = (100 * S.kept_min / S.n0).round(1)
S["changed_mean_pct"] = (100 * (S.lost_mean + S.gained_mean) / S.n0).round(1)
S.to_csv(OUT + "rotation_summary.csv", index=False)
pd.set_option("display.width", 250)
print(S.to_string(index=False))
if D.missing.any():
    print("missing runs:", D[D.missing][["dataset", "tool", "k"]].values.tolist())
