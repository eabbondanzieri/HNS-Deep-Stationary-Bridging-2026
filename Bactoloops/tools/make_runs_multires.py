"""Headless harness for the WT 96 h multi-resolution test: bin-based tunables rescaled to keep their 10 kb physical
size (5 kb: lengths x2, areas x4; 20 kb: lengths /2, areas /4, rounded, minimum 1); donut test at HiCCUPS defaults
for the resolution (5 kb p=4 w=7; 10 kb p=2 w=5; 25 kb p=1 w=3). targetRowSum scaled with the bin count."""
import os as _hns_os
def _hns_root():
    """data/code root: $HNS_ROOT, else the repository folder (marked by .hns_root)"""
    r = _hns_os.environ.get("HNS_ROOT")
    if r:
        return _hns_os.path.join(r, "")
    d = _hns_os.path.dirname(_hns_os.path.abspath(__file__))
    while not _hns_os.path.exists(_hns_os.path.join(d, ".hns_root")) and _hns_os.path.dirname(d) != d:
        d = _hns_os.path.dirname(d)
    return _hns_os.path.join(d, "")
HNS_ROOT = _hns_root()
import re, sys
ROOT = HNS_ROOT + "Bactoloops/"
MR = ROOT + "iced_matrix_conversion/multires/"
SETS = {5: dict(targetRowSum=10000000, padding=400, medianRadius=4, rollingBallRadius=20, ridgeSmoothSigma=6,
                minParticleSize=8000, maskBlurSigma=2, finalBlurSigma=2, minSmallParticleArea=16, maxSmallParticleArea=400,
                refineRadius=2, localBackgroundRadius=100, integratedExcessHalf=5, minLoopDist=6, donutPeakWidth=4,
                donutWindow=7, apaSize=19, loopCircleRadius=12, dup=6, rc=735),
        10: dict(dup=3, rc=368),
        20: dict(targetRowSum=2500000, padding=100, medianRadius=1, rollingBallRadius=5, ridgeSmoothSigma=1.5,
                 minParticleSize=500, maskBlurSigma=0.5, finalBlurSigma=0.5, minSmallParticleArea=1, maxSmallParticleArea=25,
                 refineRadius=1, localBackgroundRadius=25, integratedExcessHalf=1, minLoopDist=2, donutPeakWidth=1,
                 donutWindow=3, apaSize=5, loopCircleRadius=3, dup=2, rc=184)}
out_dir = sys.argv[1]
base = open(ROOT + "HiC_PeakFindingAlgorithm_with_symmetrization.ijm").read()
for kb, S in SETS.items():
    m = base
    for k, v in S.items():
        if k in ("dup", "rc"):
            continue
        m, n = re.subn(rf"^({k}\s*=\s*)[0-9.]+;", rf"\g<1>{v};", m, count=1, flags=re.M)
        assert n == 1, k
    for arr in ("main", "acc"):
        m = m.replace(f"abs(b1 - {arr}B1[k]) <= 3 && abs(b2 - {arr}B2[k]) <= 3", f"abs(b1 - {arr}B1[k]) <= {S['dup']} && abs(b2 - {arr}B2[k]) <= {S['dup']}")
    a = m.index('Dialog.create("HiC Loop Finding");')
    b = m.index("skipSaving      = Dialog.getCheckbox();") + len("skipSaving      = Dialog.getCheckbox();")
    name = f"wt96_repavg_{kb}kb"
    body = m[:a] + f'sourceTitle = getTitle();\nsaveDir = "{out_dir}/";\nimageIdentifier = "{name}";\nrecenterBin = {S["rc"]};\nbinSizeBp = {kb * 1000};\nskipSaving = false;' + m[b:]
    open(f"{out_dir}/run_{name}.ijm", "w").write(f'open("{MR}{name}.tif");\n' + body +
        f'\nFile.saveString(getInfo("log"), "{out_dir}/log_{name}.txt");\neval("script", "System.exit(0);");\n')
print("ok")
