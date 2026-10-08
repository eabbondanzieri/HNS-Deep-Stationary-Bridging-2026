"""
Gene start positions for every gene in the DESeq2 output, on the lab genome.

Lab genome = E. coli K-12 W3110 (GenBank AP009048.1 = RefSeq NC_007779.1)
with bases DEL_START..DEL_END (W3110 coordinates) deleted; 4,607,088 bp.
HiC, ChIP and RNA-seq data were all mapped to this genome.

Positions come from NameAndRegion.txt (from collaborators): the annotation
used to build the count tables, already in lab-genome coordinates, in
GenBank format ("start..end" or "complement(start..end)").

"Gene start" = 5' end of the annotated gene, taking strand into account
(left end for + strand genes, right end for - strand genes). This is the
start codon side, not an experimentally mapped TSS.

W3110 coordinates are also given (the deletion added back), for comparing
with W3110-based resources.

Output:  genome_positions/gene_positions.txt
         genome_positions/DESeq2_apeglm_<Norm>_<Time>_with_positions.txt
Standard library only:  python3 make_gene_positions.py
"""

import glob
import os
import re

REGION_FILE = "NameAndRegion.txt"
OUT_DIR     = "genome_positions"
GENOME_LEN  = 4607088        # lab genome length

DEL_START = 1398953          # first deleted base (W3110 coordinates, 1-based)
DEL_END   = 1438196          # last deleted base
DEL_LEN   = DEL_END - DEL_START + 1


def read_lines(path):
    # Some files have old-Mac (CR) line endings
    return open(path, newline="").read().replace("\r\n", "\n").replace("\r", "\n").splitlines()


def to_w3110(p):
    return p + DEL_LEN if p >= DEL_START else p


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    rows = []
    for line in read_lines(REGION_FILE)[1:]:
        if not line.strip():
            continue
        gid, region = line.split("\t")[:2]
        m = re.fullmatch(r"(complement\()?(\d+)\.\.(\d+)\)?", region.strip())
        if not m:
            raise ValueError(f"Unrecognized region for {gid}: {region}")
        lo, hi = int(m.group(2)), int(m.group(3))
        strand = "-" if m.group(1) else "+"
        start, end = (lo, hi) if strand == "+" else (hi, lo)
        rows.append(dict(GeneID=gid, Strand=strand,
                         Gene_start_lab=start, Gene_end_lab=end,
                         Gene_start_W3110=to_w3110(start), Gene_end_W3110=to_w3110(end),
                         Length=hi - lo + 1))

    assert all(max(r["Gene_start_lab"], r["Gene_end_lab"]) <= GENOME_LEN for r in rows)

    cols = ["GeneID", "Strand", "Gene_start_lab", "Gene_end_lab",
            "Gene_start_W3110", "Gene_end_W3110", "Length"]
    with open(f"{OUT_DIR}/gene_positions.txt", "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in rows:
            out.write("\t".join(str(r[c]) for c in cols) + "\n")

    # Add positions to each DESeq2 result file
    pos = {r["GeneID"]: r for r in rows}
    missing = set()
    for f in sorted(glob.glob("DESeq2_apeglm_*_*.txt")):
        if not re.match(r"DESeq2_apeglm_(CFU|Count|OD)_(96h|Exp)\.txt$", os.path.basename(f)):
            continue
        lines = read_lines(f)
        tag = os.path.basename(f)[:-4]
        with open(f"{OUT_DIR}/{tag}_with_positions.txt", "w") as out:
            out.write(lines[0] + "\t" + "\t".join(cols[1:]) + "\n")
            for l in lines[1:]:
                gid = l.split("\t")[0]
                r = pos.get(gid)
                if r is None:
                    missing.add(gid)
                extra = [str(r[c]) for c in cols[1:]] if r else ["NA"] * (len(cols) - 1)
                out.write(l + "\t" + "\t".join(extra) + "\n")

    print("genes with positions:", len(rows))
    print("DESeq2 genes without a position:", sorted(missing) if missing else "none")


if __name__ == "__main__":
    main()
