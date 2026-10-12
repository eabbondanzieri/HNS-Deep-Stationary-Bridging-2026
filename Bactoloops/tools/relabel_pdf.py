"""Replace the existing "Figure N" / "Supplementary Figure N" label on page 1 of a figure PDF without re-rendering the page
(replaces relabel_pdf.swift / relabel_main_pdf.swift, whose Quartz redraw split text into single-glyph objects): the old label
text is removed (text only) and the new one is set in 16-pt Arial Bold, baseline 33 pt from the top, 20 pt from the left.
Usage: python3 relabel_pdf.py in.pdf out.pdf "Supplementary Figure 5" """
import sys
import pymupdf
ARIAL_B = "/System/Library/Fonts/Supplemental/Arial Bold.ttf"
src, dst, label = sys.argv[1:4]
d = pymupdf.open(src); p = d[0]; n = 0
for b in p.get_text("dict")["blocks"]:
    for l in b.get("lines", []):
        for s in l["spans"]:
            if s["size"] > 13 and s["bbox"][1] < 80 and ("Figure" in s["text"] or s["text"].strip().isdigit()):
                p.add_redact_annot(pymupdf.Rect(s["bbox"])); n += 1
if n:
    p.apply_redactions(images=pymupdf.PDF_REDACT_IMAGE_NONE, graphics=pymupdf.PDF_REDACT_LINE_ART_NONE)
p.insert_text((20, 33), label, fontsize=16, fontname="ArialBoldLbl", fontfile=ARIAL_B)
d.save(dst, garbage=3, deflate=True)
print(f"wrote {dst} (removed {n} old label span(s))")
