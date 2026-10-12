"""Add the figure label above every page of a figure PDF without re-rendering the page (replaces add_label_pdf.swift, whose
Quartz redraw split text into single-glyph objects): each page is embedded unchanged and the page is made 42 pt taller for a
16-pt Arial Bold label, baseline 33 pt from the top, 20 pt from the left; later pages read "(continued)".
Usage: python3 add_label_pdf.py in.pdf out.pdf "Supplementary Figure 3" """
import sys
import pymupdf
BAND, SIZE, X, BASE = 42, 16, 20, 33
ARIAL_B = "/System/Library/Fonts/Supplemental/Arial Bold.ttf"
src, dst, label = sys.argv[1:4]
s = pymupdf.open(src); out = pymupdf.open()
for i, p in enumerate(s):
    w, h = p.rect.width, p.rect.height
    q = out.new_page(width=w, height=h + BAND)
    q.show_pdf_page(pymupdf.Rect(0, BAND, w, h + BAND), s, i)
    q.insert_text((X, BASE), label if i == 0 else label + " (continued)", fontsize=SIZE, fontname="ArialBoldLbl", fontfile=ARIAL_B)
out.save(dst, garbage=3, deflate=True)
print("wrote", dst)
