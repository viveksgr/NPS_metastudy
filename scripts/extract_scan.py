# -*- coding: utf-8 -*-
# extract_scan.py - pull MRI acquisition sentences from each study PDF
import io, sys, re
from pathlib import Path
import fitz  # pymupdf

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
LIT = Path("Literature")

# folders we want acquisition params for (NOT RETRIEVED rows)
TARGETS = [
    "Atlas_2012_remifentanil", "Brascher_2016_control", "Bingel_2011_remifentanil",
    "Jepma_2015_conceptualConditioning", "Jepma_2018_expectancy",
    "LopezSola_2018_prosocialMeaning", "LopezSola_2019_handholding",
    "Roy_2009_emotionModulation", "Riegner_2023_mindfulness",
    "Bingel_2006_placebo", "Choi_2011_placebo", "Ellingsen_2013_placebo",
    "Elsenbruch_2012_placebo", "Freeman_2015_placebo", "Geuter_2013_placebo",
    "Kong_2009_placebo", "Lui_2010_placebo", "Schenk_2014_placebo",
    "Theysohn_2014_placebo", "Wager_2004_study1_placebo", "Wrobel_2014_placebo",
]
# accept the folder passed on cmdline, else all targets
if len(sys.argv) > 1:
    TARGETS = sys.argv[1:]

# regex for sentences likely to hold acquisition detail
KEYS = re.compile(
    r"(repetition time|echo time|flip angle|\bTR\b|\bTE\b|voxel|isotropic|"
    r"field of view|\bFOV\b|slices?|tesla|\bT\b\s|Trio|Prisma|Skyra|Allegra|Verio|"
    r"Achieva|Signa|Magnetom|Siemens|Philips|General Electric|\bGE\b|3\s?T|1\.5\s?T|"
    r"smoothing|FWHM|full[- ]width|gaussian kernel|EPI|echo[- ]planar|gradient[- ]echo|"
    r"multiband|multi[- ]band|GRAPPA|SENSE|acceleration|matrix|in[- ]plane|"
    r"acquired|functional images|BOLD|head coil|channel coil|spiral)",
    re.IGNORECASE)

def clean(t):
    t = t.replace("­", "")  # soft hyphen
    t = re.sub(r"-\n", "", t)     # de-hyphenate line breaks
    t = re.sub(r"\s+", " ", t)
    return t

for folder in TARGETS:
    fdir = LIT / folder
    pdfs = list(fdir.glob("*.pdf"))
    print("\n" + "="*90)
    print(f"### {folder}")
    if not pdfs:
        print("  (no PDF)")
        continue
    pdf = pdfs[0]
    print(f"  file: {pdf.name}")
    try:
        doc = fitz.open(pdf)
    except Exception as e:
        print(f"  OPEN ERR: {e}")
        continue
    full = " ".join(page.get_text() for page in doc)
    doc.close()
    full = clean(full)
    # split into sentences and keep those with acquisition keywords
    sents = re.split(r"(?<=[.;])\s+", full)
    hits = []
    for s in sents:
        if KEYS.search(s) and len(s) < 600:
            # require at least one number to avoid pure prose
            if re.search(r"\d", s):
                hits.append(s.strip())
    # dedupe, keep order
    seen = set(); out = []
    for h in hits:
        k = h[:80]
        if k not in seen:
            seen.add(k); out.append(h)
    for h in out[:40]:
        print("  -", h)
