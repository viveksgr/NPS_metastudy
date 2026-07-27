# -*- coding: utf-8 -*-
# download_oa3.py - third pass
import io, sys, time, requests
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
LIT = Path("Literature")

def H(ref=None, accept="application/pdf,*/*"):
    h = {
        "User-Agent": ("Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
                       "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124.0.0.0 Safari/537.36"),
        "Accept": accept,
        "Accept-Language": "en-US,en;q=0.9",
    }
    if ref:
        h["Referer"] = ref
    return h

def get_pdf(url, label="", ref=None):
    try:
        r = requests.get(url, headers=H(ref), allow_redirects=True, timeout=40, stream=True)
        ct = r.headers.get("Content-Type","")
        data = r.content
        if r.status_code == 200 and (data[:4] == b"%PDF" or "pdf" in ct):
            print(f"  OK  {label} [{len(data)//1024} KB]")
            return data
        print(f"  FAIL {label} HTTP {r.status_code} ct={ct[:60]}")
        return None
    except Exception as e:
        print(f"  ERR  {label}: {e}")
        return None

def save(folder, fname, content):
    p = LIT / folder / fname
    p.write_bytes(content)
    print(f"  SAVED {p}")

def skip(folder):
    return bool(list((LIT/folder).glob("*.pdf")))

# ── Steininger 2025: follow DOI redirect properly ───────────────────────────
print("=== Steininger 2025 ===")
if not skip("Steininger_2025_natureExposure"):
    # DOI resolves to nature.com but PDF isn't found - check actual article ID
    r = requests.get("https://doi.org/10.1038/s41467-025-56870-x",
                     headers=H(accept="text/html,*/*"), allow_redirects=True, timeout=20)
    print(f"  DOI final URL: {r.url}  status={r.status_code}")
    # The resolved URL gives us the real article ID
    import re
    # look for the canonical article slug
    canonical = re.search(r'nature\.com/articles/([^\s"\'/?]+)', r.url)
    art_id = canonical.group(1) if canonical else "s41467-025-56870-x"
    print(f"  Article ID from redirect: {art_id}")
    c = get_pdf(f"https://www.nature.com/articles/{art_id}.pdf",
                f"NatCommun/{art_id}", ref="https://www.nature.com/")
    if c:
        save("Steininger_2025_natureExposure", f"{art_id}.pdf", c)
    time.sleep(2)

# ── Eippert 2009 via Cell Press PDF ─────────────────────────────────────────
print("\n=== Eippert 2009 (Neuron, Cell Press) ===")
if not skip("Eippert_2009_placebo"):
    c = get_pdf("https://www.cell.com/neuron/pdf/S0896-6273(09)00543-1.pdf",
                "Cell/Neuron Eippert", ref="https://www.cell.com/")
    if not c:
        c = get_pdf("https://www.cell.com/article/S0896627309005431/pdf",
                    "Cell redirect", ref="https://www.cell.com/")
    if c:
        save("Eippert_2009_placebo", "Eippert_2009_Neuron.pdf", c)
    time.sleep(2)

# ── Ellingsen 2013 via PNAS ──────────────────────────────────────────────────
print("\n=== Ellingsen 2013 (PNAS) ===")
if not skip("Ellingsen_2013_placebo"):
    urls = [
        "https://www.pnas.org/doi/pdf/10.1073/pnas.1305050110",
        "https://www.pnas.org/content/pnas/110/44/17993.full.pdf",
    ]
    for url in urls:
        c = get_pdf(url, f"PNAS Ellingsen", ref="https://www.pnas.org/")
        if c:
            save("Ellingsen_2013_placebo", "Ellingsen_2013_PNAS.pdf", c)
            break
        time.sleep(1)
    time.sleep(2)

# ── Kong 2006 via J Neurosci ─────────────────────────────────────────────────
print("\n=== Kong 2006 (J Neurosci) ===")
if not skip("Kong_2006_placebo"):
    urls = [
        "https://www.jneurosci.org/content/26/2/381.full.pdf",
        "https://www.jneurosci.org/content/jneuro/26/2/381.full.pdf",
    ]
    for url in urls:
        c = get_pdf(url, "JNeurosci Kong", ref="https://www.jneurosci.org/")
        if c:
            save("Kong_2006_placebo", "Kong_2006_JNeurosci.pdf", c)
            break
        time.sleep(1)
    time.sleep(2)

# ── Kim 2024 Sci Adv via PMC OA API ─────────────────────────────────────────
print("\n=== Kim 2024 (Sci Adv, PMC11389792) ===")
if not skip("Kim_2024_cueIntegration"):
    # Try Unpaywall locations (not just best)
    doi = "10.1126/sciadv.adk7421"
    upw = requests.get(f"https://api.unpaywall.org/v2/{doi}?email=nps_metastudy@research.org",
                       timeout=15)
    locs = upw.json().get("oa_locations",[]) if upw.status_code==200 else []
    print(f"  Unpaywall locations: {len(locs)}")
    for loc in locs:
        purl = loc.get("url_for_pdf")
        if purl:
            print(f"    trying: {purl}")
            c = get_pdf(purl, f"Sci Adv/{purl[:50]}")
            if c:
                save("Kim_2024_cueIntegration", "sciadv.adk7421.pdf", c)
                break
            time.sleep(1)
    # Also try Sci Adv open access direct URL variant
    if not skip("Kim_2024_cueIntegration"):
        c = get_pdf("https://www.science.org/doi/epdf/10.1126/sciadv.adk7421",
                    "SciAdv epdf", ref="https://www.science.org/")
        if c:
            save("Kim_2024_cueIntegration", "sciadv.adk7421.pdf", c)
    time.sleep(2)

# ── Jepma 2015 (Psychol Sci, SAGE) ──────────────────────────────────────────
print("\n=== Jepma 2015 (Psychol Sci) ===")
if not skip("Jepma_2015_conceptualConditioning"):
    urls = [
        "https://journals.sagepub.com/doi/pdf/10.1177/0956797615597658",
        "https://journals.sagepub.com/doi/reader/10.1177/0956797615597658",
    ]
    for url in urls:
        c = get_pdf(url, "SAGE Jepma 2015", ref="https://journals.sagepub.com/")
        if c:
            save("Jepma_2015_conceptualConditioning", "Jepma_2015_PsychSci.pdf", c)
            break
        time.sleep(1)
    time.sleep(2)

# ── Freeman 2015 (NeuroImage, PMC4408248) - try PubMed OA ──────────────────
print("\n=== Freeman 2015 (NeuroImage) ===")
if not skip("Freeman_2015_placebo"):
    # Elsevier open archive
    urls = [
        "https://www.sciencedirect.com/science/article/pii/S1053811915002700/pdfft?isDTMRedir=true",
        "https://pdf.sciencedirectassets.com/272508/1-s2.0-S1053811915X00079/1-s2.0-S1053811915002700/main.pdf",
    ]
    for url in urls:
        c = get_pdf(url, f"Elsevier Freeman", ref="https://www.sciencedirect.com/")
        if c:
            save("Freeman_2015_placebo", "Freeman_2015_NeuroImage.pdf", c)
            break
        time.sleep(1)
    time.sleep(2)

# ── Summary ──────────────────────────────────────────────────────────────────
print("\n=== FINAL STATUS ===")
ok, miss = [], []
for d in sorted(LIT.iterdir()):
    if d.is_dir():
        pdfs = list(d.glob("*.pdf"))
        (ok if pdfs else miss).append(d.name)
print(f"Have PDF: {len(ok)}")
print(f"Missing:  {len(miss)}")
for m in miss:
    print(f"  - {m}")
