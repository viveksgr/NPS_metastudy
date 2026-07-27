# -*- coding: utf-8 -*-
# download_oa2.py - second-pass for remaining open-access papers
import io, sys, time, requests, re
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

LIT = Path("Literature")
H = {
    "User-Agent": ("Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
                   "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124.0.0.0 Safari/537.36"),
    "Accept": "text/html,application/xhtml+xml,*/*",
    "Accept-Language": "en-US,en;q=0.9",
}
H_PDF = dict(H, Accept="application/pdf,*/*", Referer="https://www.nature.com/")

def save(folder, fname, content):
    p = LIT / folder / fname
    p.write_bytes(content)
    print(f"  SAVED {p} [{len(content)//1024} KB]")
    return p

def get_pdf(url, label="", referer=None):
    h = dict(H_PDF)
    if referer:
        h["Referer"] = referer
    try:
        r = requests.get(url, headers=h, allow_redirects=True, timeout=40)
        if r.status_code == 200 and r.content[:4] == b"%PDF":
            print(f"  OK  {label} [{len(r.content)//1024} KB]")
            return r.content
        print(f"  FAIL {label} HTTP {r.status_code} ct={r.headers.get('Content-Type','')[:50]}")
        return None
    except Exception as e:
        print(f"  ERR  {label}: {e}")
        return None

def doi_to_url(doi):
    """Follow DOI redirect to find actual article URL."""
    try:
        r = requests.get(f"https://doi.org/{doi}", headers=H, allow_redirects=True, timeout=20)
        return r.url
    except:
        return None

# ── 1. Steininger 2025 - fix DOI ────────────────────────────────────────────
print("=== Steininger 2025 ===")
folder = "Steininger_2025_natureExposure"
if not list((LIT / folder).glob("*.pdf")):
    doi = "10.1038/s41467-025-56870-x"
    url = doi_to_url(doi)
    print(f"  DOI resolved: {url}")
    if url and "nature.com" in url:
        # Extract article ID from URL
        art_id = url.rstrip("/").split("/")[-1]
        print(f"  Article ID: {art_id}")
        pdf_url = f"https://www.nature.com/articles/{art_id}.pdf"
        c = get_pdf(pdf_url, f"NatCommun/{art_id}")
        if c:
            save(folder, f"{art_id}.pdf", c)
    time.sleep(2)

# ── 2. Tinnermann 2022 (eLife) - via eLife CDN ──────────────────────────────
print("\n=== Tinnermann 2022 (eLife 74293) ===")
folder = "Tinnermann_2022_remifentanil"
if not list((LIT / folder).glob("*.pdf")):
    # eLife CDN URL pattern
    cdn_urls = [
        "https://cdn.elifesciences.org/articles/74293/elife-74293-v2.pdf",
        "https://cdn.elifesciences.org/articles/74293/elife-74293-v1.pdf",
        "https://cdn.elifesciences.org/articles/74293/elife-74293.pdf",
    ]
    for url in cdn_urls:
        c = get_pdf(url, f"eLife CDN: {url}", referer="https://elifesciences.org/")
        if c:
            save(folder, "elife-74293.pdf", c)
            break
        time.sleep(1)
    # If still not found, try via PMC9042228
    if not list((LIT / folder).glob("*.pdf")):
        # Europe PMC XML -> get PDF link
        r = requests.get("https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9042228/fullTextXML",
                         headers=H, timeout=30)
        print(f"  EuropePMC XML: {r.status_code} {len(r.content)} bytes")
    time.sleep(2)

# ── 3. Kim 2024 (Sci Adv) via Unpaywall ─────────────────────────────────────
print("\n=== Kim 2024 Sci Adv ===")
folder = "Kim_2024_cueIntegration"
if not list((LIT / folder).glob("*.pdf")):
    doi = "10.1126/sciadv.adk7421"
    # Try Unpaywall API
    upw = requests.get(f"https://api.unpaywall.org/v2/{doi}?email=nps_metastudy@research.org",
                       timeout=15)
    if upw.status_code == 200:
        data = upw.json()
        best = data.get("best_oa_location") or {}
        pdf_url = best.get("url_for_pdf")
        print(f"  Unpaywall best OA PDF: {pdf_url}")
        if pdf_url:
            c = get_pdf(pdf_url, f"Sci Adv/{doi}", referer="https://www.science.org/")
            if c:
                save(folder, "sciadv.adk7421.pdf", c)
    time.sleep(2)

# ── 4. Riegner 2023 (PMC9823141) via Unpaywall ──────────────────────────────
print("\n=== Riegner 2023 ===")
folder = "Riegner_2023_mindfulness"
if not list((LIT / folder).glob("*.pdf")):
    doi = "10.1097/j.pain.0000000000002731"
    upw = requests.get(f"https://api.unpaywall.org/v2/{doi}?email=nps_metastudy@research.org",
                       timeout=15)
    if upw.status_code == 200:
        data = upw.json()
        best = data.get("best_oa_location") or {}
        pdf_url = best.get("url_for_pdf")
        print(f"  Unpaywall best OA PDF: {pdf_url}")
        if pdf_url:
            c = get_pdf(pdf_url, f"Pain/{doi}")
            if c:
                save(folder, "riegner2023.pdf", c)
    time.sleep(2)

# ── 5. Try Unpaywall for all remaining empty folders ─────────────────────────
print("\n=== Unpaywall pass for remaining empty folders ===")
doi_map = {
    "Eippert_2009_placebo":   "10.1016/j.neuron.2009.07.014",
    "Ellingsen_2013_placebo": "10.1073/pnas.1305050110",
    "Freeman_2015_placebo":   "10.1016/j.neuroimage.2015.03.015",
    "Geuter_2013_placebo":    "10.1016/j.neuroimage.2012.11.029",   # already has PDF
    "Kong_2006_placebo":      "10.1523/JNEUROSCI.3556-05.2006",
    "Kong_2009_placebo":      "10.1016/j.neuroimage.2008.12.025",
    "Bingel_2011_remifentanil_posExpectancy": "10.1126/scitranslmed.3001244",
    "Bingel_2011_placebo":    "10.1126/scitranslmed.3001244",       # same paper
    "Jepma_2015_conceptualConditioning": "10.1177/0956797615597658",
    "LopezSola_2019_handholding": "10.1097/j.pain.0000000000001599",
    "Wager_2004_study1_placebo": "10.1126/science.1093065",
    "Wager_2004_study2_placebo": "10.1126/science.1093065",         # same paper
    "Bingel_2006_placebo":    "10.1016/j.pain.2005.08.027",
    "Elsenbruch_2012_placebo":"10.1016/j.neuroimage.2012.07.002",
    "Lui_2010_placebo":       "10.1016/j.pain.2010.09.021",
    "Schenk_2014_placebo":    "10.1016/j.pain.2013.09.024",
    "Wrobel_2014_placebo":    "10.1016/j.cortex.2014.02.023",
    "Atlas_2012_remifentanil":"10.1523/JNEUROSCI.0383-12.2012",    # already has PDF
}
downloaded_doi = {}  # doi -> bytes

for folder, doi in doi_map.items():
    fdir = LIT / folder
    if not fdir.exists():
        continue
    if list(fdir.glob("*.pdf")):
        print(f"  SKIP {folder} (has PDF)")
        continue
    if doi in downloaded_doi:
        content = downloaded_doi[doi]
        if content:
            safe = doi.replace("/","_").replace(":","_")
            save(folder, f"{safe}.pdf", content)
        continue
    upw = requests.get(f"https://api.unpaywall.org/v2/{doi}?email=nps_metastudy@research.org",
                       timeout=15)
    if upw.status_code != 200:
        print(f"  Unpaywall {doi}: {upw.status_code}")
        downloaded_doi[doi] = None
        time.sleep(1)
        continue
    data = upw.json()
    best = data.get("best_oa_location") or {}
    pdf_url = best.get("url_for_pdf")
    oa_status = data.get("oa_status","unknown")
    print(f"  {folder}: OA={oa_status}, pdf={pdf_url}")
    if pdf_url:
        c = get_pdf(pdf_url, f"Unpaywall/{doi[:40]}")
        downloaded_doi[doi] = c
        if c:
            safe = doi.replace("/","_").replace(":","_")
            save(folder, f"{safe}.pdf", c)
    else:
        downloaded_doi[doi] = None
        print(f"    No OA PDF (status={oa_status})")
    time.sleep(2)

# ── Final status ─────────────────────────────────────────────────────────────
print("\n=== Final folder status ===")
n_ok = 0
n_missing = []
for d in sorted(LIT.iterdir()):
    if d.is_dir():
        pdfs = list(d.glob("*.pdf"))
        if pdfs:
            n_ok += 1
            print(f"  OK  {d.name}")
        else:
            n_missing.append(d.name)

print(f"\n{n_ok} folders have PDF(s)")
print(f"{len(n_missing)} still need manual download:")
for m in n_missing:
    print(f"  - {m}")
