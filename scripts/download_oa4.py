# -*- coding: utf-8 -*-
# download_oa4.py - Steininger PMC lookup + Kim 2024 via EuroPMC + AAAS CDN
import io, sys, time, requests, re
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
LIT = Path("Literature")

H = {
    "User-Agent": ("Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
                   "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124.0.0.0 Safari/537.36"),
    "Accept": "text/html,*/*",
    "Accept-Language": "en-US,en;q=0.9",
}
H_PDF = dict(H, Accept="application/pdf,*/*")

def get_pdf(url, label="", ref=None):
    h = dict(H_PDF)
    if ref: h["Referer"] = ref
    try:
        r = requests.get(url, headers=h, allow_redirects=True, timeout=40)
        data = r.content
        if r.status_code == 200 and (data[:4] == b"%PDF" or "pdf" in r.headers.get("Content-Type","")):
            print(f"  OK  {label} [{len(data)//1024} KB]")
            return data
        print(f"  FAIL {label} HTTP {r.status_code} ct={r.headers.get('Content-Type','')[:60]}")
        return None
    except Exception as e:
        print(f"  ERR  {label}: {e}")
        return None

def save(folder, fname, content):
    p = LIT / folder / fname
    p.write_bytes(content)
    print(f"  SAVED {p} [{len(content)//1024} KB]")

# ── Steininger 2025 via PMC ID ───────────────────────────────────────────────
print("=== Steininger 2025 (PMC11906725) ===")
if not list((LIT/"Steininger_2025_natureExposure").glob("*.pdf")):
    # Get article metadata from EuropePMC to find correct DOI/URL
    r = requests.get("https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11906725"
                     "?format=json", headers=H, timeout=20)
    print(f"  EuropePMC metadata status: {r.status_code}")
    if r.status_code == 200:
        try:
            d = r.json()
            doi = d.get("doi") or d.get("DOI")
            title = d.get("title","")
            source = d.get("source","")
            print(f"  DOI: {doi}  source: {source}")
            print(f"  Title: {title[:80]}")
            if doi:
                # Try Nature PDF with this DOI
                art = doi.split("10.1038/")[-1] if "10.1038/" in doi else doi.replace("/","_")
                c = get_pdf(f"https://www.nature.com/articles/{art}.pdf",
                            f"NatCommun/{art}", ref="https://www.nature.com/")
                if c:
                    save("Steininger_2025_natureExposure", f"{art}.pdf", c)
        except Exception as e:
            print(f"  parse err: {e}")
    # Also try the fullTextXML endpoint to find the PDF link
    r2 = requests.get("https://www.ebi.ac.uk/europepmc/webservices/rest/"
                      "PMC11906725/fullTextXML", headers=H, timeout=30)
    print(f"  EuropePMC XML: {r2.status_code} {len(r2.content)} bytes")
    time.sleep(2)

# ── Kim 2024 Sci Adv (PMC11389792) ──────────────────────────────────────────
print("\n=== Kim 2024 Sci Adv (PMC11389792) ===")
if not list((LIT/"Kim_2024_cueIntegration").glob("*.pdf")):
    # EuroPMC metadata
    r = requests.get("https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11389792"
                     "?format=json", headers=H, timeout=20)
    print(f"  EuropePMC metadata: {r.status_code}")
    if r.status_code == 200:
        try:
            d = r.json()
            print(f"  doi={d.get('doi')}  fullTextUrlList={d.get('fullTextUrlList','')}")
        except:
            pass
    # Try AAAS Sci Adv CDN patterns
    for pattern in [
        "https://www.science.org/doi/pdf/10.1126/sciadv.adk7421",
        "https://www.science.org/doi/reader/10.1126/sciadv.adk7421",
        "https://advances.sciencemag.org/content/10/40/eadk7421.full.pdf",
    ]:
        c = get_pdf(pattern, f"SciAdv/{pattern[-30:]}", ref="https://www.science.org/")
        if c:
            save("Kim_2024_cueIntegration", "sciadv.adk7421.pdf", c)
            break
        time.sleep(1)
    time.sleep(2)

# ── Final summary + write download_list.txt ──────────────────────────────────
print("\n=== Writing manual download list ===")
missing = {
    "Bingel_2006_placebo":             ("10.1016/j.pain.2005.08.027", "Pain 2006 (paywalled)"),
    "Bingel_2011_placebo":             ("10.1126/scitranslmed.3001244", "Sci Transl Med 2011 (paywalled) — same PDF as Bingel_2011_remifentanil"),
    "Bingel_2011_remifentanil_posExpectancy": ("10.1126/scitranslmed.3001244", "Sci Transl Med 2011 (paywalled)"),
    "Choi_2011_placebo":               ("NeuroReport 2011 Choi JC et al.", "Citation unclear — verify exact journal"),
    "Eippert_2009_placebo":            ("10.1016/j.neuron.2009.07.014", "Neuron 2009 — has PMC6670627 but NCBI bot-blocked"),
    "Ellingsen_2013_placebo":          ("10.1073/pnas.1305050110", "PNAS 2013 — has PMC3816412 but NCBI bot-blocked"),
    "Elsenbruch_2012_placebo":         ("10.1016/j.neuroimage.2012.07.002", "NeuroImage 2012 (paywalled)"),
    "Freeman_2015_placebo":            ("10.1016/j.neuroimage.2015.03.015", "NeuroImage 2015 — PMC4408248 bot-blocked; try PMC or institutional access"),
    "Jepma_2015_conceptualConditioning":("10.1177/0956797615597658", "Psychol Sci 2015 — publicly listed but SAGE blocks bots"),
    "Kim_2024_cueIntegration":         ("10.1126/sciadv.adk7421", "Sci Adv 2024 (open access) — blocked bot-detection; try https://www.science.org/doi/10.1126/sciadv.adk7421"),
    "Kong_2006_placebo":               ("10.1523/JNEUROSCI.3556-05.2006", "J Neurosci 2006 — PMC6674420 bot-blocked"),
    "Kong_2009_placebo":               ("10.1016/j.neuroimage.2008.12.025", "NeuroImage 2009 — paywalled; PMC2737445 bot-blocked"),
    "Lui_2010_placebo":                ("10.1016/j.pain.2010.09.021", "Pain 2010 (paywalled)"),
    "Schenk_2014_placebo":             ("10.1016/j.pain.2013.09.024", "Pain 2014 (paywalled)"),
    "Steininger_2025_natureExposure":  ("10.1038/s41467-025-56870-x", "Nat Commun 2025 — DOI may be incorrect; verify on PubMed"),
    "Theysohn_2014_placebo":           ("Neurogastroenterol Motil ~2014, Theysohn et al.", "Citation unclear — verify exact reference"),
    "Wager_2004_study1_placebo":       ("10.1126/science.1093065", "Science 2004 (paywalled) — Wager_2004_study2 is same paper"),
    "Wager_2004_study2_placebo":       ("10.1126/science.1093065", "Science 2004 — same paper as study1"),
    "Wrobel_2014_placebo":             ("10.1016/j.cortex.2014.02.023", "Cortex 2014 (paywalled)"),
}

lines = ["MANUAL DOWNLOAD LIST\n"+"="*60+"\n",
         "Papers below could not be downloaded automatically.\n",
         "Either paywalled or publisher bot-detection blocked download.\n",
         "Use institutional access / VPN, or download from the DOI/link.\n\n"]
for folder, (doi, note) in sorted(missing.items()):
    if not list((LIT/folder).glob("*.pdf")):
        lines.append(f"FOLDER:  {folder}\n")
        lines.append(f"DOI/ref: {doi}\n")
        lines.append(f"Notes:   {note}\n\n")

(LIT/"MANUAL_DOWNLOAD_LIST.txt").write_text("".join(lines), encoding="utf-8")
print("  Wrote Literature/MANUAL_DOWNLOAD_LIST.txt")

ok = sum(1 for d in LIT.iterdir() if d.is_dir() and list(d.glob("*.pdf")))
miss = sum(1 for d in LIT.iterdir() if d.is_dir() and not list(d.glob("*.pdf")))
print(f"\n  TOTAL: {ok} folders with PDF, {miss} without")
