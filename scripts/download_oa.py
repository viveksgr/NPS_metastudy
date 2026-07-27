# -*- coding: utf-8 -*-
# download_oa.py  - retry open-access papers with correct URL patterns
import io, sys, time, requests
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

LIT = Path("Literature")

H_CHROME = {
    "User-Agent": ("Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
                   "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124.0.0.0 Safari/537.36"),
    "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
    "Accept-Language": "en-US,en;q=0.5",
}
H_NAT = dict(H_CHROME, Referer="https://www.nature.com/",
             Accept="application/pdf,text/html,*/*")

def save(folder, fname, content):
    p = LIT / folder / fname
    p.write_bytes(content)
    print(f"  SAVED {p} [{len(content)//1024} KB]")
    return p

def get(url, headers=None, label=""):
    h = headers or H_CHROME
    try:
        r = requests.get(url, headers=h, allow_redirects=True, timeout=40)
        ct = r.headers.get("Content-Type", "")
        if r.status_code == 200 and (r.content[:4] == b"%PDF" or "pdf" in ct):
            print(f"  OK  {label} [{len(r.content)//1024} KB]")
            return r.content
        print(f"  FAIL {label} HTTP {r.status_code} ct={ct[:50]}")
        return None
    except Exception as e:
        print(f"  ERR  {label}: {e}")
        return None

# ── Nature Commun / Nature Publishing Group (open access) ──────────────────
nat_papers = [
    ("Habermann_2025_control",     "s41467-025-66038-7"),
    ("Habermann_2025_predictability", "s41467-025-66038-7"),   # same paper
    ("Woo_2017_expectancy",        "ncomms14211"),
    ("Woo_2017_perceivedControl",  "ncomms14211"),             # same paper
    ("Koban_2019_socialExpectancy","s41467-019-11934-y"),
    ("Steininger_2025_natureExposure", "s41467-025-56870-x"),
    ("BotvinikNezer_2024_placebo", "s41467-024-50103-8"),
]
downloaded_nat = {}  # article_id -> bytes (avoid double-fetch)

print("=== Nature Commun PDFs ===")
for folder, art_id in nat_papers:
    if list((LIT / folder).glob("*.pdf")):
        print(f"  SKIP {folder} (already has PDF)")
        continue
    url = f"https://www.nature.com/articles/{art_id}.pdf"
    if art_id in downloaded_nat:
        content = downloaded_nat[art_id]
        print(f"  REUSE {art_id} for {folder}")
    else:
        content = get(url, H_NAT, f"NatCommun/{art_id}")
        if content:
            downloaded_nat[art_id] = content
        time.sleep(2)
    if content:
        save(folder, f"{art_id}.pdf", content)

# ── Science Advances (Kim 2024) ─────────────────────────────────────────────
print("\n=== Science Advances ===")
folder = "Kim_2024_cueIntegration"
if not list((LIT / folder).glob("*.pdf")):
    url = "https://www.science.org/doi/pdf/10.1126/sciadv.adk7421"
    c = get(url, dict(H_CHROME, Referer="https://www.science.org/"), "SciAdv/adk7421")
    if c:
        save(folder, "sciadv.adk7421.pdf", c)
    time.sleep(2)

# ── eLife (Tinnermann 2022) ─────────────────────────────────────────────────
print("\n=== eLife ===")
folder = "Tinnermann_2022_remifentanil"
if not list((LIT / folder).glob("*.pdf")):
    # Try CDN-style URL via DOI redirect
    import re
    r = requests.get("https://doi.org/10.7554/eLife.74293",
                     headers=H_CHROME, allow_redirects=True, timeout=30)
    print(f"  eLife DOI resolved to: {r.url}")
    # Try article page to find PDF link
    r2 = requests.get("https://elifesciences.org/articles/74293",
                      headers=H_CHROME, allow_redirects=True, timeout=30)
    matches = re.findall(r'(https?://[^\s"\'<>]+\.pdf)', r2.text[:8000])
    print(f"  eLife PDF links found: {matches[:5]}")
    tried = set()
    for m in matches[:3]:
        if m in tried:
            continue
        tried.add(m)
        c = get(m, H_CHROME, f"eLife/{m[:60]}")
        if c:
            save(folder, "eLife.74293.pdf", c)
            break
        time.sleep(1)
    if not list((LIT / folder).glob("*.pdf")):
        # Try PMC version (Tinnermann 2022 may have a PMC ID)
        r3 = requests.get(
            "https://www.ebi.ac.uk/europepmc/webservices/rest/search?"
            "query=DOI:10.7554/eLife.74293&format=json",
            timeout=15)
        try:
            data = r3.json()
            res = data.get("resultList",{}).get("result",[])
            pmc = res[0].get("pmcid") if res else None
            print(f"  eLife PMC ID: {pmc}")
            if pmc:
                c = get(f"https://europepmc.org/backend/ptpmcrender.fcgi?accid={pmc}&blobtype=pdf",
                        H_CHROME, f"EuropePMC/{pmc}")
                if c:
                    save(folder, "eLife.74293.pdf", c)
        except:
            pass
    time.sleep(2)

# ── PMC papers via NCBI direct ───────────────────────────────────────────────
print("\n=== PMC downloads (NCBI) ===")
pmc_papers = [
    ("Eippert_2009_placebo",   "PMC6670627"),
    ("Ellingsen_2013_placebo", "PMC3816412"),
    ("Freeman_2015_placebo",   "PMC4408248"),
    ("Kong_2006_placebo",      "PMC6674420"),
    ("Kong_2009_placebo",      "PMC2737445"),
    ("Riegner_2023_mindfulness","PMC9823141"),
]

H_PMC = dict(H_CHROME, Referer="https://www.ncbi.nlm.nih.gov/",
             Accept="application/pdf,*/*")

for folder, pmcid in pmc_papers:
    if list((LIT / folder).glob("*.pdf")):
        print(f"  SKIP {folder}")
        continue
    urls = [
        f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/pdf/",
        f"https://europepmc.org/backend/ptpmcrender.fcgi?accid={pmcid}&blobtype=pdf",
    ]
    for url in urls:
        c = get(url, H_PMC, f"{pmcid}")
        if c:
            save(folder, f"{pmcid}.pdf", c)
            break
        time.sleep(2)
    time.sleep(2)

# ── Summary ──────────────────────────────────────────────────────────────────
print("\n=== Final status ===")
for d in sorted(LIT.iterdir()):
    if d.is_dir():
        pdfs = list(d.glob("*.pdf"))
        supps = [p for p in d.glob("supp_*")]
        status = f"{len(pdfs)} PDF(s)"
        if supps:
            status += f", {len(supps)} supp"
        print(f"  {d.name:<55} {status}")
