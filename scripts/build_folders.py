# -*- coding: utf-8 -*-
# build_folders.py
# Creates per-study subfolders under Literature/, moves existing loose PDFs,
# and attempts to download main PDFs (and supplements where possible) from
# open-access sources (Europe PMC, eLife, Nat Commun, Sci Adv).
# Run from: the NPS_metastudy directory.
# Requires: pip install requests

import os, sys, shutil, time, re, json, io
import requests
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

LIT   = Path("Literature")
HEADS = {
    "User-Agent": ("Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
                   "AppleWebKit/537.36 (KHTML, like Gecko) "
                   "Chrome/124.0.0.0 Safari/537.36"),
    "Accept": "application/pdf,*/*",
}

# ── study registry ────────────────────────────────────────────────────────────
# Each entry: (folder_name, pmc_id, extra_pdf_url, existing_loose_file, shared_with_folder)
#   pmc_id           : "PMC1234567" or None
#   extra_pdf_url    : direct PDF url (e.g. eLife) or None
#   existing_loose   : filename already in Literature/ root, or None
#   shared_with      : if the PDF is the SAME paper as another folder, give that
#                      folder name so we copy rather than re-download
STUDIES = [
    # ── Non-placebo studies ──────────────────────────────────────────────────
    ("Atlas_2010_expectancy",                  "PMC2966558",  None, "atlas_2010_exp.pdf",               None),
    ("Atlas_2012_remifentanil",                "PMC3387557",  None, "atlas_2013_remi_open_hidden.pdf",  None),
    ("Brascher_2016_control",                  "PMC6601855",  None, "becker_2016_pain_control.pdf",     None),
    ("Becker_2017_reward",                     "PMC5390724",  None, "becker_2017_pain_reward.pdf",      None),
    ("Bingel_2011_remifentanil_posExpectancy", None,          None, None,                               None),  # Sci Transl Med – paywalled
    ("Woo_2015_selfregulation",                "PMC4285399",  None, "bmrk3.pdf",                        None),
    ("Habermann_2025_control",                 "PMC12627469", None, None,                               None),
    ("Habermann_2025_predictability",          "PMC12627469", None, None,          "Habermann_2025_control"),  # same paper
    ("Jepma_2015_conceptualConditioning",      None,          None, None,                               None),  # Psychol Sci – paywalled
    ("Woo_2017_expectancy",                    None,          None, None,                               None),  # Nat Commun – try PMC lookup
    ("Woo_2017_perceivedControl",              None,          None, None,              "Woo_2017_expectancy"),  # same paper
    ("Jepma_2018_expectancy",                  "PMC6768437",  None, "jepma_2018_ie2.pdf",               None),
    ("Koban_2019_socialExpectancy",            "PMC6736972",  None, "koban_2019_scebl_social_pain.pdf", None),
    ("Kober_2019_mindfulAcceptance",           "PMC7057281",  None, "kober_2019_mindful_acceptance_mrp.pdf", None),
    ("LopezSola_2018_prosocialMeaning",        "PMC6218300",  None, "lopezsola_2018_pain_meaning.pdf",  None),
    ("LopezSola_2019_handholding",             None,          None, "lopezsola_2019_handholding.pdf",   None),  # Pain – paywalled but file exists
    ("BotvinikNezer_2024_placebo",             "PMC11255344", None, None,                               None),
    ("Roy_2009_emotionModulation",             "PMC2779826",  None, "roy_emomod_2009.pdf",              None),
    ("Steininger_2025_natureExposure",         "PMC11906725", None, None,                               None),
    ("Kim_2024_cueIntegration",                "PMC11389792", None, None,                               None),
    ("Tinnermann_2022_remifentanil",           None,
        "https://elifesciences.org/articles/74293",           None,                                     None),
    ("Riegner_2023_mindfulness",               "PMC9823141",  None, None,                               None),
    # ── Placebo / Zunhammer studies ──────────────────────────────────────────
    ("Bingel_2006_placebo",                    None,          None, None,                               None),  # Pain 2006 – paywalled
    ("Bingel_2011_placebo",                    None,          None, None, "Bingel_2011_remifentanil_posExpectancy"),  # same paper
    ("Choi_2011_placebo",                      None,          None, None,                               None),  # unclear/paywalled
    ("Eippert_2009_placebo",                   "PMC6670627",  None, None,                               None),
    ("Ellingsen_2013_placebo",                 "PMC3816412",  None, None,                               None),
    ("Elsenbruch_2012_placebo",                None,          None, None,                               None),  # NeuroImage – paywalled
    ("Freeman_2015_placebo",                   "PMC4408248",  None, None,                               None),
    ("Geuter_2013_placebo",                    "PMC3578963",  None, "geuter_2013_placebo.pdf",          None),
    ("Kong_2006_placebo",                      "PMC6674420",  None, None,                               None),
    ("Kong_2009_placebo",                      "PMC2737445",  None, None,                               None),
    ("Lui_2010_placebo",                       None,          None, None,                               None),  # Pain – paywalled
    ("Schenk_2014_placebo",                    None,          None, None,                               None),  # Pain – paywalled
    ("Theysohn_2014_placebo",                  None,          None, None,                               None),  # unclear
    ("Wager_2004_study1_placebo",              None,          None, None,                               None),  # Science – paywalled
    ("Wager_2004_study2_placebo",              None,          None, None,   "Wager_2004_study1_placebo"),  # same paper
    ("Wrobel_2014_placebo",                    None,          None, None,                               None),  # Cortex – paywalled
]

# ── helpers ──────────────────────────────────────────────────────────────────

def make_folder(name):
    p = LIT / name
    p.mkdir(parents=True, exist_ok=True)
    return p

def move_loose(fname, dest_folder):
    src = LIT / fname
    if src.exists():
        dst = dest_folder / fname
        shutil.move(str(src), str(dst))
        print(f"  MOVED  {fname} → {dest_folder.name}/")
        return dst
    return None

def download_pdf(url, dest_path, label=""):
    """Download url to dest_path; return True on success."""
    try:
        r = requests.get(url, headers=HEADS, allow_redirects=True, timeout=40, stream=True)
        ctype = r.headers.get("Content-Type","")
        if r.status_code == 200 and ("pdf" in ctype or r.content[:4] == b"%PDF"):
            dest_path.write_bytes(r.content)
            kb = len(r.content)//1024
            print(f"  DL OK  {label} → {dest_path.name}  [{kb} KB]")
            return True
        else:
            print(f"  DL FAIL {label} – HTTP {r.status_code}  ctype={ctype[:60]}")
            return False
    except Exception as e:
        print(f"  DL ERR  {label} – {e}")
        return False

def try_pmc_pdf(pmcid, dest_folder, stem):
    """Try several PMC/EuropePMC URL patterns."""
    urls = [
        f"https://europepmc.org/backend/ptpmcrender.fcgi?accid={pmcid}&blobtype=pdf",
        f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/pdf/",
        f"https://europepmc.org/articles/{pmcid}/pdf/main.pdf",
    ]
    for url in urls:
        dest = dest_folder / f"{stem}.pdf"
        if download_pdf(url, dest, f"{pmcid}"):
            return dest
        time.sleep(1.5)
    return None

def try_pmc_supplement(pmcid, dest_folder, stem):
    """Fetch Europe PMC article page and extract supplement file links."""
    api = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/supplementaryFiles?format=json"
    try:
        r = requests.get(api, timeout=20)
        if r.status_code != 200:
            return []
        data = r.json()
        files = data.get("supplementaryFiles", []) or data.get("files", [])
        downloaded = []
        for f in files:
            href = f.get("url") or f.get("href") or f.get("link")
            fname = f.get("filename") or f.get("name") or Path(href).name
            if not href:
                continue
            dest = dest_folder / f"supp_{fname}"
            if download_pdf(href, dest, f"supp/{fname}"):
                downloaded.append(dest)
            time.sleep(1)
        return downloaded
    except Exception as e:
        print(f"  SUPP ERR {pmcid} – {e}")
        return []

def try_elife(elife_url, dest_folder, stem):
    """Download eLife main PDF and attempt to list supplements via JSON API."""
    article_id = elife_url.rstrip("/").split("/")[-1]
    pdf_url = f"https://elifesciences.org/articles/{article_id}.pdf"
    dest = dest_folder / f"{stem}.pdf"
    ok = download_pdf(pdf_url, dest, f"eLife/{article_id}")
    if not ok:
        return
    # supplements via eLife's public API
    try:
        api = f"https://api.elifesciences.org/articles/{article_id}"
        r = requests.get(api, headers={"Accept":"application/vnd.elife.article-vor+json;version=1"}, timeout=20)
        if r.status_code == 200:
            data = r.json()
            supp_files = []
            # walk 'additionalFiles' section
            for block in data.get("body",[]) + data.get("additionalFiles",[]):
                if isinstance(block,dict) and block.get("type") == "file":
                    u = block.get("uri") or block.get("url")
                    n = block.get("title","supp")
                    if u:
                        supp_files.append((u, n))
            for su, sn in supp_files[:10]:
                dest_s = dest_folder / f"supp_{Path(su).name}"
                download_pdf(su, dest_s, f"eLife supp/{sn}")
                time.sleep(1)
    except Exception as e:
        print(f"  eLife supp API err: {e}")

def lookup_woo2017_pmc():
    """Try to get PMC ID for Woo 2017 Nat Commun (doi 10.1038/ncomms14211)."""
    url = ("https://www.ebi.ac.uk/europepmc/webservices/rest/search?"
           "query=DOI:10.1038/ncomms14211&format=json&resultType=lite")
    try:
        r = requests.get(url, timeout=15)
        results = r.json().get("resultList",{}).get("result",[])
        for res in results:
            pmcid = res.get("pmcid") or res.get("pmcId")
            if pmcid:
                return pmcid
    except:
        pass
    return None

# ── main ─────────────────────────────────────────────────────────────────────

def main():
    results = {}  # folder -> status

    # Pre-resolve Woo 2017 PMC
    print("Looking up Woo 2017 PMC ID ...")
    woo17_pmc = lookup_woo2017_pmc()
    print(f"  Woo 2017 PMC: {woo17_pmc or 'not found'}")
    # patch into STUDIES
    patched = []
    for s in STUDIES:
        if s[0] == "Woo_2017_expectancy" and woo17_pmc:
            s = (s[0], woo17_pmc, s[2], s[3], s[4])
        patched.append(s)

    # Build a folder→file map for already-downloaded PDFs (for shared-paper copies)
    folder_pdf = {}

    for (folder, pmcid, extra_url, loose_file, shared_with) in patched:
        print(f"\n── {folder}")
        fdir = make_folder(folder)
        status = "no source"

        # 1. Move existing loose PDF
        if loose_file:
            moved = move_loose(loose_file, fdir)
            if moved:
                folder_pdf[folder] = moved
                status = "existing PDF moved"

        # 2. Copy from shared-paper folder if already downloaded
        if shared_with and folder not in folder_pdf:
            src_pdf = folder_pdf.get(shared_with)
            if src_pdf and src_pdf.exists():
                dst = fdir / src_pdf.name
                shutil.copy2(str(src_pdf), str(dst))
                folder_pdf[folder] = dst
                print(f"  COPY   from {shared_with}/")
                status = "copied (same paper)"
            else:
                status = f"shared with {shared_with} (not yet downloaded)"

        # 3. Download from eLife
        elif extra_url and "elifesciences" in extra_url and folder not in folder_pdf:
            try_elife(extra_url, fdir, folder)
            pdfs = list(fdir.glob("*.pdf"))
            if pdfs:
                folder_pdf[folder] = pdfs[0]
                status = "downloaded (eLife)"
            else:
                status = "eLife download failed"

        # 4. Download from PMC
        elif pmcid and folder not in folder_pdf:
            dl = try_pmc_pdf(pmcid, fdir, folder)
            if dl:
                folder_pdf[folder] = dl
                status = "downloaded (PMC)"
                # also try supplements
                print(f"  Trying PMC supplements for {pmcid} ...")
                supps = try_pmc_supplement(pmcid, fdir, folder)
                if supps:
                    print(f"    → {len(supps)} supplement file(s) saved")
                time.sleep(2)
            else:
                status = "PMC download failed"

        # 5. Already had a PDF from step 1
        elif folder in folder_pdf:
            # try supplements if PMC known
            if pmcid:
                print(f"  Trying PMC supplements for {pmcid} ...")
                supps = try_pmc_supplement(pmcid, fdir, folder)
                if supps:
                    print(f"    → {len(supps)} supplement file(s) saved")
                time.sleep(1)

        results[folder] = status
        time.sleep(0.5)

    # ── summary report ───────────────────────────────────────────────────────
    print("\n\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    ok  = [(f,s) for f,s in results.items() if "fail" not in s and "no source" not in s and "shared" not in s.split()[0].lower()]
    bad = [(f,s) for f,s in results.items() if "fail" in s or "no source" in s]
    shr = [(f,s) for f,s in results.items() if "shared" in s and "not yet" in s]

    print(f"\n{'✓':2} GOT PDF  ({len(ok)}):")
    for f,s in ok:
        print(f"  {f[:55]:<55} {s}")
    print(f"\n{'~':2} NEEDS MANUAL DOWNLOAD ({len(bad)+len(shr)}):")
    for f,s in bad+shr:
        print(f"  {f[:55]:<55} {s}")

    # write STATUS.txt into Literature/
    with open(LIT/"STATUS.txt","w",encoding="utf-8") as fh:
        fh.write("Literature folder - download status\n")
        fh.write("="*60+"\n\n")
        for f,s in sorted(results.items()):
            fh.write(f"{f:<55} {s}\n")
    print(f"\nStatus log → Literature/STATUS.txt")

if __name__ == "__main__":
    main()
