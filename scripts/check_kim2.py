# -*- coding: utf-8 -*-
import requests, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
H = {"User-Agent": "Mozilla/5.0 Chrome/124", "Accept": "application/pdf,*/*",
     "Referer": "https://www.science.org/"}
for doi_suf in ["ado8230", "adk7421"]:
    r = requests.get(f"https://www.science.org/doi/pdf/10.1126/sciadv.{doi_suf}",
                     headers=H, allow_redirects=True, timeout=30)
    ct = r.headers.get("Content-Type","")
    print(f"{doi_suf}: {r.status_code} {ct[:40]} first4={r.content[:4]}")
    if r.content[:4] == b"%PDF":
        fname = f"Literature/Kim_2024_cueIntegration/sciadv.{doi_suf}.pdf"
        open(fname,"wb").write(r.content)
        print(f"  SAVED {len(r.content)//1024} KB")
        break
