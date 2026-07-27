# -*- coding: utf-8 -*-
import requests, re, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

r = requests.get("https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11389792/fullTextXML", timeout=30)
print("Status:", r.status_code, "Size:", len(r.content))
xml = r.text
doi_m = re.search(r'pub-id-type=["\']doi["\']>([^<]+)<', xml)
print("DOI:", doi_m.group(1) if doi_m else "not found")
print("First 500:", xml[:500])
