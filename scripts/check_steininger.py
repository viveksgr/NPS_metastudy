# -*- coding: utf-8 -*-
import requests, re, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

r = requests.get("https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11906725/fullTextXML", timeout=30)
xml = r.text

# Find DOI
doi_m = re.search(r'pub-id-type=["\']doi["\']>([^<]+)<', xml)
print("DOI:", doi_m.group(1) if doi_m else "not found")

# Find any Nature article ID
nat_m = re.findall(r's41467-\d{3}-\d{5}-\w', xml)
print("Nature IDs:", nat_m[:5])

print("\nFirst 1500 chars of XML:")
print(xml[:1500])
