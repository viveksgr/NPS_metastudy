# -*- coding: utf-8 -*-
# Verify Kim 2024 DOI vs PMC ID discrepancy
import requests, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

# Check what adk7421 DOI resolves to
r = requests.get("https://doi.org/10.1126/sciadv.adk7421",
                 headers={"User-Agent":"Mozilla/5.0","Accept":"text/html"},
                 allow_redirects=True, timeout=20)
print("adk7421 resolved URL:", r.url, "status:", r.status_code)
import re
title = re.search(r"<title>([^<]+)", r.text)
print("  page title:", title.group(1)[:80] if title else "n/a")

# Check ado8230
r2 = requests.get("https://doi.org/10.1126/sciadv.ado8230",
                  headers={"User-Agent":"Mozilla/5.0","Accept":"text/html"},
                  allow_redirects=True, timeout=20)
print("ado8230 resolved URL:", r2.url)
title2 = re.search(r"<title>([^<]+)", r2.text)
print("  page title:", title2.group(1)[:80] if title2 else "n/a")
