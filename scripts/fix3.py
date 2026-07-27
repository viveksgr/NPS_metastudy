# -*- coding: utf-8 -*-
import io, sys, openpyxl
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
FN="Study_summary_completed.xlsx"
wb=openpyxl.load_workbook(FN); ws=wb["Studies"]
idx={ws.cell(r,3).value:r for r in range(2,ws.max_row+1)}
r=idx["Jepma_2018_expectancy_study1"]
ws.cell(r,18).value=None
ws.cell(r,19).value=None
wb.save(FN)
wb2=openpyxl.load_workbook(FN, read_only=True); s=wb2["Studies"]
for row in s.iter_rows(min_row=2, values_only=True):
    if str(row[2])=="Jepma_2018_expectancy_study1":
        print(f"Jepma1 -> class={row[17]!r} cal={row[18]!r}")
