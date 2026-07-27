# -*- coding: utf-8 -*-
import io, sys, openpyxl
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
FN="Study_summary_completed.xlsx"
wb=openpyxl.load_workbook(FN); ws=wb["Studies"]
idx={ws.cell(r,3).value:r for r in range(2,ws.max_row+1)}

# Brascher calibration
r=idx["Brascher_2016_control"]; ws.cell(r,19,"Calibrated")

# LopezSola 2018: N was stuck in Population ("29 (women)"); repair
r=idx["LopezSola_2018_prosocialMeaning"]
ws.cell(r,12,"29"); ws.cell(r,13,"Healthy women")

# Jepma study1 intentionally blank -> clear moderator guesses
r=idx["Jepma_2018_expectancy_study1"]
ws.cell(r,18,None); ws.cell(r,19,None)

wb.save(FN)
print("fixed Brascher cal, LopezSola N, Jepma1 blanked")
# show the three
wb2=openpyxl.load_workbook(FN, read_only=True); s=wb2["Studies"]
for row in s.iter_rows(min_row=2, values_only=True):
    if str(row[2]) in ("Brascher_2016_control","LopezSola_2018_prosocialMeaning","Jepma_2018_expectancy_study1"):
        print(f"  {row[2][:34]:35} N={str(row[11])[:14]:15} Pop={str(row[12])[:14]:15} class={row[17]} cal={row[18]}")
