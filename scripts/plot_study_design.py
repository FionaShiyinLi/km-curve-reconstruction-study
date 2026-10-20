#!/usr/bin/env python3
"""Assemble Figure 1 from actual saved native-PDF inputs and dataset assignments.
Inputs: data/run_manifest.csv, data/input_figures/; outputs: build/figures/Figure1_*.
Run: python scripts/plot_study_design.py (repository root).
The original vector crop/placement code is retained. No curves are simulated.
"""
from collections import Counter
from copy import deepcopy
import csv
import hashlib
from io import BytesIO
import json
from pathlib import Path
from matplotlib import get_data_path
from pypdf import PdfReader, PdfWriter, Transformation
from pypdf.generic import RectangleObject
from reportlab.lib.colors import HexColor
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen import canvas
ROOT=Path(__file__).resolve().parents[1]
FIGURES=ROOT/'build/figures'
SCENARIOS=('S1','S2','S3','S4')

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def prepare():
    FIGURES.mkdir(parents=True,exist_ok=True)
    with (ROOT/'data/run_manifest.csv').open() as f:runs=list(csv.DictReader(f))
    counts={};cases=[];sources={}
    for scenario in SCENARIOS:
        part=[r for r in runs if r['scenario']==scenario]
        primary={r['dataset_id'] for r in part if r['primary64']=='True'}
        repeated={r['dataset_id'] for r in part if r['repeatability']=='True'}
        followup={r['dataset_id'] for r in part if r['experiment']=='followup'}
        alternate={r['dataset_id'] for r in part if r['experiment']=='survminer'}
        assert repeated<=primary and alternate<=followup<=primary
        counts[scenario]=dict(primary=len(primary),repeat=len(repeated),followup=len(followup),survminer=len(alternate))
        dataset=min(alternate)
        pair={}
        for name,experiment in [('original','primary'),('alternative','survminer')]:
            row=next(r for r in part if r['dataset_id']==dataset and r['experiment']==experiment and r['condition']=='PDF' and r['complete']=='True' and r['is_first_run']=='True')
            pair[name]=row['input_figure'];sources[pair[name]]=sha(ROOT/pair[name])
        n=int(row['n_per_arm'])
        cases.append(dict(scenario=scenario,dataset=dataset,n_per_arm=n,hazard_pattern='Proportional' if scenario in ('S1','S3') else 'Crossing',inputs=pair))
    font_dir = Path(get_data_path()) / "fonts/ttf"
    pdfmetrics.registerFont(TTFont("FigureSans", str(font_dir / "DejaVuSans.ttf")))
    pdfmetrics.registerFont(TTFont("FigureSans-Bold", str(font_dir / "DejaVuSans-Bold.ttf")))
    width, height = 468, 564
    plot_width, plot_height = 194, 66
    columns = {"original": 38, "alternative": 267}
    buffer = BytesIO()
    c = canvas.Canvas(buffer, pagesize=(width, height), invariant=1)
    c.setTitle("Simulation scenarios, plotting templates, and dataset subsets")
    for name, title, colors in (("original", "Original template (ReportLab)", ("#1f77b4", "#ff7f0e")),
                                 ("alternative", "R template (survminer)", ("#F8766D", "#00BFC4"))):
        center = columns[name] + plot_width / 2
        c.setFillColor(HexColor("#202020")); c.setFont("FigureSans-Bold", 9)
        c.drawCentredString(center, 551, title)
        for offset, group, color in ((-64, "Group A", colors[0]), (16, "Group B", colors[1])):
            x = center + offset
            c.setStrokeColor(HexColor(color)); c.setLineWidth(1.1)
            c.line(x, 538, x + 15, 538); c.line(x + 8, 535, x + 8, 541)
            c.setFillColor(HexColor("#333333")); c.setFont("FigureSans", 8)
            c.drawString(x + 20, 535, group)

    # Native-PDF crop boxes, in bottom-left page coordinates. ReportLab's
    # panel spans the data limits; ggplot2's panel includes its 5% expansion.
    crops = {"original": (64.3, 100.3, 410.9, 281.3),
             "alternative": (57.152, 117.390, 426.519, 279.324)}
    placements = []
    for index, case in enumerate(cases):
        title_y, bottom = 520 - index * 100, 441 - index * 100
        c.setFillColor(HexColor("#202020")); c.setFont("FigureSans-Bold", 9)
        c.drawString(27, title_y, f"{case['scenario']}   {case['n_per_arm']:,} per arm  |  {case['hazard_pattern']} hazards")
        c.setFont("FigureSans", 7.5)
        c.drawRightString(461, title_y, f"Example: {case['dataset'].replace('PP32_', '').replace('_', ' ')}")
        for name in ("original", "alternative"):
            x = columns[name]
            c.setFillColor(HexColor("#202020")); c.setFont("FigureSans", 7.5)
            for t in (0, 12, 24, 36):
                fraction = ((64.8 + 345.6 * t / 36 - crops[name][0]) / (crops[name][2] - crops[name][0])
                            if name == "original" else (t / 36 + 0.05) / 1.1)
                xx = x + plot_width * fraction
                c.setStrokeColor(HexColor("#333333")); c.setLineWidth(0.4)
                c.line(xx, bottom, xx, bottom - 2.5)
                c.drawCentredString(xx, bottom - 10, str(t))
            for s in (0, 0.5, 1):
                fraction = ((100.8 + 180 * s - crops[name][1]) / (crops[name][3] - crops[name][1])
                            if name == "original" else (s + 0.05) / 1.1)
                yy = bottom + plot_height * fraction
                c.line(x - 2.5, yy, x, yy)
                c.drawRightString(x - 5, yy - 2.4, f"{s:.1f}")
            placements.append({"scenario": case["scenario"], "template": name,
                               "source": case["inputs"][name], "crop_points": crops[name],
                               "destination_points": (x, bottom, plot_width, plot_height)})
    c.saveState(); c.translate(10, 330); c.rotate(90); c.setFont("FigureSans", 8)
    c.drawCentredString(0, 0, "Survival probability"); c.restoreState()
    c.setFillColor(HexColor("#202020")); c.setFont("FigureSans", 8)
    for x in columns.values(): c.drawCentredString(x + plot_width / 2, 116, "Time (months)")

    c.setFont("FigureSans-Bold", 9); c.drawString(14, 101, "Dataset allocation")
    edges = (14, 135, 199, 255, 311, 379, 461)
    headings = (("Scenario", ""), ("Participants", "per arm"), ("Primary", "datasets"),
                ("Repeat", "datasets"), ("Follow-up", "datasets"), ("Survminer", "datasets"))
    c.setStrokeColor(HexColor("#333333")); c.setLineWidth(0.6); c.line(14, 95, 461, 95)
    c.setFont("FigureSans-Bold", 7)
    for i, lines in enumerate(headings):
        center = (edges[i] + edges[i + 1]) / 2
        c.drawCentredString(center, 84, lines[0]); c.drawCentredString(center, 75, lines[1])
    c.setLineWidth(0.4); c.line(14, 71, 461, 71)
    table_rows = []
    for i, case in enumerate(cases):
        scenario = case["scenario"]
        values = [f"{scenario}: {case['hazard_pattern']}", f"{case['n_per_arm']:,}"]
        values += [str(counts[scenario][key]) for key in ("primary", "repeat", "followup", "survminer")]
        table_rows.append(values)
    totals = {key: sum(row[key] for row in counts.values()) for key in ("primary", "repeat", "followup", "survminer")}
    assert totals == {"primary": 64, "repeat": 24, "followup": 32, "survminer": 16}
    table_rows.append(["Total", ""] + [str(v) for v in totals.values()])
    for i, values in enumerate(table_rows):
        y = 60 - i * 12
        if i == 4: c.line(14, y + 9, 461, y + 9)
        c.setFont("FigureSans-Bold" if i == 4 else "FigureSans", 7.5)
        for j, value in enumerate(values):
            if j == 0: c.drawString(edges[j] + 3, y, value)
            else: c.drawCentredString((edges[j] + edges[j + 1]) / 2, y, value)
    c.setLineWidth(0.6); c.line(14, 6, 461, 6); c.showPage(); c.save()

    page = PdfReader(BytesIO(buffer.getvalue())).pages[0]
    for placement in placements:
        source = deepcopy(PdfReader(ROOT / placement["source"]).pages[0])
        x0, y0, x1, y1 = placement["crop_points"]
        x, y, w, h = placement["destination_points"]
        source.cropbox = RectangleObject((x0, y0, x1, y1))
        transform = Transformation().translate(-x0, -y0).scale(w / (x1 - x0), h / (y1 - y0)).translate(x, y)
        page.merge_transformed_page(source, transform, expand=False)
    writer = PdfWriter(); writer.add_page(page)
    writer.add_metadata({"/Title": "Simulation scenarios, plotting templates, and dataset subsets"})
    target = FIGURES / "Figure1_study_design.pdf"
    with target.open("wb") as handle: writer.write(handle)
    with (FIGURES / "Figure1_dataset_subsets.csv").open("w", newline="") as handle:
        out = csv.writer(handle); out.writerow(["scenario", "participants_per_arm", "primary", "repeat", "followup", "survminer"])
        out.writerows(table_rows)
    audit = {"selection_rule": "Lowest dataset ID in each scenario among the 16 datasets with both plotting templates; no accuracy or appearance selection",
             "cases": cases, "counts": counts, "totals": totals,
             "display": "Cropped vector curve regions from actual inputs, resized with common external tick labels and legends. Risk tables omitted from thumbnails; no curve or censoring data regenerated.",
             "placements": placements, "sources_sha256": sources,
             "outputs_sha256": {name: sha(FIGURES / name) for name in ("Figure1_study_design.pdf", "Figure1_dataset_subsets.csv")}}
    (FIGURES / "Figure1_study_design.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps({"cases": [c["dataset"] for c in cases], "totals": totals, "source_checks": "passed"}, indent=2))

if __name__=="__main__":prepare()
