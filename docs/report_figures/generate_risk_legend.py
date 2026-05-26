"""Generate a stand-alone legend figure showing the four-level risk
classification (Action / Review / Normal / No Data) used by both the
patient and clinician reports.

Produces: output/risk_legend.svg
"""

from pathlib import Path
from html import escape

OUT = Path(__file__).parent / "output" / "risk_legend.svg"

# ---- colors (mirroring src/reports/templates/base.html palette) ----
COL = {
    "action":     "#C62828",
    "action_bg":  "#FFEBEE",
    "review":     "#E65100",
    "review_bg":  "#FFF3E0",
    "normal":     "#2E7D32",
    "normal_bg":  "#E8F5E9",
    "nodata":     "#455A64",
    "nodata_bg":  "#ECEFF1",
    "text":       "#263238",
    "text_light": "#607D8B",
    "border":     "#CFD8DC",
    "title":      "#1A2730",
}
FONT = ("-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, "
        "Helvetica, Arial, sans-serif")

# ---- canvas ----
W = 900
PAD = 32
ROW_H = 92
ROW_GAP = 14
TITLE_BLOCK = 70   # space taken by title at top
N_ROWS = 4
H = TITLE_BLOCK + N_ROWS * ROW_H + (N_ROWS - 1) * ROW_GAP + PAD * 2 - PAD

# Inside each row:
ICON_BOX = 64           # square icon area on the left
LABEL_X = PAD + ICON_BOX + 18
DEF_X = LABEL_X + 130    # text column for the definition
DEF_W = W - DEF_X - PAD


def eye_svg(cx: int, cy: int, color: str) -> str:
    """Stylized eye glyph used for the Review category in the patient
    report. Centered at (cx, cy), sized about 28 px wide."""
    return (
        f'<g transform="translate({cx-14},{cy-10})">'
        f'<path d="M2 12s3.5-7 10-7 10 7 10 7-3.5 7-10 7-10-7-10-7z"'
        f' fill="none" stroke="{color}" stroke-width="2.2"'
        f' stroke-linecap="round" stroke-linejoin="round"/>'
        f'<circle cx="12" cy="12" r="3" fill="none" stroke="{color}"'
        f' stroke-width="2.2"/>'
        f'</g>'
    )


def category_row(y, label, accent, bg, definition, icon_kind):
    """One row: icon box on the left, label + definition on the right."""
    parts = []
    # row background — light tint of the accent color
    parts.append(
        f'<rect x="{PAD}" y="{y}" width="{W - 2*PAD}" height="{ROW_H}" '
        f'rx="8" fill="{bg}" stroke="{accent}" stroke-width="1.2"/>'
    )
    # icon box
    icon_x = PAD + 14
    icon_y = y + (ROW_H - ICON_BOX) / 2
    parts.append(
        f'<rect x="{icon_x}" y="{icon_y}" width="{ICON_BOX}" '
        f'height="{ICON_BOX}" rx="10" fill="{accent}"/>'
    )
    cx = icon_x + ICON_BOX / 2
    cy = icon_y + ICON_BOX / 2
    if icon_kind == "exclamation":
        parts.append(
            f'<text x="{cx}" y="{cy + 18}" text-anchor="middle" '
            f'font-family="{FONT}" font-size="44" font-weight="700" '
            f'fill="#FFFFFF">!</text>'
        )
    elif icon_kind == "eye":
        parts.append(eye_svg(int(cx), int(cy), "#FFFFFF"))
    elif icon_kind == "check":
        parts.append(
            f'<text x="{cx}" y="{cy + 14}" text-anchor="middle" '
            f'font-family="{FONT}" font-size="40" font-weight="700" '
            f'fill="#FFFFFF">✓</text>'
        )
    elif icon_kind == "question":
        parts.append(
            f'<text x="{cx}" y="{cy + 16}" text-anchor="middle" '
            f'font-family="{FONT}" font-size="40" font-weight="700" '
            f'fill="#FFFFFF">?</text>'
        )
    # category label
    parts.append(
        f'<text x="{LABEL_X}" y="{y + 36}" font-family="{FONT}" '
        f'font-size="22" font-weight="700" fill="{accent}">'
        f'{escape(label)}</text>'
    )
    # short definition (one or two lines)
    line_y = y + 36
    for i, line in enumerate(definition):
        parts.append(
            f'<text x="{DEF_X}" y="{line_y + i*20}" font-family="{FONT}" '
            f'font-size="14" fill="{COL["text"]}">'
            f'{escape(line)}</text>'
        )
    return "\n".join(parts)


def main():
    rows = [
        dict(label="Action", accent=COL["action"], bg=COL["action_bg"],
             icon_kind="exclamation",
             definition=[
                 "The patient's genotype calls for a clear change in",
                 "prescribing — different drug, dose, or contraindication.",
             ]),
        dict(label="Review", accent=COL["review"], bg=COL["review_bg"],
             icon_kind="eye",
             definition=[
                 "Standard prescribing may still apply, but the patient's",
                 "genotype warrants dose adjustment or extra monitoring.",
             ]),
        dict(label="Normal", accent=COL["normal"], bg=COL["normal_bg"],
             icon_kind="check",
             definition=[
                 "The patient's genotype produces a standard / expected",
                 "phenotype — no pharmacogenomic action is needed.",
             ]),
        dict(label="No Data", accent=COL["nodata"], bg=COL["nodata_bg"],
             icon_kind="question",
             definition=[
                 "The pipeline could not generate a result for this gene —",
                 "missing coverage or no actionable phenotype call.",
             ]),
    ]

    body = []
    # Title
    body.append(
        f'<text x="{PAD}" y="{PAD + 20}" font-family="{FONT}" '
        f'font-size="20" font-weight="700" fill="{COL["title"]}">'
        f'Four-level risk classification used by the PGx reports</text>'
    )
    body.append(
        f'<text x="{PAD}" y="{PAD + 44}" font-family="{FONT}" '
        f'font-size="12" fill="{COL["text_light"]}">'
        f'Applied uniformly to every gene and drug; drives section grouping '
        f'in both the patient and clinician reports.</text>'
    )

    y = TITLE_BLOCK + PAD
    for r in rows:
        body.append(category_row(y, **r))
        y += ROW_H + ROW_GAP

    svg = (
        f'<?xml version="1.0" encoding="UTF-8"?>\n'
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'width="{W}" height="{H}" viewBox="0 0 {W} {H}">\n'
        f'<rect width="100%" height="100%" fill="#FFFFFF"/>\n'
        + "\n".join(body) + "\n</svg>\n"
    )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(svg, encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
