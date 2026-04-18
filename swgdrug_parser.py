"""
SWGDRUG Library Parser
======================
Parses the SWGDRUG .txt library, exports all compounds to CSV,
and generates a styled mass spectrum plot for each compound.

Requirements:
    pip install matplotlib pandas

Usage:
    python swgdrug_parser.py
"""

import re
import csv
import os
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import FancyArrowPatch
import matplotlib.patheffects as pe
import pandas as pd

# ─── CONFIG ────────────────────────────────────────────────────────────────────
INPUT_FILE  = r"C:\Users\txbar\Downloads\SWGDRUG 314.txt"
OUTPUT_DIR  = r"C:\Users\txbar\Downloads"
CSV_OUTPUT  = os.path.join(OUTPUT_DIR, "SWGDRUG_compounds.csv")

# Plot style constants
BG_COLOR    = "#0d1117"
PANEL_COLOR = "#161b22"
SPINE_COLOR = "#30363d"
TICK_COLOR  = "#8b949e"
LABEL_COLOR = "#e6edf3"
GRID_COLOR  = "#21262d"
BAR_COLOR   = "#58a6ff"
BASE_COLOR  = "#388bfd"
ANNOT_COLOR = "#ffa657"
TITLE_COLOR = "#f0f6fc"
SUB_COLOR   = "#8b949e"

TOP_N_LABEL = 10       # how many peaks to annotate with m/z labels
DPI         = 150
FIG_W, FIG_H = 14, 6
# ───────────────────────────────────────────────────────────────────────────────


def sanitize_filename(name: str) -> str:
    """Remove characters that are invalid in Windows filenames."""
    return re.sub(r'[\\/*?:"<>|]', "_", name).strip()


def parse_swgdrug(filepath: str) -> list[dict]:
    """
    Parse the SWGDRUG text file into a list of compound dicts.
    Each dict has metadata keys + a 'peaks' list of (mz, rel_ab) tuples.
    """
    compounds = []
    current   = {}
    in_peaks  = False
    peak_buf  = []
    num_peaks = 0

    meta_keys = {"Name", "Formula", "MW", "ExactMass", "CASNO", "ID", "Comment"}

    with open(filepath, "r", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.rstrip()

            # ── start of a new compound ──────────────────────────────────────
            if line.startswith("Name:"):
                if current:                        # save previous
                    current["peaks"] = peak_buf
                    compounds.append(current)
                current   = {}
                in_peaks  = False
                peak_buf  = []
                num_peaks = 0
                current["Name"] = line.split(":", 1)[1].strip()
                continue

            # ── metadata key:value lines ─────────────────────────────────────
            matched = False
            for key in meta_keys - {"Name"}:
                if line.startswith(f"{key}:"):
                    current[key] = line.split(":", 1)[1].strip()
                    matched = True
                    break

            if matched:
                continue

            # ── Num peaks triggers peak-reading mode ─────────────────────────
            if line.startswith("Num peaks:"):
                num_peaks = int(line.split(":", 1)[1].strip())
                in_peaks  = True
                continue

            # ── peak data (space-separated pairs on one or more lines) ───────
            if in_peaks and line.strip():
                tokens = line.split()
                # each pair is mz rel_ab
                for i in range(0, len(tokens) - 1, 2):
                    try:
                        mz    = int(tokens[i])
                        rel_ab = float(tokens[i + 1])
                        peak_buf.append((mz, rel_ab))
                    except ValueError:
                        pass

    # save last compound
    if current:
        current["peaks"] = peak_buf
        compounds.append(current)

    return compounds


def write_csv(compounds: list[dict], csv_path: str) -> None:
    """Write all compound metadata (no peaks) to a CSV file."""
    meta_fields = ["Name", "Formula", "MW", "ExactMass", "CASNO", "ID", "Comment", "Num_Peaks"]
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=meta_fields, extrasaction="ignore")
        writer.writeheader()
        for c in compounds:
            row = {k: c.get(k, "") for k in meta_fields}
            row["Num_Peaks"] = len(c.get("peaks", []))
            writer.writerow(row)
    print(f"  ✓ CSV written → {csv_path}  ({len(compounds)} compounds)")


def plot_spectrum(compound: dict, output_dir: str) -> None:
    """Render a dark-themed mass spectrum bar plot and save as PNG."""
    name   = compound.get("Name", "Unknown")
    formula = compound.get("Formula", "")
    mw     = compound.get("MW", "")
    casno  = compound.get("CASNO", "")
    peaks  = compound.get("peaks", [])

    if not peaks:
        print(f"  ⚠  No peaks for {name!r} — skipped")
        return

    mz_vals  = [p[0] for p in peaks]
    rel_vals = [p[1] for p in peaks]
    max_ab   = max(rel_vals) if rel_vals else 1

    # normalise to 0-999 (already in that scale, but just in case)
    rel_norm = [v / max_ab * 999 for v in rel_vals]

    # ── figure ────────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H), facecolor=BG_COLOR)
    ax.set_facecolor(PANEL_COLOR)

    # subtle grid
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, color=GRID_COLOR, linewidth=0.6, linestyle="--", alpha=0.7)
    ax.xaxis.grid(False)

    # ── bars ──────────────────────────────────────────────────────────────────
    # colour-map: taller bars get lighter blue tones
    norm_frac = [v / 999 for v in rel_norm]
    colors = [
        plt.cm.Blues(0.35 + 0.65 * f)     # range from pale to vivid blue
        for f in norm_frac
    ]

    ax.bar(mz_vals, rel_norm, width=1.2, color=colors,
           linewidth=0, zorder=3, snap=False)

    # ── base line ─────────────────────────────────────────────────────────────
    ax.axhline(0, color=SPINE_COLOR, linewidth=1)

    # ── annotations: top N peaks ──────────────────────────────────────────────
    # pick the top-N by rel_ab
    sorted_idx = sorted(range(len(rel_norm)), key=lambda i: rel_norm[i], reverse=True)
    label_set  = set(sorted_idx[:TOP_N_LABEL])

    for i in label_set:
        mz  = mz_vals[i]
        rab = rel_norm[i]
        ax.text(
            mz, rab + 18, str(mz_vals[i]),
            ha="center", va="bottom",
            fontsize=7.5, fontweight="bold",
            color=ANNOT_COLOR,
            path_effects=[pe.withStroke(linewidth=1.5, foreground=PANEL_COLOR)],
            clip_on=True,
        )

    # ── axes ──────────────────────────────────────────────────────────────────
    ax.set_xlim(max(0, min(mz_vals) - 5), max(mz_vals) + 10)
    ax.set_ylim(-30, 1100)

    ax.set_xlabel("m/z", color=LABEL_COLOR, fontsize=12, labelpad=8,
                  fontfamily="monospace")
    ax.set_ylabel("Relative Abundance", color=LABEL_COLOR, fontsize=12, labelpad=8,
                  fontfamily="monospace")

    ax.tick_params(colors=TICK_COLOR, labelsize=9, length=4, width=0.8)
    for spine in ax.spines.values():
        spine.set_color(SPINE_COLOR)
        spine.set_linewidth(0.8)

    ax.yaxis.set_major_formatter(ticker.FuncFormatter(
        lambda x, _: f"{int(x)}" if x >= 0 else ""))
    ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True, nbins=16))

    # ── title block ───────────────────────────────────────────────────────────
    subtitle_parts = []
    if formula: subtitle_parts.append(f"Formula: {formula}")
    if mw:      subtitle_parts.append(f"MW: {mw}")
    if casno:   subtitle_parts.append(f"CAS: {casno}")
    subtitle_parts.append(f"Peaks: {len(peaks)}")
    subtitle = "   |   ".join(subtitle_parts)

    fig.text(0.5, 0.97, name,
             ha="center", va="top",
             fontsize=16, fontweight="bold",
             color=TITLE_COLOR, fontfamily="monospace")
    fig.text(0.5, 0.925, subtitle,
             ha="center", va="top",
             fontsize=9, color=SUB_COLOR, fontfamily="monospace")

    # ── SWGDRUG watermark ─────────────────────────────────────────────────────
    ax.text(0.98, 0.97, "SWGDRUG Library",
            transform=ax.transAxes,
            ha="right", va="top",
            fontsize=8, color=SPINE_COLOR,
            fontstyle="italic", fontfamily="monospace")

    # ── base intensity reference line at 999 ──────────────────────────────────
    ax.axhline(999, color=SPINE_COLOR, linewidth=0.5, linestyle=":", alpha=0.5)
    ax.text(ax.get_xlim()[0] + 1, 1003, "999",
            color=TICK_COLOR, fontsize=7, va="bottom", fontfamily="monospace")

    # ── tight layout & save ───────────────────────────────────────────────────
    fig.subplots_adjust(top=0.86, bottom=0.12, left=0.07, right=0.97)

    safe_name = sanitize_filename(name)
    out_path  = os.path.join(output_dir, f"{safe_name}.png")
    fig.savefig(out_path, dpi=DPI, facecolor=BG_COLOR, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ {safe_name}.png")


# ── MAIN ───────────────────────────────────────────────────────────────────────
def main():
    print(f"\n{'─'*60}")
    print("  SWGDRUG Library Parser")
    print(f"{'─'*60}")
    print(f"  Input : {INPUT_FILE}")
    print(f"  Output: {OUTPUT_DIR}\n")

    if not os.path.isfile(INPUT_FILE):
        raise FileNotFoundError(f"Cannot find input file:\n  {INPUT_FILE}")

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("Parsing compounds …")
    compounds = parse_swgdrug(INPUT_FILE)
    print(f"  → {len(compounds)} compounds found\n")

    print("Writing CSV …")
    write_csv(compounds, CSV_OUTPUT)
    print()

    print(f"Generating {len(compounds)} mass spectrum plots …")
    for i, cmpd in enumerate(compounds, 1):
        print(f"  [{i}/{len(compounds)}] ", end="")
        plot_spectrum(cmpd, OUTPUT_DIR)

    print(f"\n{'─'*60}")
    print("  Done! All files saved to:")
    print(f"  {OUTPUT_DIR}")
    print(f"{'─'*60}\n")


if __name__ == "__main__":
    main()