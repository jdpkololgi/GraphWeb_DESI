# GraphWeb Plot Style Guide

Canonical visual style for every figure produced in this repo (and in `Illustris`,
where applicable) — matplotlib, seaborn, or anything built on top of them, in
both `.py` scripts and `.ipynb` notebooks.

**For agents (Claude Code, Cursor, etc.):** import `shared/plot_style.py` and
call `apply_style()` once per session/notebook, then use `finalize_axes(...)`
on every axes object before saving or displaying. Do not hand-roll colors,
fonts, or rcParams that duplicate or contradict this file — extend
`shared/plot_style.py` instead and update this doc in the same change.

This style derives from the FLATS conference deck (Eurostile, true-black,
magenta/blue/red accents) and the canonical four-class cosmic-web palette
confirmed against the deck on 2026-06-17.

---

## 1. Color system

### 1.1 Canvas

| Role | Value | rcParam(s) |
|---|---|---|
| Background (figure + axes) | `#000000` (true black) | `figure.facecolor`, `axes.facecolor`, `savefig.facecolor` |
| Primary text / ticks / axis labels | `#F2F2F2` (off-white) | `text.color`, `axes.labelcolor`, `xtick.color`, `ytick.color` |
| Spines (left/bottom only — top/right hidden) | `#F2F2F2` at reduced weight | `axes.edgecolor` |
| Gridlines (off by default; opt in per-plot) | `#F2F2F2` at `alpha=0.15` | set manually via `ax.grid(...)` |

Off-white rather than pure `#FFFFFF` is a deliberate choice — pure white on
true black causes visible halation/glow on projectors and is harsher to read
for a 15-minute talk. If you disagree, it's one constant
(`plot_style.TEXT_COLOR`) — change it there, not per-plot.

### 1.2 Cosmic-web environment colors (categorical — confirmed from FLATS deck)

Ordered low → high density. **This is a qualitative/categorical palette, not
a perceptually continuous one** — do not interpolate between these colors for
continuous quantities (density fields, posteriors, R² heatmaps). Use them
only for discrete class identity: scatter-by-predicted-class, confusion
matrices, class-fraction bar charts, legend swatches.

| Class | Hex | Swatch |
|---|---|---|
| Void | `#A1FCDD` | mint |
| Wall | `#4E84F7` | blue |
| Filament | `#EB336F` | magenta-pink |
| Cluster | `#F5C144` | gold |

Always use this exact order (Void, Wall, Filament, Cluster) for axis ticks,
legend entries, and confusion-matrix rows/columns — it matches the deck and
the established density ordering, and keeping it fixed avoids silently
reordered legends between plots.

### 1.3 General accent colors (non-categorical — deck UI accents)

For single-series highlights that are *not* encoding cosmic-web class — e.g.
an R² curve, a loss curve, a single annotation callout, a TARP diagonal
reference line.

| Name | Hex |
|---|---|
| Magenta | `#FF006E` |
| Blue | `#3A86FF` |
| Red | `#D62828` |

Do not reuse these for class encoding — that's what §1.2 is for. Mixing the
two palettes in one categorical legend is the most common way this style
guide gets violated; keep them in separate semantic roles.

### 1.4 Continuous data (colormaps)

Neither palette above is perceptually uniform, so neither is appropriate for
heatmaps, density fields, or posterior surfaces. Use:

- **Sequential** (density fields, posterior density, attention weights):
  matplotlib `magma` or `viridis` — both render well on true black.
- **Diverging** (residuals, predicted-minus-true): matplotlib `coolwarm` or
  `RdBu_r`, anchored at zero.

Always add a labeled colorbar when a colormap encodes a continuous quantity
— an unlabeled colorbar is treated as a style-guide violation (§3).

---

## 2. Typography

### 2.1 Font: IBM Plex Sans (everywhere — titles, labels, ticks, legends, body)

Eurostile (the FLATS deck title font) has no usable Greek glyph coverage,
which matters here because eigenvalue notation (λ₁, λ₂, λ₃) is central to
every other plot in this project. **IBM Plex Sans** was selected as the
replacement because it has genuine native monotonic Greek support (since
v3.0), is fully open-source (OFL), and has a clean, engineered, geometric
character that's the closest widely-available cousin to Eurostile's
technical feel without being a literal knockoff.

Use one family throughout, per the decision to keep this simple — don't
introduce a second body font.

| Weight | Use |
|---|---|
| IBM Plex Sans Bold | Titles |
| IBM Plex Sans Regular | Axis labels, tick labels, legend text, annotations |
| IBM Plex Sans Italic | Captions, secondary annotations (sparing use) |

### 2.2 Greek letters and math notation — use mathtext, not literal Unicode

Matplotlib renders anything inside `$...$` (e.g. `r'$\lambda_1$'`,
`r'$\lambda_2 = \lambda_1 + \mathrm{softplus}(\cdot)$'`) through its own
internal **mathtext** engine, which is independent of whatever font is set
for regular text. This means eigenvalue notation renders correctly
regardless of the body font's own glyph coverage — IBM Plex Sans's native
Greek glyphs are a fallback for plain (non-mathtext) Greek text only (e.g. a
literal "λ" typed directly into a legend string outside `$...$`), not the
primary mechanism.

`plot_style.apply_style()` sets `mathtext.fontset = 'dejavusans'` (bundled
with matplotlib, zero extra install, visually compatible sans-serif) so this
works out of the box on any machine, including bare NERSC nodes with no
fonts installed.

**Rule: always write eigenvalues and other Greek-letter quantities as
mathtext** (`r'$\lambda_1$'`, not `'λ1'` or `'lambda_1'`) for both correctness
and visual consistency.

### 2.3 Font installation (one-time, per machine)

IBM Plex Sans is not a default system font on NERSC login/compute nodes or
guaranteed on a fresh Cursor/Mac setup, so it's bundled in-repo rather than
assumed installed — same pattern as the model/scaler paths in
`shared/config_paths.py`.

```bash
mkdir -p assets/fonts
curl -L -o assets/fonts/IBMPlexSans-Regular.ttf \
  "https://raw.githubusercontent.com/IBM/plex/master/packages/plex-sans/fonts/complete/ttf/IBMPlexSans-Regular.ttf"
curl -L -o assets/fonts/IBMPlexSans-Bold.ttf \
  "https://raw.githubusercontent.com/IBM/plex/master/packages/plex-sans/fonts/complete/ttf/IBMPlexSans-Bold.ttf"
curl -L -o assets/fonts/IBMPlexSans-Italic.ttf \
  "https://raw.githubusercontent.com/IBM/plex/master/packages/plex-sans/fonts/complete/ttf/IBMPlexSans-Italic.ttf"
```

`plot_style.register_fonts()` (called automatically by `apply_style()`) loads
these from `assets/fonts/` if present, and falls back to DejaVu Sans with a
single `warnings.warn(...)` if they're missing — it will never hard-crash a
notebook over a missing font file, but the fallback is visually obvious
(different typeface) so it won't go unnoticed.

---

## 3. Mandatory elements on every figure

Every figure produced in this repo must have all of the following. An
agent generating or editing plotting code should treat a missing item below
as a bug, not a style nitpick:

1. **A descriptive title** — what the plot shows, not just a variable name
   (`"Eigenvalue regression: predicted vs. true λ₂"`, not `"Figure 3"`).
2. **Axis labels with units** where units exist (`"Comoving distance [Mpc/h]"`,
   not `"distance"`).
3. **A legend** whenever more than one series, class, or colored element
   appears — including confusion-matrix-style categorical plots, where the
   class color key counts as the legend.
4. **A labeled colorbar** whenever a colormap encodes a continuous quantity
   (§1.4).

### Size hierarchy (set globally by `apply_style()`, don't override per-plot
unless there's a specific reason)

| Element | Size | Weight |
|---|---|---|
| Title | 18pt | Bold |
| Axis labels | 14pt | Regular |
| Tick labels | 12pt | Regular |
| Legend text | 12pt | Regular |

### Figure defaults

- Default figure size: `(8, 6)` inches.
- Screen/notebook DPI: 150. Saved-figure DPI: 300.
- Top and right spines hidden; left and bottom kept at `TEXT_COLOR`.
- Gridlines off by default; if a plot benefits from them (e.g. TARP coverage
  diagonal), turn on explicitly with `alpha=0.15`.
- Legend: filled box (`facecolor=BACKGROUND`, `framealpha=0.85`,
  `edgecolor=GRID_COLOR`) rather than frameless — keeps it legible over busy
  scatter/heatmap content.
- Always save with `facecolor=fig.get_facecolor()` (or just rely on
  `savefig.facecolor` rcParam) so the true-black background survives export
  — the matplotlib default is to save figures on a white background, which
  silently breaks this entire style guide if not set.

---

## 4. Usage

```python
from shared.plot_style import apply_style, finalize_axes, COSMIC_WEB_COLORS, CLASS_ORDER
import matplotlib.pyplot as plt

apply_style()  # call once per script / notebook session

fig, ax = plt.subplots()
for cls in CLASS_ORDER:
    mask = labels == cls
    ax.scatter(x[mask], y[mask], color=COSMIC_WEB_COLORS[cls], label=cls.capitalize(), s=8)

finalize_axes(
    ax,
    title=r"Predicted vs. true $\lambda_2$ by cosmic-web class",
    xlabel=r"True $\lambda_2$",
    ylabel=r"Predicted $\lambda_2$",
)
fig.savefig("figure.png")
```

`finalize_axes` sets the title/labels, calls `ax.legend()` if any labeled
artists exist, and warns (rather than failing) if a legend was expected but
no labeled artists were found — that warning is the signal to go add
`label=...` to whatever was plotted, not to ignore it.

---

## 5. Provenance

- Background, font choice, and rationale: science discussion, 2026-06-17.
- Cosmic-web class hex values: confirmed directly against the FLATS deck
  legend swatches (slide 9) by Dakshesh, 2026-06-17.
- General accent colors and Eurostile-as-title-font precedent: FLATS deck
  visual identity (see `Illustris/SCIENCE_LOG.md`, 2026-06-15 entry).
