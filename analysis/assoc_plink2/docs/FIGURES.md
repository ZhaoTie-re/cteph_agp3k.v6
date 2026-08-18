# Figures — `assoc_plink2`

Conventions shared by every figure. Per-figure detail lives in the `.md` sitting beside each PNG,
written by `scripts/figure_doc.py`; this file is the layer above that — what is true of all of them.

**19 figures**: 9 scan (3 cohorts × 3 models), 3 × 3 per-peak figures for the three genome-wide peaks,
1 cross-cohort comparison. There is no separate peak figure — a scan's peaks are panels of its own
scan figure, on the Manhattan's axis.

---

## Where they are

```
results/<cohort>/figures/
  01.scan/         scan.additive.png   scan.dominant.png   scan.recessive.png
  02.regional/     regional.<peak_id>.png
  03.finemap/      finemap.<peak_id>.png
  04.conditional/  conditional.<peak_id>.png
results/_comparison/figures/
                   cohort_compare.png   cohort_manhattan.png
```

`02`–`04` exist only for **genome-wide** peaks of the **additive** scan. The suggestive tier is
represented entirely by panels (b) and (d) of each scan figure, and by the tables.

## The plot block is a fixed physical size

Every figure is authored at **journal double-column width, 183 mm (7.20 in)**, with the `paper` type
scale (base 8 pt, labels 8 pt, ticks 7 pt, panel letters 9 pt).

`caption_block(plot_h=…)` fixes the **plot block** — panel titles, the axes, and the x tick labels and
x-axis label beneath them — at a stated height in inches. The canvas then **grows downward** to hold
whatever caption the figure needs. Crop `plot_h` inches off the top of the PNG and you have the
publication figure, identical every time regardless of how long the caption is.

| figure | layout | plot block (in) | plot / canvas |
|---|---|---|---|
| `scan.<model>` | **2×2**. Left column: (a) Manhattan over (c) its peaks, on a **shared genomic x-axis**. Right column: (b) square QQ, (d) effect/MAF | 7.20 × 6.0 | 72 % |
| `regional.<peak>` | 4 stacked | 7.20 × 6.6 | 74 % |
| `finemap.<peak>` | Left column: (a) association over (b) PIP, **shared genomic x**. Right: (c) cumulative mass, spanning both rows | 7.20 × 3.6 | 67 % |
| `conditional.<peak>` | 1×2 (1.6 : 1) | 7.20 × 2.9 | 65 % |
| `cohort_compare` | (a) two sub-panels on one cohort axis; (b) forest, both full width | 7.20 × 4.4 | 66 % |

**Two genomic panels stacked on a shared axis is this component's recurring shape** — the scan figure
and the fine-mapping figure both use it, because in both cases the lower panel answers a question about
a position in the upper one. Side by side, the reader has to re-locate every position by eye.

**The plot is the deliverable.** These figures are laid out as slides, so the plot block is 64–74 % of
the canvas. It used to be 30–53 %: the caption carried the interpretation, the symbol glossary and the
estimator, and `scan.additive.png` came out 16.2 in tall of which 10.2 in was text. Those three blocks
now live in the sidecar `.md`, which already held the value table, the reading order and the limits.
Nothing was deleted — see *The three layers of explanation* below.

(a) and (b) of the scan figure share an x-axis so they align to the pixel, and **both print chromosome
numbers**; only (b) carries the axis title. The QQ panel is forced square (`set_box_aspect(1)`) with one
limit on both axes — observed and expected are the same quantity in the same units, so any other box
draws the identity line at an angle that is not 45°.

Two consequences shape the panels:

- **`savefig.bbox` is `None`, not `'tight'`.** A tight bbox trims the canvas at save time, so the
  saved image would not be the figure size and a fixed plot block could not be guaranteed. Margins
  are deterministic instead.
- **Narrow panels carry no centred title.** At ~2.2 in a centred title collides with the flush-left
  panel letter, and where a legend sits above the axes it collides with that instead. Panels are
  labelled a/b/c and described in the caption, as journals do.

## Colour is assigned by meaning

`plot_style.py` defines **four semantic layers**. A mark takes its colour from the layer it belongs to,
never from what looked good in one panel, so one meaning is never two colours across 19 figures.

| layer | constant | hex | what it marks |
|---|---|---|---|
| data | `DATA` / `DATA_DARK` | `#3C6E9F` / `#1B3B5F` | the measurement itself — points, bars, the estimate |
| accent | `ACCENT` | `#B0413E` | what must not be missed — the genome-wide line, a genome-wide peak, a gene label |
| reference | `REFERENCE` | `#8A9199` | what the measurement is judged against — null lines, identity lines, the suggestive threshold |
| neutral | `NEUTRAL` | `#D4D8DC` | structure carrying no information — grid, concentration bands, fills |

Two ordered/discrete scales derive from the same family:

- **`COHORT_RAMP`** `#1B3B5F → #3C6E9F → #8FB4D4` — narrow → intermediate → full. Light-to-dark follows
  the nesting, so "more samples" always reads as "lighter".
- **`SERIES`** — five hue-distinct Wong-safe colours for unordered discrete levels (credible sets,
  conditioning rounds).

All are Wong (2011) colourblind-safe hues taken down ~15 % in saturation from the previous set: at
7.2 in the fully-saturated blue and vermilion vibrated against the grid.

**LD *r*² keeps its own five-step scale**, because it is a binned continuous variable and cannot be
expressed in four semantic layers:

| bin | < 0.2 | 0.2–0.4 | 0.4–0.6 | 0.6–0.8 | 0.8–1.0 | unknown |
|---|---|---|---|---|---|---|
| colour | blue | light blue | green | orange | red | **grey** |

Grey means *not in that panel at all*. A variant known to be below a source's reporting floor is
drawn in the lowest bin, not grey — that distinction is the whole reason a truncated panel is
informative. See METHODS §10c.

**Peak tiers** (`TIER_STYLE`, defined once): genome-wide = `ACCENT` diamonds; suggestive = `DATA`
circles.

**Cross-cohort states** (`CALLED_STYLE`, defined once) encode what ONE cohort made of a variant, with
shape and fill so colour stays free for the cohort. **Three** states, not two:

| `called_peak` | marker |
|---|---|
| `genome_wide` | filled diamond |
| `suggestive` | filled circle, smaller |
| `not_a_peak` | open circle |
| `not_in_call_set` | no marker; the row is annotated |

Collapsing `suggestive` and `not_a_peak` into one "open" marker threw away the distinction the reader
most wants — whether that cohort saw anything there at all. The estimate and its interval are drawn in
every case.

**Peak lead** in regional and fine-mapping figures: purple diamond (`LEAD_COLOR`), a fifth colour
reserved for exactly one thing.

## Collisions are prevented by measurement, not by tuned offsets

A fixed offset works until the data changes shape. Four helpers in `plot_style.py` measure the rendered
figure instead, and every script routes through them:

| helper | the defect class it closes |
|---|---|
| `gene_labels` | locus annotation in the **gwaslab** idiom (see below) — a label band above the data, rotated italic text, and an L-shaped leader arm. |
| `repel_from_centre` | separates 1-D positions by pushing overlapping pairs apart **from their midpoint**, then enforces the gap directly against the axis ends. |
| `fit_left_margin` | the left margin is measured from the widest y tick label and y-axis label. A hardcoded `left=0.088` clipped `intermediate` off the canvas edge; `caption_block(left='auto')` cannot. Twin axes are excluded, so a right-hand `N_eff` axis cannot inflate the left margin. |
| `value_labels` | a point's value label is offset off its marker **and grows the axis to make room for itself**. λ = 1.132 used to print past the top of its panel; `dx` additionally nudges a label clear of a vertical element such as a stem. |
| `spread_labels` | text annotations whose rendered boxes overlap are separated along one axis, greedily, on measured boxes. Replaces the `k % 2` alternating-row hacks, which staggered labels whether or not they collided and still overlapped when three landed together. |
| `place_legend` | above the axes if a single row fits the panel width, otherwise inside the corner holding the fewest **actual data points**. The measurement decides the *placement*: wrapping an over-wide strip to a second row above the axes only moves the overflow into the vertical, where it is clipped. |


It measures the axes box, so it **must run after `caption_block`**, which resizes the figure and
re-runs `subplots_adjust`. Asking for bands at 0.70 and 0.92 before that got both at ~0.55, overlapping.

Plus `thin_tick_labels`, which blanks tick labels that would overlap their kept neighbour — chr19–22
collide at any readable size on a whole-genome axis.

### Locus annotation follows gwaslab, `anno_style = expand`

`gene_labels()` takes gwaslab's own style names. `params.AnnoStyle` selects one and reaches the
figures through the pipeline, so the choice is recorded in `run_manifest.json`:

| style | positions | rotation |
|---|---|---|
| `right` | greedy left-to-right; gwaslab's own default | 40° |
| `expand` | symmetric repel | 90° |
| `tight` | as `expand`, with a short arm | 90° |
| **`auto`** | **symmetric repel — this component's default** | **measured: 0°, else 40°, else 90°** |

`auto` takes the least intrusive rotation the labels actually fit in, measured against the axis width
with the same `w·cos θ + h·sin θ` footprint the repel uses, so the choice can never drift from the
spacing. Always-vertical is right for eleven names on a peak panel and needlessly heavy for the two on
a Manhattan. **`expand` is unchanged and one parameter away** — `params.AnnoStyle = 'expand'` restores
the always-vertical look with no other change.

### Labels live OUTSIDE the data area

`gene_labels()` reserves its strip by **shrinking the axes box**, then draws above it with
`clip_on=False`. The data limits are never inflated to make room for text.

They used to be, and it was the worst thing about these figures: to hold a band of 90° labels the peak
panel's ceiling was pushed to `max(-log10 P) × 1.95`, so the peaks occupied roughly a third of the
panel and the rest was reserved whitespace. With the labels moved into a strip outside the axes the
data fills most of it — and the remainder is not slack, it is the two threshold lines: the axis has to
contain the suggestive line below the lowest peak and the genome-wide line above the highest, and
where no peak clears the genome-wide line that line alone sets the ceiling.

Consequences the callers handle, via two helpers:

- `align_panel_tops(*axes)` — only the labelled panel shrinks, so its unlabelled neighbour in the same
  gridspec row must be brought down to match or the row stops reading as a row.
- `strip_pad(ax)` — the panel letter lives in the title slot, which is attached to the axes, so after
  the strip is reserved it would print *below* the gene names it introduces. One pad per row, the
  larger of the two, keeps the letters on one line.

**One band per panel, not one per tier.** The peak panel's genome-wide and suggestive names go into a
single `gene_labels()` call, sorted by position, with the tier carried by **colour and weight alone** —
bold `ACCENT` against regular `REFERENCE`. Two bands cost two strips (0.98 in of text above a 1.00 in
axes) and put two text orientations in one panel; merged, the strip is 0.76 in, there is one
orientation, and the repel sees every label at once so red and grey names are spaced against each other
rather than in two independent passes. `gene_labels()` still stacks a second band above a first if a
caller asks for one, but nothing here does.

**Row proportions follow from the strip.** The labelled panel gives up its strip, so its gridspec row
needs the larger `height_ratios` share (1.5 against 1.0) or the axes ends up half the height of the
unlabelled row above it. With that, `align_panel_tops` brings the row's summary panel down to match and
both land near 1.9 in; without the extra share it produced two 1.0 in panels under a 0.98 in strip,
which is the arrangement that looked wrong.

The minimum separation is `max(measured footprint, repel_force × span)`, `repel_force` = 0.03 as in
gwaslab. gwaslab spaces on `repel_force × span` alone, which overlaps as soon as the labels are wider
than that fraction of the axis; the measured footprint is the floor that makes a collision impossible,
and `repel_force` on top lets a caller spread labels *further* than they strictly need.

`gene_labels()` **measures its own result** and prints a warning to stderr if any two labels touch —
the spacing depends on rendered text metrics, so it cannot be verified by reading the code. All nine
scan figures currently report zero.

### Geometry

Matched against the user's own install
(`gwaslab/viz/viz_aux_annotate_plot.py`, v4.x):

- labels sit in **one band above the data**, not at a fixed offset from each point;
- text is **italic and rotated** — 40° normally, **90°** when dense. The switch is a *measurement*:
  if the labels laid side by side would not fit the axis width, they go vertical;
- the leader is an **L-shaped arm**, `arc,angleA=0,armA=0,angleB=90,armB=<px>` — vertical up from the
  peak, then horizontal to the text. The vertical segment is what marks the peak's true position, so
  the horizontal segment can absorb any amount of displacement without the reader losing the anchor;
- positions are repelled **symmetrically from the midpoint of each overlapping pair**, then forced to
  honour the minimum gap within `[lo, hi]`.

Two things this got wrong before they were measured, both worth keeping written down:

1. A greedy left-to-right sweep only ever pushes **right**, so a dense cluster drags every later label
   with it — the chr2–6 run ended up shifted past chr6, arms pointing at peaks it did not belong to.
2. The horizontal footprint of rotated text is `w·cos θ + h·sin θ`. Dropping the `h·sin θ` term
   understates a 40° label by its own line height, which was enough for two names in the dominant
   band to print on top of each other.

`Annotation.get_window_extent()` unions in the **arrow patch**, so a collision check that uses it
measures the arms, not the text. Use `Text.get_window_extent(annotation, renderer=…)`.

Where a label must sit over a populated region anyway (the threshold labels at the right edge of a
Manhattan), it carries an opaque white backing rather than being moved somewhere less meaningful.

## The three layers of explanation

| layer | where | what it carries |
|---|---|---|
| in-figure caption | below the panels of every PNG | **three things only**: the finding in one bold sentence, one line per panel, one `Data:` line. Nothing else — the figure is being laid out as a slide |
| sidecar `.md` | `<figure>.md` | the same panels expanded, **plus the concrete numbers for that rendering**, full statistics tables a figure cannot hold, a numbered reading order, and an explicit "what this figure does NOT establish" section |
| `METHODS.md` | `docs/` | the reasoning the figures deliberately do not editorialise about |

Every sidecar carries the same eight sections in the same order: the question · panels ·
**interpretation** · values in this rendering · full statistics · how to read it · what it does NOT
establish · symbols and **model**. The two bold ones moved out of the PNG.

A symbol that appears in a figure is defined in its sidecar. The
definitions come from one place — `SYMBOL_DEFS` in `scripts/plot_style.py`, as `(label, meaning)` — so
they cannot drift between a figure, its sidecar and METHODS §2. Each is rendered on its **own line**:
run together into a paragraph, six definitions read as an undifferentiated wall.

## Rules the figures follow

1. **No downsampling.** Manhattan and QQ draw all ~4.8–5.1 M analysed variants. The scatter is
   rasterised so the file stays a few MB while axes, ticks and text remain vector. Thinning would
   distort the QQ's lower arm, which is where λ is read.
2. **Nothing overlaps.** Panel letters sit in the title slot (`panel_tag`), where matplotlib reserves
   space, so they cannot land on data. Captions are laid out from renderer-measured text, and the
   figure grows rather than letting axes overrun a tall caption.
3. **Every symbol is defined where it is used**, so a reader never has to infer whether β is an effect
   size or an artefact measure.
4. **Each scan is annotated with its own result.** A Manhattan labels the genome-wide loci of *its*
   model, by gene symbol. Labelling a dominant or recessive scan with the additive peak list says
   nothing about the scan being drawn.
5. **A deflated scan says so in its own caption.** λ < 0.9 fires an explicit sentence in the caption
   rather than leaving the caveat to METHODS.
6. **Deterministic output.** Where points would otherwise be jittered, the offsets are computed from
   rank (`np.linspace`), never sampled, so re-running produces the same figure.

## What each figure answers

| figure | the one question |
|---|---|
| `scan.<model>.png` | where does **this** scan show association, what are its peaks, and is it calibrated enough to believe? |
| | *Peaks and labels are that model's own.* **Both genomic panels are labelled**: (a) carries the genome-wide loci, (c) those plus the `params.LabelSuggestive` smallest-*P* suggestive ones. a full suggestive tier does not fit a 7.2 in axis, so the cap is stated in the caption and the figure never implies the unlabelled peaks do not exist. |
| `regional.<peak>.png` | is the LD structure producing this peak a population property or an artefact of this sample? |
| `finemap.<peak>.png` | does the posterior concentrate on a few variants, or does LD spread it out? |
| `conditional.<peak>.png` | does this peak carry one signal, or more than one? |
| `cohort_compare.png` | how do calibration and effect size move as the ancestry filter is relaxed — for **every** genome-wide lead in **all three** cohorts? |
| | *Additive only.* (a) is two sub-panels on one cohort axis — λ_GC as a deviation from 1 over the case/control composition — rather than a dual-axis plot; the stacked bar is what shows that `narrow` keeps 95 % of the cases against only 67 % of the controls. |
| `cohort_manhattan.png` | what does the genome look like in each of the three nested sample sets, side by side? |
| | *Additive only.* Each row is one cohort on the **same 2.5 : 1 grid as the scan figure's top row** — Manhattan left, QQ right — so each Manhattan is **4.08 × 1.77 in, the identical box `scan.<model>.png` gives its panel (a)**, and a peak has the same shape on both figures. Matching the aspect by construction, not by arithmetic, is the point: a full-width row is 6.50 in across and squashes the same data to aspect 4.6 against 2.3. The header carries composition (cases, controls, N_eff), the QQ carries calibration (λ_GC) — the scan figure's division of labour, and it keeps λ_GC off the row twice. |
| | **Five things are shared and each is a correctness requirement**: one cumulative offset map (per-cohort offsets would put chr16 at three different x), one −log₁₀ *P* limit (so a taller peak is a stronger peak, not a rescaled axis), one data height per row (a shared limit on unequal heights is not a shared scale — `plot_style.equalise_row_heights`), one fixed y-tick locator (rows differing by hundredths of an inch otherwise flip to a different tick step), and one QQ limit. `align_panel_tops` is **not** used here: `MQ.qq` sets `box_aspect(1)`, so it would compare the Manhattan to a squared box and shrink the *Manhattan* onto it. Chromosome numbers are on the bottom row only. The genome-wide loci are named **once**, above the top row, over the **union** across cohorts, in `INK` rather than `ACCENT` — red is left to the per-row diamonds, which mean "this cohort called it". Dotted guides drop through all three panels so a locus can be read *down* the figure. |

## Deliberate omissions

- **No power / minimum-detectable-effect figure.** Removed by decision (METHODS §13). With no
  replication cohort available, its only load-bearing statement was that every peak sits at the
  detection boundary — arithmetically the same statement as "the peak barely cleared the threshold",
  which *P* already says.
- **No winner's-curse-corrected effect.** Removed by scope decision. The cross-cohort forest carries
  the equivalent caution by showing each estimate in all three nested sample sets.
- **No observed-against-null panel.** The excess of suggestive variants over *M*α restates the
  inflation λ_GC already reports, and *M*α assumes the *M* tests are independent — LD makes that false,
  so the ratio has no defensible scale.
- **No per-peak figures for the suggestive tier.** 128 regional plots of noise would assert findings
  the tier cannot support.
- **No LDSC panel.** Not usable at N_eff ≈ 1,400 (METHODS §7).
