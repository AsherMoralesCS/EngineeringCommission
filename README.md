# Engineering Traverse Solver

Engineering 1035 desktop calculator for closed traverses with two omitted quantities. The student application is **one runnable file**, `traverse_solver.py`; regression tests are separate. Requires Python 3.10+ and Tkinter.

```bash
python traverse_solver.py
```

Tkinter normally ships with desktop Python (Debian/Ubuntu: `python3-tk`). The solver and input/report files work offline using the standard library. **Plotting requires Matplotlib**:

```bash
python -m pip install matplotlib
```

If Matplotlib is absent, incompatible, or fails during initialization, the app still opens with numerical solving and exports available. It shows an installation message in the plot pane. There is one plotting backend and no duplicate Canvas renderer.

## Using the calculator

- Set **Sides** to 3–20. Spinner arrows update the table immediately; typed counts apply on Enter or when leaving the field. Existing rows are preserved. Scroll the table with the wheel or scrollbar; drag the horizontal divider to give it more room.
- Enter positive lengths and mark omitted values with the **Unknown** checkboxes. Case detection uses those flags directly. Known bearings accept azimuth (`125.5`), DMS (`125°30'00"` or `125:30:00`), or quadrantal text (`N 30° E`). All internal bearings are clockwise from North: 0° North, 90° East.
- Set initial Easting/Northing and units. Switching **m/ft** relabels values and retains results; it does not convert distances. Changes to measurements or starting coordinates invalidate the previous solution.
- **Load Example** supplies valid rectangle data for each case. **Solve Traverse** shows omitted values, perimeter, area, calculation details, and a scrollable table with side components and station coordinates.
- Multiple solutions appear in the selector. Simple polygons are listed first; crossed or degenerate solutions are explicitly flagged and have no area. They remain available because they satisfy the equations, but require additional engineering constraints to select.
- The solved plot uses equal scales, station labels, arrows, and local coordinates relative to the supplied origin. Known sides are blue and solved sides red. **Enlarge Plot** opens a larger view; Matplotlib's toolbar provides zoom, pan, and image saving. No assumed unsolved geometry is plotted.
- **Save Inputs / Load Inputs** round-trip JSON, preserving text, unknown flags, side count, coordinates, and units, including incomplete drafts. **Export Results** writes all candidates as a formatted TXT report or CSV. Exported reports also include calculations and coordinates.

The window leaves room for desktop panels and title bars based on the reported screen dimensions. Both the input/result split and the result/plot split are resizable. The status bar is reserved at the bottom; text and calculation tables scroll.

## Numerical interpretation

`ΔE = L sin(θ)` and `ΔN = L cos(θ)`. Four independent case functions return candidate lists; a separate explanation function builds the derivation. `solve_traverse` returns a `SolveResult` dataclass.

- **Case I:** direct closure vector and `atan2(E, N)`.
- **Case II:** positive ray/circle intersections, retaining all isolated positive roots.
- **Case III:** a 2 × 2 length system. Parallel directions are rejected. If the unit-vector determinant has magnitude below **0.01**, a warning and the **2-norm condition number** appear in the results. Check bearings carefully: small input errors can cause large length errors.
- **Case IV:** circle/circle intersections, including mirror solutions and tangent deduplication. Equal, opposite free bearings are reported as infinitely ambiguous.

Lengths, squared quantities, and dimensionless matrix diagnostics have separately named tolerances in `NumericalPolicy`. Length residuals and deduplication scale with traverse size; a machine-epsilon bound handles square-root roundoff at tangency. Area/intersection checks operate on local coordinates, before adding large global offsets.

A negative solved length suggests checking for a reversed bearing (180°), the bearing quadrant, side order, or a mistyped length. Errors give these checks explicitly.

**Closure is imposed by solving the equations.** Numerical residuals appear in Calculation Details as substitution checks, never as a survey accuracy headline or a “Perfect” precision rating. This app does not adjust measured traverses: Bowditch/compass-rule adjustment and DMD area are outside this omitted-measurement workflow. Area uses the coordinate/shoelace method and is unavailable for crossed, touching, retraced, or degenerate boundaries.

## Tests

```bash
python -m unittest -v test_traverse_solver
```

Covers all four examples, nonrectangular cases, tangent/multiple solutions, singular/impossible data, conditioning warnings, actionable negative-length errors, large coordinate offsets, strict ASCII numeric parsing, DMS/cardinal formatting, input persistence, and both report formats.

With a desktop display:

```bash
python -m unittest -v test_gui
```

GUI regression tests cover 1366×768 screen fit, footer visibility, 20-row mouse scrolling, preserving results when relabeling units, input save/load actions, local plots at large coordinate offsets, crossed-solution labels, and recovery from Matplotlib initialization errors. GUI tests are skipped when no display is available.
