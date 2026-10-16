# Engineering Traverse Solver

A self-contained Engineering 1035 desktop calculator for closed traverses with two omitted quantities. All application code and regression tests are in **`traverse_solver.py`**.

## Run

Requires **Python 3.10+**, Tkinter, and a graphical desktop session:

```bash
python traverse_solver.py
```

Tkinter ships with most desktop Python installations. On Debian/Ubuntu, install `python3-tk` if needed. The solver works offline and requires no NumPy or SciPy. Matplotlib is optional (`python -m pip install matplotlib`); without it, the application draws an equal-scale plot using Tkinter Canvas.

1. Choose 3–20 sides and click **Generate Traverse**. Generating starts a fresh input table.
2. Enter positive known lengths and bearings. Mark exactly two quantities **Unknown** using the status selectors; disabled input values are ignored.
3. Use azimuths clockwise from North, or select **Quadrantal Bearing** and use the N/S and E/W controls. Switching formats converts valid known bearings automatically.
4. Enter angles as decimal degrees (`125.5`), DMS (`125°30'00"`), or colon-separated DMS (`125:30:00`). In azimuth mode, full quadrantal text such as `N 30 E` also works.
5. Optionally enter initial coordinates and select metres or feet. The unit selector labels the supplied numbers; it does not convert them.
6. Click **Solve Traverse**. Review solved quantities, closure, area, the complete side table, coordinates, and calculation details in the tabs.
7. If multiple solutions exist, review each with the solution selector. Extra engineering constraints are required to select between them. Every solution is included in TXT/CSV exports.
8. Use **Plot Traverse** or **Export Results**. Matplotlib's plot toolbar also supports image saving. Editing inputs invalidates previous results and closes their plot.

**Load Example** provides an 80 × 100 rectangle with omitted values for each of the four cases. Select the case beside the button, load it, then solve. Case IV demonstrates two valid solutions; one crosses itself, so its polygon area is unavailable.

## Mathematics and limits

The program uses `ΔE = L sin(θ)` and `ΔN = L cos(θ)` with azimuth clockwise from North. It preserves full floating-point precision internally.

- **Case I:** direct closure vector, length by `hypot`, bearing by `atan2(E, N)`.
- **Case II:** ray/circle intersection, retaining all distinct positive-length roots.
- **Case III:** a 2 × 2 linear system, rejecting parallel or nearly parallel directions and nonpositive results.
- **Case IV:** exact circle/circle intersection, enumerating both mirror solutions and deduplicating tangent solutions. This analytical strategy replaces iterative root searches and needs no initial guesses. Concentric equal circles are explicitly reported as infinitely ambiguous.

Each candidate is substituted into both closure equations and must meet a length-scaled tolerance. Area uses local coordinates to reduce cancellation from large initial coordinate offsets. Crossed, touching, overlapping, and degenerate polygons produce an area warning rather than a misleading area. No traverse adjustment is applied.

Closure is imposed by solving the omitted values. “Perfect / Numerical Closure” describes the mathematical residual and does **not** establish measurement accuracy. Very nearly parallel bearings and quantities smaller than the numerical tolerance may be indeterminate at floating-point precision.

## Verification

Run the embedded numerical tests without a display or optional packages:

```bash
python traverse_solver.py --self-test
```

The suite covers all four examples, 240 rotated/scaled pentagon problems, decimal/DMS and quadrant conversion, multiple roots, tangency, impossible/singular/indeterminate data, area validity, coordinate offsets, and both report formats.

With a graphical display available, exercise the actual GUI widgets and plotting/export paths:

```bash
python traverse_solver.py --gui-test
```

The GUI smoke test covers all four cases in both input formats, alternative solution selection, Matplotlib and built-in plotting, TXT/CSV export dialogs, 20-side generation, clear, and invalid input. On a headless Linux machine, a virtual X display is needed for this check; numerical tests remain display-independent.
