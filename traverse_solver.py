#!/usr/bin/env python3
"""Engineering 1035 closed-traverse calculator. Run with Python 3.10+.

No computational packages are required. Matplotlib is optional. Use --self-test
for the embedded numerical tests, or --gui-test with a display for a GUI smoke test.
Azimuths are measured clockwise from North; vectors are (Easting, Northing).
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass, replace
import io
import math
from pathlib import Path
import re
import sys
import unittest

CASES = (
    "Case I: Length & Bearing of One Side Omitted",
    "Case II: Length of One Side & Bearing of Another Omitted",
    "Case III: Length of Two Sides Omitted",
    "Case IV: Bearing of Two Sides Omitted",
)


class TraverseError(ValueError):
    """A user-correctable input or geometric problem."""


def number(text, label="Value"):
    try:
        value = float(text)
    except (ValueError, TypeError):
        raise TraverseError(f"{label} must be a number.") from None
    if not math.isfinite(value):
        raise TraverseError(f"{label} must be finite.")
    return value


def parse_angle(text):
    """Parse decimal degrees, symbol DMS, or colon/space-separated DMS."""
    text = str(text).strip().replace("º", "°").replace("′", "'").replace("″", '"')
    decimal = re.fullmatch(r"([+]?\d+(?:\.\d*)?|[+]?\.\d+)\s*°?", text)
    if decimal:
        return number(decimal[1], "Angle")
    patterns = (
        r"(\d+)\s*°\s*(\d+(?:\.\d*)?)\s*'\s*(?:(\d+(?:\.\d*)?)\s*\")?",
        r"(\d+)\s*[:\s]\s*(\d+(?:\.\d*)?)\s*[:\s]\s*(\d+(?:\.\d*)?)",
    )
    for pattern in patterns:
        match = re.fullmatch(pattern, text)
        if match:
            degrees, minutes, seconds = (float(x or 0) for x in match.groups())
            if not (0 <= minutes < 60 and 0 <= seconds < 60):
                raise TraverseError("DMS minutes and seconds must each be in [0, 60).")
            return degrees + minutes / 60 + seconds / 3600
    raise TraverseError("Use decimal degrees (125.5) or DMS (125°30'00\").")


def parse_bearing(text):
    text = str(text).strip().upper()
    match = re.fullmatch(r"([NS])\s*(.*?)\s*([EW])", text)
    if match:
        north_south, angle, east_west = match.groups()
        angle = parse_angle(angle)
        if not 0 <= angle <= 90:
            raise TraverseError("Quadrantal angles must be between 0° and 90°.")
        return {("N", "E"): angle, ("S", "E"): 180 - angle,
                ("S", "W"): 180 + angle, ("N", "W"): 360 - angle}[
                    north_south, east_west] % 360
    angle = parse_angle(text)
    if not 0 <= angle <= 360:
        raise TraverseError("Azimuth must be between 0° and 360° (360° equals North).")
    return angle % 360


def dms(angle):
    # Round total seconds first so 59.999 seconds carries into the next minute.
    hundredths = round(angle * 3600 * 100)
    degrees, remainder = divmod(hundredths, 360000)
    minutes, seconds = divmod(remainder, 6000)
    return f"{degrees}°{minutes:02d}'{seconds / 100:05.2f}\""


def quadrant_parts(azimuth):
    angle = azimuth % 360
    if angle <= 90:
        return "N", angle, "E"
    if angle <= 180:
        return "S", 180 - angle, "E"
    if angle <= 270:
        return "S", angle - 180, "W"
    return "N", 360 - angle, "W"


def quadrantal(azimuth):
    ns, angle, ew = quadrant_parts(azimuth)
    return f"{ns} {dms(angle)} {ew}"


@dataclass(frozen=True)
class TraverseSide:
    length: float | None
    azimuth: float | None


def detect_case(sides):
    lengths = [i for i, side in enumerate(sides) if side.length is None]
    bearings = [i for i, side in enumerate(sides) if side.azimuth is None]
    if len(lengths) == len(bearings) == 1:
        return 0 if lengths == bearings else 1
    if len(lengths) == 2 and not bearings:
        return 2
    if len(bearings) == 2 and not lengths:
        return 3
    raise TraverseError("Mark exactly two unknown quantities: length and bearing on one side; "
                        "length and bearing on different sides; two lengths; or two bearings.")


def unit_vector(angle):
    angle = math.radians(angle)
    return math.sin(angle), math.cos(angle)


def vector(side):
    east, north = unit_vector(side.azimuth)
    return side.length * east, side.length * north


def azimuth(east, north):
    result = math.degrees(math.atan2(east, north)) % 360
    # A tiny negative angle can round to exactly 360 during modulo reduction.
    return 0.0 if result == 360 else result


def display_azimuth(angle):
    return f"{round(angle % 360, 4) % 360:.4f}"


@dataclass
class Solution:
    sides: list[TraverseSide]
    coordinates: list[tuple[float, float]]  # (E, N), including final endpoint
    components: list[tuple[float, float]]
    perimeter: float
    sum_e: float
    sum_n: float
    closure: float
    area: float | None
    warnings: list[str]
    tolerance: float

    @property
    def precision(self):
        if self.closure <= self.tolerance:
            return "Perfect / Numerical Closure"
        return f"1 : {self.perimeter / self.closure:,.0f}"


def polygon_area(points, tolerance):
    """Return unsigned shoelace area, rejecting crossed/touching polygons."""
    # Work relative to the origin to avoid cancellation with large site coordinates.
    origin_e, origin_n = points[0]
    points = [(e - origin_e, n - origin_n) for e, n in points[:-1]]
    count = len(points)
    extent = max((math.hypot(e, n) for e, n in points), default=1)
    cross_tolerance = tolerance * max(extent, 1)

    def cross(a, b, c):
        return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])

    def on_segment(a, b, p):
        return (abs(cross(a, b, p)) <= cross_tolerance
                and min(a[0], b[0]) - tolerance <= p[0] <= max(a[0], b[0]) + tolerance
                and min(a[1], b[1]) - tolerance <= p[1] <= max(a[1], b[1]) + tolerance)

    for i in range(count):
        a, b = points[i], points[(i + 1) % count]
        previous = points[(i - 1) % count]
        # Consecutive segments may share an endpoint, but must not overlap.
        if (abs(cross(previous, a, b)) <= cross_tolerance
                and (previous[0] - a[0]) * (b[0] - a[0])
                + (previous[1] - a[1]) * (b[1] - a[1]) > tolerance ** 2):
            return None, ["Area unavailable: adjacent sides overlap or retrace."]
        for j in range(i + 1, count):
            if j == i + 1 or (i == 0 and j == count - 1):
                continue
            c, d = points[j], points[(j + 1) % count]
            ab_c, ab_d, cd_a, cd_b = cross(a, b, c), cross(a, b, d), cross(c, d, a), cross(c, d, b)
            proper_cross = ((ab_c > cross_tolerance and ab_d < -cross_tolerance
                             or ab_c < -cross_tolerance and ab_d > cross_tolerance)
                            and (cd_a > cross_tolerance and cd_b < -cross_tolerance
                                 or cd_a < -cross_tolerance and cd_b > cross_tolerance))
            if proper_cross or any((on_segment(a, b, c), on_segment(a, b, d),
                                    on_segment(c, d, a), on_segment(c, d, b))):
                return None, ["Area unavailable: nonadjacent sides intersect or touch. "
                              "Closure alone does not establish a simple polygon."]
    area = abs(math.fsum(points[i][0] * points[(i + 1) % count][1]
                         - points[(i + 1) % count][0] * points[i][1]
                         for i in range(count))) / 2
    if area <= cross_tolerance:
        return None, ["Area unavailable: the traverse is degenerate or collinear."]
    return area, []


class TraverseSolver:
    """Analytical solvers enumerate every isolated, positive-length solution."""

    def solve(self, sides, initial_e=0.0, initial_n=0.0):
        if not 3 <= len(sides) <= 20:
            raise TraverseError("Number of sides must be from 3 to 20.")
        initial_e = number(initial_e, "Initial Easting")
        initial_n = number(initial_n, "Initial Northing")
        for i, side in enumerate(sides, 1):
            if side.length is not None and (not math.isfinite(side.length) or side.length <= 0):
                raise TraverseError(f"Side {i}: length must be a positive finite number.")
            if side.azimuth is not None and (not math.isfinite(side.azimuth) or not 0 <= side.azimuth < 360):
                raise TraverseError(f"Side {i}: azimuth must be normalized to [0°, 360°).")
        case = detect_case(sides)
        known = [vector(s) for s in sides if s.length is not None and s.azimuth is not None]
        known_e = math.fsum(v[0] for v in known)
        known_n = math.fsum(v[1] for v in known)
        target_e, target_n = -known_e, -known_n
        scale = max(math.fsum(s.length or 0 for s in sides), 1e-100)
        tolerance = max(scale * 1e-10, 1e-110)
        details = [CASES[case], "", "Azimuth is clockwise from North.",
                   "ΔE = L sin(θ); ΔN = L cos(θ)",
                   f"Fully known sides: ΣΔE = {known_e:.12g}; ΣΔN = {known_n:.12g}",
                   f"Unknown-side vectors must sum to E = {target_e:.12g}, N = {target_n:.12g}.",
                   "Closure equations: ΣΔE = 0; ΣΔN = 0.", ""]
        length_ids = [i for i, s in enumerate(sides) if s.length is None]
        bearing_ids = [i for i, s in enumerate(sides) if s.azimuth is None]
        candidates = []
        if case == 0:
            i = length_ids[0]
            length = math.hypot(target_e, target_n)
            if length <= tolerance:
                raise TraverseError("The omitted side has zero length; its bearing is indeterminate.")
            candidate = list(sides)
            candidate[i] = TraverseSide(length, azimuth(target_e, target_n))
            candidates.append(candidate)
            details += [f"Side {i + 1}: L = √(({target_e:.12g})² + ({target_n:.12g})²)",
                        f"θ = atan2({target_e:.12g}, {target_n:.12g}) mod 360°."]
        elif case == 1:
            i, j = length_ids[0], bearing_ids[0]
            ue, un = unit_vector(sides[i].azimuth)
            radius = sides[j].length
            # Intersect a ray x*u with the circle |target - x*u| = radius.
            projection = target_e * ue + target_n * un
            perpendicular = target_e * un - target_n * ue
            normalized_disc = (radius / scale) ** 2 - (perpendicular / scale) ** 2
            details += [f"Side {i + 1} unknown length x; side {j + 1} unknown bearing φ.",
                        f"{ue:.12g} x + {radius:.12g} sin(φ) = {target_e:.12g}",
                        f"{un:.12g} x + {radius:.12g} cos(φ) = {target_n:.12g}",
                        f"Eliminate φ: x² - 2({projection:.12g})x + "
                        f"({math.hypot(target_e, target_n):.12g})² - ({radius:.12g})² = 0.",
                        "Intersect the positive-length ray with the known-length circle; retain all positive roots."]
            if normalized_disc < -1e-12:
                raise TraverseError("No solution: the bearing-unknown side cannot reach the known-bearing ray.")
            if abs(normalized_disc) <= 32 * math.ulp(max((radius / scale) ** 2, (perpendicular / scale) ** 2)):
                normalized_disc = 0.0
            root = scale * math.sqrt(max(0, normalized_disc))
            for length in (projection + root, projection - root):
                if length <= tolerance:
                    continue
                candidate = list(sides)
                candidate[i] = replace(sides[i], length=length)
                candidate[j] = replace(sides[j], azimuth=azimuth(target_e - length * ue, target_n - length * un))
                candidates.append(candidate)
        elif case == 2:
            i, j = length_ids
            ae, an = unit_vector(sides[i].azimuth)
            be, bn = unit_vector(sides[j].azimuth)
            determinant = ae * bn - be * an
            details += [f"{ae:.12g} L{i + 1} + {be:.12g} L{j + 1} = {target_e:.12g}",
                        f"{an:.12g} L{i + 1} + {bn:.12g} L{j + 1} = {target_n:.12g}",
                        f"2 × 2 determinant = {determinant:.12g}; solve using Cramer's rule."]
            if abs(determinant) < 1e-10:
                raise TraverseError("Singular or ill-conditioned system: the two known bearings are parallel "
                                    "or nearly parallel. Missing lengths cannot be determined reliably.")
            first = (target_e * bn - be * target_n) / determinant
            second = (ae * target_n - target_e * an) / determinant
            if min(first, second) <= tolerance:
                raise TraverseError("No physical solution: closure requires a zero or negative side length.")
            candidate = list(sides)
            candidate[i] = replace(sides[i], length=first)
            candidate[j] = replace(sides[j], length=second)
            candidates.append(candidate)
        else:
            i, j = bearing_ids
            a, b = sides[i].length, sides[j].length
            distance = math.hypot(target_e, target_n)
            details += [f"{a:.12g} sin(θ{i + 1}) + {b:.12g} sin(θ{j + 1}) = {target_e:.12g}",
                        f"{a:.12g} cos(θ{i + 1}) + {b:.12g} cos(θ{j + 1}) = {target_n:.12g}",
                        "Use exact circle-circle intersection to enumerate both mirror solutions.",
                        "This geometric method needs no initial guesses and cannot miss an isolated root.",
                        "Let R = |target|, x = (R² + a² - b²)/(2R), h = √(a² - x²).",
                        "First side vector = x·unit(target) ± h·perpendicular(unit(target)).",
                        "Second side vector = target − first; θ = atan2(E, N) mod 360°."]
            if distance <= tolerance:
                if abs(a - b) <= tolerance:
                    raise TraverseError("Infinitely many solutions: equal opposite sides can rotate freely. "
                                        "An additional bearing or engineering constraint is required.")
                raise TraverseError("No solution: unequal side lengths cannot sum to a zero vector.")
            if distance > a + b + tolerance or distance < abs(a - b) - tolerance:
                raise TraverseError("No solution: the known lengths violate the triangle inequality.")
            r, aa, bb = distance / scale, a / scale, b / scale
            x = (r * r + aa * aa - bb * bb) / (2 * r)
            height_sq = aa * aa - x * x
            if height_sq < -1e-12:
                raise TraverseError("No stable circle intersection exists for these inputs.")
            if abs(height_sq) <= 32 * math.ulp(max(aa * aa, x * x)):
                height_sq = 0.0
            height = math.sqrt(max(0, height_sq)) * scale
            x *= scale
            ue, un = target_e / distance, target_n / distance
            for sign in (1, -1):
                east, north = x * ue - sign * height * un, x * un + sign * height * ue
                candidate = list(sides)
                candidate[i] = replace(sides[i], azimuth=azimuth(east, north))
                candidate[j] = replace(sides[j], azimuth=azimuth(target_e - east, target_n - north))
                candidates.append(candidate)
        solutions = []
        for candidate in candidates:
            if any(self._same(candidate, old.sides) for old in solutions):
                continue
            solution = self._evaluate(candidate, initial_e, initial_n)
            if solution.closure <= solution.tolerance:
                solutions.append(solution)
        if not solutions:
            raise TraverseError("No valid positive-length solution satisfies closure within tolerance.")
        details += ["", f"Verified {len(solutions)} distinct solution(s) by substitution into both closure equations."]
        for index, solution in enumerate(solutions, 1):
            details += [f"", f"Solution {index}:"]
            for i, (original, solved) in enumerate(zip(sides, solution.sides), 1):
                if original.length is None:
                    details.append(f"Side {i} length = {solved.length:.12g}")
                if original.azimuth is None:
                    details.append(f"Side {i} bearing = {solved.azimuth:.12g}° = {quadrantal(solved.azimuth)}")
            details.append(f"ΣΔE = {solution.sum_e:.12g}; ΣΔN = {solution.sum_n:.12g}; "
                           f"closure error = {solution.closure:.12g}")
        return case, solutions, "\n".join(details)

    @staticmethod
    def _same(first, second):
        return all(math.isclose(a.length, b.length, rel_tol=1e-9, abs_tol=1e-110)
                   and abs((a.azimuth - b.azimuth + 180) % 360 - 180) < 1e-7
                   for a, b in zip(first, second))

    @staticmethod
    def _evaluate(sides, initial_e, initial_n):
        components = [vector(s) for s in sides]
        perimeter = math.fsum(s.length for s in sides)
        tolerance = max(perimeter * 1e-9, 1e-110)
        sum_e = math.fsum(e for e, _ in components)
        sum_n = math.fsum(n for _, n in components)
        # Compute local coordinates first; translate only for the displayed site coordinates.
        local = [(0.0, 0.0)]
        for i in range(len(sides)):
            local.append((math.fsum(e for e, _ in components[:i + 1]),
                          math.fsum(n for _, n in components[:i + 1])))
        area, warnings = polygon_area(local, tolerance)
        coordinates = [(initial_e + e, initial_n + n) for e, n in local]
        if any(not math.isfinite(v) for point in coordinates for v in point):
            raise TraverseError("Coordinates exceed the supported floating-point range.")
        return Solution(sides, coordinates, components, perimeter, sum_e, sum_n,
                        math.hypot(sum_e, sum_n), area, warnings, tolerance)


def example(case):
    """80 × 100 rectangle, with the required values genuinely omitted."""
    sides = [TraverseSide(80, 0), TraverseSide(100, 90),
             TraverseSide(80, 180), TraverseSide(100, 270)]
    if case == 0:
        sides[3] = TraverseSide(None, None)
    elif case == 1:
        sides[0] = replace(sides[0], length=None)
        sides[2] = replace(sides[2], azimuth=None)
    elif case == 2:
        sides[0] = replace(sides[0], length=None)
        sides[1] = replace(sides[1], length=None)
    else:
        sides[0] = replace(sides[0], azimuth=None)
        sides[1] = replace(sides[1], azimuth=None)
    return sides


TABLE_HEADERS = ("Side", "Length", "Bearing", "Azimuth (°)", "ΔN / Latitude", "ΔE / Departure",
                 "Start Northing", "Start Easting", "End Northing", "End Easting")


def calculation_rows(solution):
    for i, side in enumerate(solution.sides):
        start_e, start_n = solution.coordinates[i]
        end_e, end_n = solution.coordinates[i + 1]
        east, north = solution.components[i]
        yield (i + 1, side.length, quadrantal(side.azimuth), side.azimuth,
               north, east, start_n, start_e, end_n, end_e)


def report_rows(inputs, case, solutions, details, unit, selected=0):
    """One shared report model keeps TXT and CSV equally complete."""
    yield ["ENGINEERING TRAVERSE SOLVER", "Engineering 1035"]
    yield ["Detected Case", CASES[case]]
    yield ["Distance unit", unit, "Area unit", unit + "²"]
    yield ["Selected solution", selected + 1]
    yield ["Input Data"]
    yield ["Side", "Length", "Azimuth (°)"]
    for i, side in enumerate(inputs, 1):
        yield [i, "Unknown" if side.length is None else side.length,
               "Unknown" if side.azimuth is None else side.azimuth]
    for index, solution in enumerate(solutions, 1):
        yield []
        yield ["Solution", index]
        yield ["Solved Unknowns"]
        for i, (original, solved) in enumerate(zip(inputs, solution.sides), 1):
            if original.length is None:
                yield [f"Side {i} Length", solved.length, unit]
            if original.azimuth is None:
                yield [f"Side {i} Bearing", solved.azimuth, quadrantal(solved.azimuth)]
        yield ["Traverse Calculations"]
        yield list(TABLE_HEADERS)
        yield from calculation_rows(solution)
        yield ["Coordinates"]
        yield ["Point", "Northing", "Easting"]
        for i, (east, north) in enumerate(solution.coordinates, 1):
            yield [f"P{i}" if i <= len(inputs) else "Return to P1", north, east]
        yield ["Perimeter", solution.perimeter, unit]
        yield ["Area", solution.area if solution.area is not None else "Unavailable", unit + "²"]
        yield ["ΣΔN / Northing misclosure", solution.sum_n, unit]
        yield ["ΣΔE / Easting misclosure", solution.sum_e, unit]
        yield ["Closure Error", solution.closure, unit]
        yield ["Relative Precision", solution.precision]
        for warning in solution.warnings:
            yield ["Warning", warning]
    yield []
    yield ["CALCULATION DETAILS"]
    for line in details.splitlines():
        yield [line]
    yield ["Note", "Closure is imposed by the solver; numerical closure is not a measurement accuracy assessment."]
    if len(solutions) > 1:
        yield ["Ambiguity", "Additional engineering constraints are needed to select a solution."]


def export_report(path, inputs, case, solutions, details, unit, selected=0):
    rows = report_rows(inputs, case, solutions, details, unit, selected)
    with open(path, "w", encoding="utf-8-sig" if str(path).lower().endswith(".csv") else "utf-8", newline="") as stream:
        if str(path).lower().endswith(".csv"):
            csv.writer(stream).writerows(rows)
        else:
            for row in rows:
                stream.write("\t".join(str(value) for value in row) + "\n")


# Tk imports are deferred so the mathematical tests run on headless machines.
class TraverseApp:
    def __init__(self, root):
        import tkinter as tk
        from tkinter import ttk, font
        families = {name.lower(): name for name in font.families(root)}
        self.ui_font = next((families[name.lower()] for name in
                            ("Segoe UI", "DejaVu Sans", "Liberation Sans", "Helvetica")
                            if name.lower() in families), font.nametofont("TkDefaultFont").actual("family"))
        self.mono_font = font.nametofont("TkFixedFont").actual("family")
        self.tk, self.ttk, self.root = tk, ttk, root
        self.rows = []
        self.solutions = []
        self.inputs = []
        self.plot_window = None
        self.plot_canvas = None
        self.figure = None
        self.suspend = False
        self.last_format = "Azimuth"
        root.title("Engineering 1035 · Traverse Solver")
        root.geometry("1250x900")
        root.minsize(900, 680)
        root.configure(bg="#edf2f7")
        style = ttk.Style(root)
        style.theme_use("clam")
        style.configure("TFrame", background="#edf2f7")
        style.configure("TLabel", background="#edf2f7", foreground="#18314b", font=(self.ui_font, 10))
        style.configure("TButton", padding=(12, 7), font=(self.ui_font, 10))
        style.configure("Accent.TButton", background="#126d81", foreground="white")
        style.configure("Treeview", rowheight=28, font=(self.ui_font, 10))
        style.configure("Treeview.Heading", font=(self.ui_font, 10, "bold"))
        style.configure("Case.TLabel", font=(self.ui_font, 11, "bold"), foreground="#126d81")
        header = tk.Frame(root, bg="#18314b", padx=24, pady=15)
        header.pack(fill="x")
        tk.Label(header, text="ENGINEERING TRAVERSE SOLVER", bg="#18314b", fg="white",
                 font=(self.ui_font, 21, "bold")).pack(anchor="w")
        tk.Label(header, text="ENGINEERING 1035  /  Missing Length & Bearing Analysis", bg="#18314b",
                 fg="#a9dce5", font=(self.ui_font, 11)).pack(anchor="w", pady=(4, 0))
        main = ttk.Frame(root, padding=16)
        main.pack(fill="both", expand=True)
        controls = ttk.Frame(main)
        controls.pack(fill="x")
        self.count = tk.StringVar(value="4")
        self.format = tk.StringVar(value="Azimuth")
        self.unit = tk.StringVar(value="m")
        self.initial_n, self.initial_e = tk.StringVar(value="0"), tk.StringVar(value="0")
        for column, (label, variable, width) in enumerate((
            ("Number of sides", self.count, 7), ("Bearing format", self.format, 19),
            ("Initial Northing", self.initial_n, 14), ("Initial Easting", self.initial_e, 14),
            ("Distance unit", self.unit, 7))):
            cell = ttk.Frame(controls, padding=(0, 0, 14, 8))
            cell.grid(row=0, column=column, sticky="w")
            ttk.Label(cell, text=label).pack(anchor="w", pady=(0, 4))
            if variable is self.format or variable is self.unit:
                values = ("Azimuth", "Quadrantal Bearing") if variable is self.format else ("m", "ft")
                widget = ttk.Combobox(cell, textvariable=variable, values=values, state="readonly", width=width)
                widget.pack()
                if variable is self.format:
                    widget.bind("<<ComboboxSelected>>", self.change_format)
            elif variable is self.count:
                ttk.Spinbox(cell, textvariable=variable, from_=3, to=20, width=width).pack()
            else:
                ttk.Entry(cell, textvariable=variable, width=width).pack()
        ttk.Button(controls, text="Generate Traverse", command=self.generate).grid(row=0, column=5, padx=5)
        ttk.Label(main, text="Angles: 125.5 or 125°30'00\" (also 125:30:00). "
                  "Azimuth: 0° North, 90° East. Mark omitted values Unknown. Units label inputs; they do not convert them.").pack(anchor="w", pady=(2, 8))
        self.case_label = ttk.Label(main, style="Case.TLabel", wraplength=1160)
        self.case_label.pack(anchor="w", pady=(0, 8))
        root.bind("<Configure>", lambda e: self.case_label.configure(wraplength=max(600, e.width - 50))
                  if e.widget is root else None)
        table_frame = ttk.Frame(main)
        table_frame.pack(fill="x")
        self.input_canvas = tk.Canvas(table_frame, height=200, bg="white", highlightthickness=0)
        scroll = ttk.Scrollbar(table_frame, orient="vertical", command=self.input_canvas.yview)
        self.input_canvas.configure(yscrollcommand=scroll.set)
        self.input_canvas.pack(side="left", fill="both", expand=True)
        scroll.pack(side="right", fill="y")
        self.input_body = ttk.Frame(self.input_canvas, padding=8)
        self.input_window = self.input_canvas.create_window((0, 0), window=self.input_body, anchor="nw")
        self.input_body.bind("<Configure>", lambda _: self.input_canvas.configure(scrollregion=self.input_canvas.bbox("all")))
        self.input_canvas.bind("<Configure>", lambda e: self.input_canvas.itemconfigure(self.input_window, width=e.width))
        toolbar = ttk.Frame(main, padding=(0, 10))
        toolbar.pack(fill="x")
        ttk.Button(toolbar, text="Solve Traverse", style="Accent.TButton", command=self.solve).pack(side="left", padx=(0, 6))
        ttk.Button(toolbar, text="Clear", command=self.clear).pack(side="left", padx=6)
        self.example_case = tk.StringVar(value="Case I")
        ttk.Combobox(toolbar, textvariable=self.example_case, values=("Case I", "Case II", "Case III", "Case IV"),
                     state="readonly", width=9).pack(side="left", padx=(16, 4))
        ttk.Button(toolbar, text="Load Example", command=self.load_example).pack(side="left", padx=4)
        self.plot_button = ttk.Button(toolbar, text="Plot Traverse", command=self.plot, state="disabled")
        self.plot_button.pack(side="right", padx=(6, 0))
        self.export_button = ttk.Button(toolbar, text="Export Results", command=self.export, state="disabled")
        self.export_button.pack(side="right")
        solution_bar = ttk.Frame(main)
        solution_bar.pack(fill="x", pady=(0, 7))
        self.solution_label = ttk.Label(solution_bar, text="Solve a traverse to view its results.")
        self.solution_label.pack(side="left")
        self.solution_choice = ttk.Combobox(solution_bar, state="disabled", width=16)
        self.solution_choice.pack(side="right")
        self.solution_choice.bind("<<ComboboxSelected>>", lambda _: self.show_solution())
        tabs = ttk.Notebook(main)
        tabs.pack(fill="both", expand=True)
        self.summary = self.text_tab(tabs, "Results")
        calculations = ttk.Frame(tabs)
        tabs.add(calculations, text="Calculation Table")
        self.calculation_table = self.tree(calculations, TABLE_HEADERS)
        coordinates = ttk.Frame(tabs)
        tabs.add(coordinates, text="Coordinates")
        self.coordinate_table = self.tree(coordinates, ("Point", "Northing", "Easting"))
        self.details_box = self.text_tab(tabs, "Calculation Details")
        self.status = tk.StringVar(value="Ready · Load an example or enter 3–20 sides.")
        ttk.Label(main, textvariable=self.status).pack(anchor="w", pady=(8, 0))
        for variable in (self.initial_n, self.initial_e, self.unit, self.count):
            variable.trace_add("write", self.invalidate)
        self.generate()

    def text_tab(self, tabs, title):
        frame = self.ttk.Frame(tabs)
        tabs.add(frame, text=title)
        text = self.tk.Text(frame, wrap="word", height=10, width=60, font=(self.mono_font, 11), bg="white", fg="#18314b",
                            relief="flat", padx=16, pady=12, state="disabled")
        scrollbar = self.ttk.Scrollbar(frame, command=text.yview)
        text.configure(yscrollcommand=scrollbar.set)
        scrollbar.pack(side="right", fill="y")
        text.pack(fill="both", expand=True)
        return text

    def tree(self, frame, headers):
        tree = self.ttk.Treeview(frame, columns=headers, show="headings")
        for name in headers:
            tree.heading(name, text=name)
            tree.column(name, width=155 if name == "Bearing" else 135, anchor="e", stretch=False)
        vertical = self.ttk.Scrollbar(frame, command=tree.yview)
        horizontal = self.ttk.Scrollbar(frame, orient="horizontal", command=tree.xview)
        tree.configure(yscrollcommand=vertical.set, xscrollcommand=horizontal.set)
        tree.grid(row=0, column=0, sticky="nsew")
        vertical.grid(row=0, column=1, sticky="ns")
        horizontal.grid(row=1, column=0, sticky="ew")
        frame.rowconfigure(0, weight=1)
        frame.columnconfigure(0, weight=1)
        return tree

    @staticmethod
    def set_text(widget, value):
        widget.configure(state="normal")
        widget.delete("1.0", "end")
        widget.insert("1.0", value)
        widget.configure(state="disabled")

    def invalidate(self, *_):
        if self.suspend:
            return
        self.solutions = []
        self.plot_button.configure(state="disabled")
        self.export_button.configure(state="disabled")
        self.solution_choice.configure(state="disabled", values=())
        self.solution_choice.set("")
        self.solution_label.configure(text="Inputs changed — solve to refresh results.")
        self.set_text(self.summary, "Enter known values and mark exactly two quantities Unknown.")
        self.set_text(self.details_box, "")
        for tree in (self.calculation_table, self.coordinate_table):
            tree.delete(*tree.get_children())
        self.close_plot()
        self.update_case()

    def update_case(self):
        try:
            flags = [TraverseSide(None if r["ls"].get() == "Unknown" else 1,
                                  None if r["bs"].get() == "Unknown" else 0) for r in self.rows]
            self.case_label.configure(text="Detected: " + CASES[detect_case(flags)])
        except TraverseError as error:
            self.case_label.configure(text="Detected case: " + str(error))

    def generate(self):
        from tkinter import messagebox
        try:
            count = int(self.count.get())
            if not 3 <= count <= 20:
                raise ValueError
        except ValueError:
            messagebox.showerror("Invalid Input", "Number of sides must be an integer from 3 to 20.", parent=self.root)
            return
        self.suspend = True
        for child in self.input_body.winfo_children():
            child.destroy()
        self.rows = []
        for column, label in enumerate(("Side", "Length", "Length Status", "N/S", "Angle / Azimuth", "E/W", "Bearing Status")):
            self.ttk.Label(self.input_body, text=label, font=(self.ui_font, 10, "bold")).grid(row=0, column=column, sticky="w", padx=7, pady=5)
        self.input_body.columnconfigure(4, weight=1)
        for i in range(count):
            row = {name: self.tk.StringVar(value=value) for name, value in
                   (("length", ""), ("ls", "Known"), ("angle", ""), ("bs", "Known"), ("ns", "N"), ("ew", "E"))}
            self.ttk.Label(self.input_body, text=f"{i + 1:02d}").grid(row=i + 1, column=0, padx=7)
            widgets = [self.ttk.Entry(self.input_body, textvariable=row["length"], width=15),
                       self.ttk.Combobox(self.input_body, textvariable=row["ls"], values=("Known", "Unknown"), width=10, state="readonly"),
                       self.ttk.Combobox(self.input_body, textvariable=row["ns"], values=("N", "S"), width=4, state="readonly"),
                       self.ttk.Entry(self.input_body, textvariable=row["angle"], width=25),
                       self.ttk.Combobox(self.input_body, textvariable=row["ew"], values=("E", "W"), width=4, state="readonly"),
                       self.ttk.Combobox(self.input_body, textvariable=row["bs"], values=("Known", "Unknown"), width=10, state="readonly")]
            row["widgets"] = widgets
            for column, widget in enumerate(widgets, 1):
                widget.grid(row=i + 1, column=column, sticky="ew", padx=7, pady=3)
            self.rows.append(row)
            for name in ("length", "ls", "angle", "bs", "ns", "ew"):
                row[name].trace_add("write", lambda *_, r=row: self.row_changed(r))
            self.row_changed(row)
        self.suspend = False
        self.invalidate()
        self.input_canvas.yview_moveto(0)
        self.status.set(f"Generated {count} sides. Enter lengths in {self.unit.get()}.")

    def row_changed(self, row):
        widgets = row["widgets"]
        widgets[0].configure(state="disabled" if row["ls"].get() == "Unknown" else "normal")
        bearing_known = row["bs"].get() == "Known"
        widgets[3].configure(state="normal" if bearing_known else "disabled")
        for widget in (widgets[2], widgets[4]):
            widget.configure(state="readonly" if bearing_known and self.format.get() == "Quadrantal Bearing" else "disabled")
        self.invalidate()

    def change_format(self, *_):
        new_format = self.format.get()
        if new_format == self.last_format:
            return
        # Convert transactionally: an invalid row must not silently change its meaning.
        converted = []
        try:
            for i, row in enumerate(self.rows, 1):
                if not row["angle"].get().strip() or row["bs"].get() == "Unknown":
                    converted.append(None)
                    continue
                text = row["angle"].get()
                if self.last_format == "Quadrantal Bearing":
                    text = f'{row["ns"].get()} {text} {row["ew"].get()}'
                converted.append(parse_bearing(text))
        except TraverseError as error:
            from tkinter import messagebox
            self.format.set(self.last_format)
            messagebox.showerror("Invalid Bearing", f"Side {i}: {error}\nCorrect this bearing before switching formats.", parent=self.root)
            return
        self.suspend = True
        for row, angle in zip(self.rows, converted):
            if angle is not None:
                if new_format == "Quadrantal Bearing":
                    ns, value, ew = quadrant_parts(angle)
                    row["ns"].set(ns)
                    row["ew"].set(ew)
                    row["angle"].set(f"{value:.12g}")
                else:
                    row["angle"].set(f"{angle:.12g}")
            elif row["bs"].get() == "Unknown":
                row["angle"].set("")
            self.row_changed(row)
        self.last_format = new_format
        self.suspend = False
        self.invalidate()

    def clear(self):
        self.suspend = True
        self.initial_n.set("0")
        self.initial_e.set("0")
        for row in self.rows:
            for key, value in (("length", ""), ("angle", ""), ("ls", "Known"), ("bs", "Known"), ("ns", "N"), ("ew", "E")):
                row[key].set(value)
        self.suspend = False
        self.invalidate()
        self.status.set("Cleared · Enter a new traverse.")

    def load_example(self):
        index = ("Case I", "Case II", "Case III", "Case IV").index(self.example_case.get())
        self.count.set("4")
        self.generate()
        self.suspend = True
        self.initial_n.set("0")
        self.initial_e.set("0")
        for row, side in zip(self.rows, example(index)):
            row["ls"].set("Unknown" if side.length is None else "Known")
            row["length"].set("" if side.length is None else str(side.length))
            row["bs"].set("Unknown" if side.azimuth is None else "Known")
            if side.azimuth is not None:
                ns, angle, ew = quadrant_parts(side.azimuth)
                row["ns"].set(ns)
                row["ew"].set(ew)
                row["angle"].set(str(angle if self.format.get() == "Quadrantal Bearing" else side.azimuth))
        self.suspend = False
        self.invalidate()
        self.status.set(f"Loaded {self.example_case.get()}: 80 × 100 rectangle. Click Solve Traverse.")

    def read_inputs(self):
        if self.count.get().strip() != str(len(self.rows)):
            raise TraverseError("Click Generate Traverse to apply the changed number of sides.")
        sides = []
        for i, row in enumerate(self.rows, 1):
            try:
                length = None if row["ls"].get() == "Unknown" else number(row["length"].get(), "Length")
                angle = None
                if row["bs"].get() == "Known":
                    text = row["angle"].get()
                    if self.format.get() == "Quadrantal Bearing":
                        text = f'{row["ns"].get()} {text} {row["ew"].get()}'
                    angle = parse_bearing(text)
                sides.append(TraverseSide(length, angle))
            except TraverseError as error:
                raise TraverseError(f"Side {i}: {error}") from None
        return sides

    def solve(self):
        from tkinter import messagebox
        try:
            inputs = self.read_inputs()
            case, solutions, details = TraverseSolver().solve(inputs,
                number(self.initial_e.get(), "Initial Easting"), number(self.initial_n.get(), "Initial Northing"))
        except (TraverseError, OverflowError, ArithmeticError) as error:
            self.status.set("Cannot solve · Check inputs and geometric constraints.")
            messagebox.showerror("Cannot Solve", str(error), parent=self.root)
            return
        self.inputs, self.case, self.solutions, self.details = inputs, case, solutions, details
        self.solution_choice.configure(values=[f"Solution {i + 1}" for i in range(len(solutions))], state="readonly")
        self.solution_choice.current(0)
        self.solution_label.configure(text=(f"{len(solutions)} valid solutions — use engineering constraints to choose."
                                           if len(solutions) > 1 else "One valid solution · Closure verified"))
        self.plot_button.configure(state="normal")
        self.export_button.configure(state="normal")
        self.show_solution()
        self.status.set(f"Solved · {len(solutions)} solution(s), each verified against both closure equations.")

    def show_solution(self):
        if not self.solutions:
            return
        index = max(0, self.solution_choice.current())
        solution = self.solutions[index]
        unit = self.unit.get()
        lines = [CASES[self.case], f"Solution {index + 1} of {len(self.solutions)}", "", "SOLVED VALUES"]
        for i, (original, solved) in enumerate(zip(self.inputs, solution.sides), 1):
            if original.length is None:
                lines.append(f"Side {i} Length: {solved.length:.4f} {unit}")
            if original.azimuth is None:
                lines.append(f"Side {i} Bearing: {display_azimuth(solved.azimuth)}°  /  {quadrantal(solved.azimuth)}")
        area = f"{solution.area:,.4f} {unit}²" if solution.area is not None else "Unavailable — see warning"
        lines += ["", "TRAVERSE SUMMARY", f"Perimeter: {solution.perimeter:,.4f} {unit}",
                  f"Traverse Area: {area}", f"ΣΔE / Easting misclosure: {solution.sum_e:.10g} {unit}",
                  f"ΣΔN / Northing misclosure: {solution.sum_n:.10g} {unit}",
                  f"Linear Closure Error: {solution.closure:.10g} {unit}",
                  f"Relative Precision: {solution.precision}", "",
                  "Closure is imposed by solving the omitted values. Numerical closure does not measure survey accuracy."]
        if len(self.solutions) > 1:
            lines += ["", "MULTIPLE SOLUTIONS: Review each solution using the selector above. "
                      "Both satisfy the supplied equations; extra engineering constraints are needed to choose."]
        lines += ["\nWARNING: " + warning for warning in solution.warnings]
        self.set_text(self.summary, "\n".join(lines))
        self.set_text(self.details_box, self.details)
        self.calculation_table.delete(*self.calculation_table.get_children())
        for row in calculation_rows(solution):
            values = [f"{value:.4f}" if isinstance(value, float) else value for value in row]
            values[3] = display_azimuth(row[3])
            self.calculation_table.insert("", "end", values=values)
        self.coordinate_table.delete(*self.coordinate_table.get_children())
        for i, (east, north) in enumerate(solution.coordinates):
            name = f"P{i + 1}" if i < len(solution.sides) else "Return to P1"
            self.coordinate_table.insert("", "end", values=(name, f"{north:.4f}", f"{east:.4f}"))
        if self.plot_window is not None:
            self.draw_plot()

    def close_plot(self):
        if self.plot_window is not None:
            self.plot_window.destroy()
            self.plot_window = None
        self.plot_canvas = None
        self.figure = None

    def plot(self):
        if not self.solutions:
            return
        if self.plot_window is None:
            self.plot_window = self.tk.Toplevel(self.root)
            self.plot_window.title("Solved Closed Traverse")
            self.plot_window.geometry("780x620")
            self.plot_window.protocol("WM_DELETE_WINDOW", self.close_plot)
            try:
                from matplotlib.figure import Figure
                from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
            except ImportError:
                self.ttk.Label(self.plot_window, text="Matplotlib unavailable — using built-in plot. "
                               "Install matplotlib for zoom and image export.").pack(fill="x")
                self.plot_canvas = self.tk.Canvas(self.plot_window, bg="white", width=760, height=560)
                self.plot_canvas.pack(fill="both", expand=True)
                self.plot_canvas.bind("<Configure>", lambda _: self.draw_plot())
            else:
                self.figure = Figure(figsize=(7.8, 5.6), dpi=100, layout="constrained")
                self.plot_canvas = FigureCanvasTkAgg(self.figure, master=self.plot_window)
                self.plot_canvas.get_tk_widget().pack(fill="both", expand=True)
                NavigationToolbar2Tk(self.plot_canvas, self.plot_window)
        self.draw_plot()
        self.plot_window.lift()

    def draw_plot(self):
        if not self.solutions or self.plot_canvas is None:
            return
        solution = self.solutions[max(0, self.solution_choice.current())]
        # Include the computed return point, then explicitly close to the start.
        points = solution.coordinates + [solution.coordinates[0]]
        if self.figure is not None:
            self.figure.clear()
            axes = self.figure.add_subplot(111)
            axes.plot([e for e, _ in points], [n for _, n in points], "o-", color="#126d81", linewidth=2)
            axes.scatter(*points[0], color="#d35e35", s=100, zorder=4, label="Start / P1")
            for i, (east, north) in enumerate(points[:-2], 1):
                axes.annotate(f"P{i}", (east, north), xytext=(7, 7), textcoords="offset points")
            for start, end in zip(points, points[1:]):
                axes.annotate("", xy=end, xytext=start, arrowprops={"arrowstyle": "->", "color": "#126d81"})
            axes.set(title=f"Solved Closed Traverse · Solution {self.solution_choice.current() + 1}",
                     xlabel=f"Easting ({self.unit.get()})", ylabel=f"Northing ({self.unit.get()})")
            axes.set_aspect("equal", adjustable="datalim")
            axes.grid(True, alpha=0.3)
            axes.legend()
            self.plot_canvas.draw()
        else:
            canvas = self.plot_canvas
            canvas.delete("all")
            width, height = max(canvas.winfo_width(), 300), max(canvas.winfo_height(), 300)
            min_e, max_e = min(e for e, _ in points), max(e for e, _ in points)
            min_n, max_n = min(n for _, n in points), max(n for _, n in points)
            scale = min((width - 120) / max(max_e - min_e, 1e-100), (height - 120) / max(max_n - min_n, 1e-100))
            center_e, center_n = (min_e + max_e) / 2, (min_n + max_n) / 2
            mapped = [(width / 2 + (e - center_e) * scale, height / 2 - (n - center_n) * scale) for e, n in points]
            for x in range(40, width, 40):
                canvas.create_line(x, 40, x, height - 40, fill="#edf2f7")
            for y in range(40, height, 40):
                canvas.create_line(40, y, width - 40, y, fill="#edf2f7")
            for a, b in zip(mapped, mapped[1:]):
                canvas.create_line(*a, *b, fill="#126d81", width=2, arrow="last")
            for i, (x, y) in enumerate(mapped[:-2]):
                canvas.create_oval(x - 4, y - 4, x + 4, y + 4, fill="#d35e35" if i == 0 else "#126d81")
                canvas.create_text(x + 12, y - 12, text=f"P{i + 1}")
            canvas.create_text(width / 2, 20, text=f"Solved Closed Traverse · E → / N ↑ ({self.unit.get()})")
            canvas.create_text(width / 2, height - 18, text=f"Easting: {min_e:.4f} … {max_e:.4f}   "
                                      f"Northing: {min_n:.4f} … {max_n:.4f} · Equal scale")

    def export(self):
        from tkinter import filedialog, messagebox
        if not self.solutions:
            return
        path = filedialog.asksaveasfilename(parent=self.root, title="Export all solutions", defaultextension=".txt",
                                          initialfile="traverse_results.txt", filetypes=(("Text report", "*.txt"), ("CSV report", "*.csv")))
        if not path:
            return
        try:
            export_report(path, self.inputs, self.case, self.solutions, self.details,
                          self.unit.get(), max(0, self.solution_choice.current()))
        except OSError as error:
            messagebox.showerror("Export Failed", str(error), parent=self.root)
            return
        self.status.set(f"Exported all {len(self.solutions)} solution(s) to {Path(path).name}.")


class SolverTests(unittest.TestCase):
    def test_all_examples(self):
        for case in range(4):
            with self.subTest(case=case):
                detected, solutions, details = TraverseSolver().solve(example(case), 1500, 2200)
                self.assertEqual(detected, case)
                self.assertIn("Closure equations", details)
                self.assertTrue(solutions)
                for solution in solutions:
                    self.assertLess(solution.closure, 1e-8)
                    self.assertAlmostEqual(solution.perimeter, 360)
                    self.assertAlmostEqual(solution.coordinates[-1][0], 1500)
                    self.assertAlmostEqual(solution.coordinates[-1][1], 2200)
                self.assertAlmostEqual(solutions[0].area, 8000)

    def test_bearings_and_dms(self):
        for text, expected in (("N 30° E", 30), ("S 45 E", 135), ("S 20 W", 200),
                               ("N 70 W", 290), ("125°30'00\"", 125.5), ("125:30:00", 125.5),
                               ("125 30 00", 125.5), ("360°", 0), ("N 0 W", 0)):
            self.assertAlmostEqual(parse_bearing(text), expected)
        for text in ("-1", "361", "nan", "inf", "N 91 E", "N 30 N", "20°60'00\"", "20:10:60", "", "1e3"):
            with self.subTest(text=text), self.assertRaises(TraverseError):
                parse_bearing(text)
        self.assertEqual(azimuth(-1e-20, 1), 0)
        self.assertEqual(display_azimuth(359.999999999), "0.0000")
        self.assertEqual(dms(12 + 59 / 60 + 59.999 / 3600), "13°00'00.00\"")
        for angle in range(360):
            self.assertAlmostEqual(parse_bearing(quadrantal(angle)), angle)

    def test_case_ii_two_positive_roots(self):
        # target=(0, 5), unknown northward length, other length=2 -> x=3 or 7.
        sides = [TraverseSide(None, 0), TraverseSide(2, None), TraverseSide(5, 180)]
        _, solutions, _ = TraverseSolver().solve(sides)
        self.assertEqual(len(solutions), 2)
        self.assertEqual({round(s.sides[0].length) for s in solutions}, {3, 7})

    def test_case_ii_no_intersection_and_tangent(self):
        with self.assertRaises(TraverseError):
            TraverseSolver().solve([TraverseSide(None, 0), TraverseSide(2, None), TraverseSide(5, 270)])
        _, solutions, _ = TraverseSolver().solve([TraverseSide(None, 0), TraverseSide(5, None),
                                                  TraverseSide(math.sqrt(34), azimuth(-5, -3))])
        self.assertEqual(len(solutions), 1)
        self.assertAlmostEqual(solutions[0].sides[0].length, 3, places=6)

    def test_case_iv_mirrors_and_tangent(self):
        _, solutions, _ = TraverseSolver().solve(example(3))
        self.assertEqual(len(solutions), 2)
        sides = [TraverseSide(3, None), TraverseSide(2, None), TraverseSide(5, 180)]
        self.assertEqual(len(TraverseSolver().solve(sides)[1]), 1)

    def test_impossible_singular_and_indeterminate(self):
        datasets = [
            [TraverseSide(None, 0), TraverseSide(None, 180), TraverseSide(3, 90)],
            [TraverseSide(1, None), TraverseSide(1, None), TraverseSide(5, 90)],
            [TraverseSide(1, None), TraverseSide(1, None), TraverseSide(3, 0), TraverseSide(3, 180)],
            [TraverseSide(None, None), TraverseSide(3, 0), TraverseSide(3, 180)],
            [TraverseSide(None, 0), TraverseSide(None, 90), TraverseSide(5, 45)],
            [TraverseSide(None, None), TraverseSide(-2, 90), TraverseSide(3, 180)],
            [TraverseSide(None, None), TraverseSide(float("nan"), 90), TraverseSide(3, 180)],
        ]
        for sides in datasets:
            with self.subTest(sides=sides), self.assertRaises(TraverseError):
                TraverseSolver().solve(sides)
        for sides in (example(0)[:2], [TraverseSide(1, 0)] * 4, [TraverseSide(None, None)] * 4):
            with self.assertRaises(TraverseError):
                TraverseSolver().solve(sides)

    def test_area_and_large_coordinate_offset(self):
        _, solutions, _ = TraverseSolver().solve(example(0), 1e12, -1e12)
        self.assertAlmostEqual(solutions[0].area, 8000)
        area, warnings = polygon_area([(0, 0), (2, 2), (0, 2), (2, 0), (0, 0)], 1e-9)
        self.assertIsNone(area)
        self.assertTrue(warnings)

    def test_report_export(self):
        import tempfile
        sides = example(3)
        case, solutions, details = TraverseSolver().solve(sides)
        with tempfile.TemporaryDirectory() as directory:
            for extension in ("txt", "csv"):
                path = Path(directory) / f"report.{extension}"
                export_report(path, sides, case, solutions, details, "ft")
                content = path.read_text(encoding="utf-8-sig")
                for required in ("Input Data", "Solved Unknowns", "Coordinates", "Perimeter", "Area",
                                 "Closure Error", "Relative Precision", "CALCULATION DETAILS", "ft²"):
                    self.assertIn(required, content)
                if extension == "csv":
                    rows = list(csv.reader(io.StringIO(content)))
                    self.assertIn(["Solution", "2"], rows)

    def test_rotated_scaled_traverses(self):
        # A nonrectangular pentagon exercises arbitrary bearings and all unknown placements.
        points = [(0, 0), (3, 5), (9, 6), (12, 2), (5, -3)]
        complete = []
        for first, second in zip(points, points[1:] + points[:1]):
            e, n = second[0] - first[0], second[1] - first[1]
            complete.append(TraverseSide(math.hypot(e, n), azimuth(e, n)))
        for scale in (1e-5, 1, 1e7):
            for rotation in (0, 33, 181, 299):
                sides = [TraverseSide(s.length * scale, (s.azimuth + rotation) % 360) for s in complete]
                for case in range(4):
                    for i in range(5):
                        j = (i + 1) % 5
                        inputs = list(sides)
                        if case == 0:
                            inputs[i] = TraverseSide(None, None)
                        elif case == 1:
                            inputs[i] = replace(inputs[i], length=None)
                            inputs[j] = replace(inputs[j], azimuth=None)
                        elif case == 2:
                            inputs[i] = replace(inputs[i], length=None)
                            inputs[j] = replace(inputs[j], length=None)
                        else:
                            inputs[i] = replace(inputs[i], azimuth=None)
                            inputs[j] = replace(inputs[j], azimuth=None)
                        _, solutions, _ = TraverseSolver().solve(inputs)
                        self.assertTrue(any(TraverseSolver._same(s.sides, sides) for s in solutions))
                        for solution in solutions:
                            self.assertLessEqual(solution.closure, solution.tolerance)


def gui_smoke_test(root, app):
    """Exercise actual widgets, both bearing modes, plotting, and export dialogs."""
    from unittest.mock import patch
    import tempfile
    for mode in ("Azimuth", "Quadrantal Bearing"):
        app.format.set(mode)
        app.change_format()
        for case_name in ("Case I", "Case II", "Case III", "Case IV"):
            app.example_case.set(case_name)
            app.load_example()
            app.solve()
            assert app.solutions, (mode, case_name)
            assert len(app.calculation_table.get_children()) == 4
            assert len(app.coordinate_table.get_children()) == 5
            app.plot()
            root.update()
            if len(app.solutions) > 1:
                app.solution_choice.current(1)
                app.show_solution()
                root.update()
            with tempfile.TemporaryDirectory() as directory:
                for extension in ("txt", "csv"):
                    target = str(Path(directory) / f"gui.{extension}")
                    with patch("tkinter.filedialog.asksaveasfilename", return_value=target):
                        app.export()
                    assert Path(target).stat().st_size > 100
            app.close_plot()
    app.count.set("20")
    app.generate()
    assert len(app.rows) == 20
    app.clear()
    assert not app.solutions
    with patch("tkinter.messagebox.showerror") as error:
        app.solve()
        assert error.called
    # Exercise the built-in plot fallback even when Matplotlib is installed.
    app.example_case.set("Case I")
    app.load_example()
    app.solve()
    with patch.dict(sys.modules, {"matplotlib.figure": None}):
        app.plot()
        root.update()
        assert app.figure is None
    app.close_plot()
    print("GUI smoke test passed: four cases, two bearing formats, solution selection, "
          "Matplotlib and fallback plots, TXT/CSV dialogs, 20 rows, clear, invalid input.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true", help="Run embedded mathematical regression tests")
    parser.add_argument("--gui-test", action="store_true", help="Run GUI smoke test (requires a display)")
    args = parser.parse_args()
    if args.self_test:
        suite = unittest.defaultTestLoader.loadTestsFromTestCase(SolverTests)
        return 0 if unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful() else 1
    try:
        import tkinter as tk
    except ImportError:
        print("Tkinter is required. Install your Python distribution's Tk support (Ubuntu: python3-tk).", file=sys.stderr)
        return 1
    try:
        root = tk.Tk()
    except tk.TclError as error:
        print(f"Cannot open the desktop display: {error}\nRun this app in a graphical desktop session. "
              "Use --self-test for headless numerical validation.", file=sys.stderr)
        return 1
    app = TraverseApp(root)
    if args.gui_test:
        try:
            gui_smoke_test(root, app)
        finally:
            root.destroy()
    else:
        root.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
