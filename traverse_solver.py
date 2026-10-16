#!/usr/bin/env python3
"""Engineering 1035 omitted-measurement solver. Run: python traverse_solver.py.

The desktop app is distributed as one file. Matplotlib enables solved plots;
without it, numerical solving, input persistence, and report export still work.
Tests live separately in test_traverse_solver.py.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, replace
import csv
import json
from pathlib import Path
import math
import re
import sys

# ----------------------- Parsing and numerical core -----------------------

CASES = (
    "Case I: Length & Bearing of One Side Omitted",
    "Case II: Length of One Side & Bearing of Another Omitted",
    "Case III: Length of Two Sides Omitted",
    "Case IV: Bearing of Two Sides Omitted",
)


class TraverseError(ValueError):
    """An input or geometric problem that the student can correct."""


_DECIMAL = r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?"


def number(text, label="Value"):
    token = str(text).strip()
    if not re.fullmatch(_DECIMAL, token):
        raise TraverseError(f"{label} must be an ASCII decimal number (no underscores).")
    value = float(token)
    if not math.isfinite(value):
        raise TraverseError(f"{label} must be finite.")
    return value


def parse_angle(text):
    text = str(text).strip().replace("º", "°").replace("′", "'")
    text = text.replace("″", '\"').replace("''", '\"').replace("’", "'")
    if re.fullmatch(_DECIMAL + r"\s*°?", text):
        return number(text.rstrip("°").strip(), "Angle")
    patterns = (
        r"([0-9]+)\s*°\s*([0-9]+(?:\.[0-9]*)?)\s*'\s*"
        r"(?:([0-9]+(?:\.[0-9]*)?)\s*\")?",
        r"([0-9]+)\s*[:\s]\s*([0-9]+(?:\.[0-9]*)?)\s*[:\s]\s*"
        r"([0-9]+(?:\.[0-9]*)?)",
    )
    for pattern in patterns:
        match = re.fullmatch(pattern, text)
        if match:
            degrees, minutes, seconds = (float(part or 0) for part in match.groups())
            if not (minutes < 60 and seconds < 60):
                raise TraverseError("DMS minutes and seconds must be in [0, 60).")
            return degrees + minutes / 60 + seconds / 3600
    raise TraverseError("Use decimal degrees (125.5) or DMS (125°30'00\").")


def parse_bearing(text):
    text = str(text).strip().upper()
    match = re.fullmatch(r"([NS])\s*(.*?)\s*([EW])", text)
    if match:
        ns, angle_text, ew = match.groups()
        angle = parse_angle(angle_text)
        if not 0 <= angle <= 90:
            raise TraverseError("Quadrantal angles must be between 0° and 90°.")
        return {("N", "E"): angle, ("S", "E"): 180 - angle,
                ("S", "W"): 180 + angle, ("N", "W"): 360 - angle}[ns, ew] % 360
    angle = parse_angle(text)
    if not 0 <= angle <= 360:
        raise TraverseError("Azimuth must be between 0° and 360°.")
    return angle % 360


def quadrantal(angle):
    # Round before choosing a quadrant so carrying seconds gives canonical cardinals.
    ticks = round((angle % 360) * 360000) % 129600000
    if ticks <= 32400000:
        ns, ticks, ew = "N", ticks, "E"
    elif ticks <= 64800000:
        ns, ticks, ew = "S", 64800000 - ticks, "E"
    elif ticks < 97200000:
        ns, ticks, ew = "S", ticks - 64800000, "W"
    else:
        ns, ticks, ew = "N", 129600000 - ticks, "W"
    degrees, rest = divmod(ticks, 360000)
    minutes, seconds = divmod(rest, 6000)
    return f"{ns} {degrees}°{minutes:02d}'{seconds / 100:05.2f}\" {ew}"


def display_azimuth(angle):
    return f"{round(angle % 360, 4) % 360:.4f}"


@dataclass(frozen=True)
class TraverseSide:
    length: float | None
    azimuth: float | None


@dataclass(frozen=True)
class NumericalPolicy:
    scale: float
    relative: float = 1e-10

    @property
    def length(self):
        return self.scale * self.relative

    @property
    def squared_length(self):
        return self.scale * self.length

    @property
    def roundoff_squared(self):
        return 64 * sys.float_info.epsilon * self.scale**2

    # Unit-vector determinants are dimensionless; length tolerances do not apply.
    singular = 64 * sys.float_info.epsilon
    conditioning_warning = 0.01


@dataclass
class Solution:
    sides: list[TraverseSide]
    components: list[tuple[float, float]]
    local_coordinates: list[tuple[float, float]]
    coordinates: list[tuple[float, float]]
    perimeter: float
    sum_e: float
    sum_n: float
    closure: float
    area: float | None
    warnings: list[str]
    tolerance: float

    @property
    def geometry_label(self):
        return "Simple polygon" if self.area is not None else "Crossed or degenerate"


@dataclass
class SolveResult:
    case: int
    solutions: list[Solution]
    details: str
    condition_number: float | None = None


def case_from_unknowns(length_ids, bearing_ids):
    lengths, bearings = set(length_ids), set(bearing_ids)
    if len(lengths) == len(bearings) == 1:
        return 0 if lengths == bearings else 1
    if len(lengths) == 2 and not bearings:
        return 2
    if len(bearings) == 2 and not lengths:
        return 3
    raise TraverseError("Mark exactly two unknowns: one length and one bearing, "
                        "two lengths, or two bearings.")


def detect_case(sides):
    return case_from_unknowns([i for i, s in enumerate(sides) if s.length is None],
                              [i for i, s in enumerate(sides) if s.azimuth is None])


def unit_vector(angle):
    return math.sin(math.radians(angle)), math.cos(math.radians(angle))


def vector(side):
    east, north = unit_vector(side.azimuth)
    return side.length * east, side.length * north


def vector_azimuth(east, north):
    angle = math.degrees(math.atan2(east, north)) % 360
    return 0.0 if angle == 360 else angle


def polygon_area(points, policy):
    """Shoelace area, or None for intersecting, touching, or retraced boundaries."""
    vertices = points[:-1]
    count = len(vertices)
    tol, cross_tol = policy.length, policy.squared_length

    def cross(a, b, c):
        return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])

    def on_segment(a, b, p):
        return (abs(cross(a, b, p)) <= cross_tol
                and min(a[0], b[0]) - tol <= p[0] <= max(a[0], b[0]) + tol
                and min(a[1], b[1]) - tol <= p[1] <= max(a[1], b[1]) + tol)

    def opposite(a, b):
        return (a > cross_tol and b < -cross_tol) or (a < -cross_tol and b > cross_tol)

    for i, a in enumerate(vertices):
        b, previous = vertices[(i + 1) % count], vertices[i - 1]
        retraces = sum((p - q) * (r - q) for p, q, r in zip(previous, a, b)) > tol**2
        if abs(cross(previous, a, b)) <= cross_tol and retraces:
            return None
        for j in range(i + 1, count):
            if j == i + 1 or (i == 0 and j == count - 1):
                continue
            c, d = vertices[j], vertices[(j + 1) % count]
            crosses = opposite(cross(a, b, c), cross(a, b, d))
            crosses &= opposite(cross(c, d, a), cross(c, d, b))
            if crosses or any((on_segment(a, b, c), on_segment(a, b, d),
                               on_segment(c, d, a), on_segment(c, d, b))):
                return None
    area = abs(math.fsum(a[0] * b[1] - b[0] * a[1]
                         for a, b in zip(vertices, vertices[1:] + vertices[:1]))) / 2
    return area if area > cross_tol else None


def solve_case_i(sides, target, policy):
    i = next(i for i, s in enumerate(sides) if s.length is None)
    length = math.hypot(*target)
    if length <= policy.length:
        raise TraverseError("The omitted side has zero length; its bearing is indeterminate.")
    result = list(sides)
    result[i] = TraverseSide(length, vector_azimuth(*target))
    return [result]


def solve_case_ii(sides, target, policy):
    i = next(i for i, s in enumerate(sides) if s.length is None)
    j = next(i for i, s in enumerate(sides) if s.azimuth is None)
    east, north = target
    ue, un = unit_vector(sides[i].azimuth)
    projection, perpendicular = east * ue + north * un, east * un - north * ue
    discriminant = sides[j].length**2 - perpendicular**2
    if discriminant < -policy.squared_length:
        raise TraverseError("No solution: the known length cannot reach the known-bearing ray.")
    if abs(discriminant) <= policy.roundoff_squared:
        discriminant = 0.0
    root = math.sqrt(max(0, discriminant))
    candidates = []
    for length in (projection + root, projection - root):
        if length <= policy.length:
            continue
        candidate = list(sides)
        candidate[i] = replace(sides[i], length=length)
        candidate[j] = replace(sides[j], azimuth=vector_azimuth(
            east - length * ue, north - length * un))
        candidates.append(candidate)
    return candidates


def _length_system(sides):
    i, j = [i for i, s in enumerate(sides) if s.length is None]
    ae, an = unit_vector(sides[i].azimuth)
    be, bn = unit_vector(sides[j].azimuth)
    return i, j, ae, an, be, bn, ae * bn - be * an


def solve_case_iii(sides, target, policy):
    i, j, ae, an, be, bn, determinant = _length_system(sides)
    if abs(determinant) <= policy.singular:
        raise TraverseError("Singular system: the two bearings are parallel or opposite; "
                            "the missing lengths cannot be uniquely determined.")
    east, north = target
    first = (east * bn - be * north) / determinant
    second = (ae * north - east * an) / determinant
    candidate = list(sides)
    candidate[i], candidate[j] = replace(sides[i], length=first), replace(sides[j], length=second)
    return [candidate]


def solve_case_iv(sides, target, policy):
    i, j = [i for i, s in enumerate(sides) if s.azimuth is None]
    a, b, distance = sides[i].length, sides[j].length, math.hypot(*target)
    if distance <= policy.length:
        if abs(a - b) <= policy.length:
            raise TraverseError("Infinitely many solutions: equal opposite sides can rotate "
                                "freely; supply another bearing or engineering constraint.")
        raise TraverseError("No solution: unequal lengths cannot sum to a zero vector.")
    if distance > a + b + policy.length or distance < abs(a - b) - policy.length:
        raise TraverseError("No solution: known lengths violate the triangle inequality.")
    along = (distance**2 + a**2 - b**2) / (2 * distance)
    height_sq = a**2 - along**2
    if height_sq < -policy.squared_length:
        raise TraverseError("No stable circle intersection exists for these inputs.")
    if abs(height_sq) <= policy.roundoff_squared:
        height_sq = 0.0
    height, ue, un = math.sqrt(max(0, height_sq)), target[0] / distance, target[1] / distance
    candidates = []
    for sign in (1, -1):
        east, north = along * ue - sign * height * un, along * un + sign * height * ue
        candidate = list(sides)
        candidate[i] = replace(sides[i], azimuth=vector_azimuth(east, north))
        candidate[j] = replace(sides[j], azimuth=vector_azimuth(
            target[0] - east, target[1] - north))
        candidates.append(candidate)
    return candidates


def evaluate_solution(sides, initial_e, initial_n):
    components = [vector(side) for side in sides]
    perimeter = math.fsum(side.length for side in sides)
    policy = NumericalPolicy(perimeter)
    local = [(0.0, 0.0)] + [(math.fsum(e for e, _ in components[:i]),
                              math.fsum(n for _, n in components[:i]))
                             for i in range(1, len(sides) + 1)]
    coordinates = [(initial_e + e, initial_n + n) for e, n in local]
    if any(not math.isfinite(value) for point in coordinates for value in point):
        raise TraverseError("Coordinates exceed the supported floating-point range.")
    area = polygon_area(local, policy)
    warnings = [] if area is not None else [
        "Crossed or degenerate boundary: sides intersect, touch, or retrace; area unavailable."]
    sum_e, sum_n = local[-1]
    return Solution(sides, components, local, coordinates, perimeter, sum_e, sum_n,
                    math.hypot(sum_e, sum_n), area, warnings, policy.length)


def explain_result(sides, case, target, solutions, condition_number=None):
    east, north = target
    lines = [CASES[case], "Azimuth is clockwise from North.",
             "ΔE = L sin(θ); ΔN = L cos(θ); ΣΔE = 0; ΣΔN = 0.",
             f"Unknown vectors must sum to E={east:.12g}, N={north:.12g}."]
    length_ids = [i for i, side in enumerate(sides) if side.length is None]
    bearing_ids = [i for i, side in enumerate(sides) if side.azimuth is None]
    if case == 0:
        lines += ["L = hypot(E, N); θ = atan2(E, N) modulo 360°."]
    elif case == 1:
        i, j = length_ids[0], bearing_ids[0]
        ue, un = unit_vector(sides[i].azimuth)
        radius = sides[j].length
        projection = east * ue + north * un
        lines += [f"{ue:.12g}x + {radius:.12g}sin(φ) = {east:.12g}",
                  f"{un:.12g}x + {radius:.12g}cos(φ) = {north:.12g}",
                  f"x² - 2({projection:.12g})x + {east**2 + north**2 - radius**2:.12g} = 0.",
                  "Retain every positive root; recover φ from the remaining vector."]
    elif case == 2:
        i, j, ae, an, be, bn, determinant = _length_system(sides)
        lines += [f"{ae:.12g}L{i+1} + {be:.12g}L{j+1} = {east:.12g}",
                  f"{an:.12g}L{i+1} + {bn:.12g}L{j+1} = {north:.12g}",
                  f"Determinant = {determinant:.12g}; solve by Cramer's rule.",
                  f"2-norm condition number = {condition_number:.6g}."]
    else:
        i, j = bearing_ids
        lines += [f"{sides[i].length:.12g}sin(θ{i+1}) + "
                  f"{sides[j].length:.12g}sin(θ{j+1}) = {east:.12g}",
                  f"{sides[i].length:.12g}cos(θ{i+1}) + "
                  f"{sides[j].length:.12g}cos(θ{j+1}) = {north:.12g}",
                  "R = hypot(E,N); x = (R²+a²-b²)/(2R); h = sqrt(a²-x²).",
                  "First vector = x·unit(target) ± h·perpendicular(unit(target));",
                  "second vector = target − first. Convert both using atan2(E,N)."]
    for index, solution in enumerate(solutions, 1):
        lines += ["", f"Solution {index}: {solution.geometry_label}"]
        for i, (original, solved) in enumerate(zip(sides, solution.sides), 1):
            if original.length is None:
                lines.append(f"Side {i} length = {solved.length:.12g}")
            if original.azimuth is None:
                lines.append(f"Side {i} bearing = {solved.azimuth:.12g}° "
                             f"({quadrantal(solved.azimuth)})")
        lines += [f"Substitution residual: ΣΔE={solution.sum_e:.12g}; "
                  f"ΣΔN={solution.sum_n:.12g}; norm={solution.closure:.12g}.", *solution.warnings]
    lines += ["", "Closure is imposed by solving the equations. The residual checks numerical "
              "substitution only; it does not assess measurement quality."]
    return "\n".join(lines)


def solve_traverse(sides, initial_e=0.0, initial_n=0.0):
    if not 3 <= len(sides) <= 20:
        raise TraverseError("Number of sides must be from 3 to 20.")
    initial_e = number(initial_e, "Initial Easting")
    initial_n = number(initial_n, "Initial Northing")
    for index, side in enumerate(sides, 1):
        if side.length is not None and (not math.isfinite(side.length) or side.length <= 0):
            raise TraverseError(f"Side {index}: length must be a positive finite number.")
        if side.azimuth is not None and (not math.isfinite(side.azimuth)
                                         or not 0 <= side.azimuth < 360):
            raise TraverseError(f"Side {index}: azimuth must be normalized to [0°, 360°).")
    case = detect_case(sides)
    known = [vector(s) for s in sides if s.length is not None and s.azimuth is not None]
    target = -math.fsum(e for e, _ in known), -math.fsum(n for _, n in known)
    policy = NumericalPolicy(math.fsum(s.length or 0 for s in sides))
    candidates = (solve_case_i, solve_case_ii, solve_case_iii, solve_case_iv)[case](
        sides, target, policy)
    solutions, condition_number = [], None
    if case == 2:
        _, _, ae, an, be, bn, determinant = _length_system(sides)
        condition_number = (1 + abs(ae * be + an * bn)) / abs(determinant)
    for candidate in candidates:
        if any(s.length <= policy.length for s in candidate):
            continue
        solution = evaluate_solution(candidate, initial_e, initial_n)
        duplicate = any(all(math.hypot(a[0] - b[0], a[1] - b[1]) <= solution.tolerance
                            for a, b in zip(solution.components, old.components))
                        for old in solutions)
        if solution.closure <= solution.tolerance and not duplicate:
            if case == 2 and abs(determinant) < policy.conditioning_warning:
                solution.warnings.append("Ill-conditioned bearings: small input changes can "
                                         "greatly change lengths "
                                         f"(condition number {condition_number:.6g}).")
            solutions.append(solution)
    if not solutions:
        raise TraverseError("Closure requires a zero or negative solved length, or inconsistent "
                            "inputs. Check for a reversed bearing (180°), the bearing quadrant, "
                            "side order, and entered lengths.")
    solutions.sort(key=lambda solution: solution.area is None)
    return SolveResult(case, solutions, explain_result(sides, case, target, solutions,
                                                      condition_number), condition_number)


def example(case):
    sides = [TraverseSide(80, 0), TraverseSide(100, 90),
             TraverseSide(80, 180), TraverseSide(100, 270)]
    omissions = (([3], [3]), ([0], [2]), ([0, 1], []), ([], [0, 1]))
    lengths, bearings = omissions[case]
    return [TraverseSide(None if i in lengths else s.length,
                         None if i in bearings else s.azimuth) for i, s in enumerate(sides)]


# ---------------------- Input model and report files ----------------------

@dataclass
class DraftSide:
    length: str = ""
    bearing: str = ""
    length_unknown: bool = False
    bearing_unknown: bool = False


@dataclass
class InputModel:
    rows: list[DraftSide]
    initial_e: str = "0"
    initial_n: str = "0"
    unit: str = "m"

    def sides(self):
        sides = []
        for i, row in enumerate(self.rows, 1):
            try:
                length = None if row.length_unknown else number(row.length, "Length")
                bearing = None if row.bearing_unknown else parse_bearing(row.bearing)
                sides.append(TraverseSide(length, bearing))
            except TraverseError as error:
                raise TraverseError(f"Side {i}: {error}") from None
        return sides

    def case(self):
        return case_from_unknowns(
            [i for i, row in enumerate(self.rows) if row.length_unknown],
            [i for i, row in enumerate(self.rows) if row.bearing_unknown])

    def resize(self, count):
        if type(count) is not int or not 3 <= count <= 20:
            raise TraverseError("Number of sides must be an integer from 3 to 20.")
        self.rows = (self.rows + [DraftSide() for _ in range(count)])[:count]

    @classmethod
    def from_example(cls, case):
        return cls([DraftSide("" if s.length is None else str(s.length),
                              "" if s.azimuth is None else str(s.azimuth),
                              s.length is None, s.azimuth is None) for s in example(case)])

    def save(self, path):
        data = {"version": 1, "rows": [asdict(row) for row in self.rows],
                "initial_e": self.initial_e, "initial_n": self.initial_n, "unit": self.unit}
        Path(path).write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")

    @classmethod
    def load(cls, path):
        data = json.loads(Path(path).read_text(encoding="utf-8"))
        valid = isinstance(data, dict) and type(data.get("version")) is int
        if not valid or data["version"] != 1:
            raise TraverseError("This is not a supported traverse input file (version 1).")
        rows = data.get("rows")
        if not isinstance(rows, list) or not 3 <= len(rows) <= 20:
            raise TraverseError("Input file must contain 3–20 rows.")
        fields = {"length": str, "bearing": str,
                  "length_unknown": bool, "bearing_unknown": bool}
        for row in rows:
            if not isinstance(row, dict) or set(row) != set(fields):
                raise TraverseError("Invalid side fields in input file.")
            if any(type(row[key]) is not kind for key, kind in fields.items()):
                raise TraverseError("Side values must be text and unknown flags must be booleans.")
        if data.get("unit") not in ("m", "ft") or any(
                not isinstance(data.get(key), str) for key in ("initial_e", "initial_n")):
            raise TraverseError("Invalid units or starting coordinates in input file.")
        return cls([DraftSide(**row) for row in rows], data["initial_e"],
                   data["initial_n"], data["unit"])


def calculation_rows(solution):
    for i, side in enumerate(solution.sides):
        start_e, start_n = solution.coordinates[i]
        end_e, end_n = solution.coordinates[i + 1]
        east, north = solution.components[i]
        yield (i + 1, side.length, quadrantal(side.azimuth), side.azimuth,
               north, east, start_n, start_e, end_n, end_e)


TABLE_HEADERS = ("Side", "Length", "Bearing", "Azimuth", "ΔN", "ΔE",
                 "Start N", "Start E", "End N", "End E")


def solution_summary(model, result, selected):
    solution = result.solutions[selected]
    unit = model.unit
    lines = [CASES[result.case], f"Solution {selected + 1}: {solution.geometry_label}",
             "", "SOLVED VALUES"]
    for i, (row, side) in enumerate(zip(model.rows, solution.sides), 1):
        if row.length_unknown:
            lines.append(f"Side {i} length: {side.length:.4f} {unit}")
        if row.bearing_unknown:
            lines.append(f"Side {i} bearing: {display_azimuth(side.azimuth)}° / "
                         f"{quadrantal(side.azimuth)}")
    area = f"{solution.area:,.4f} {unit}²" if solution.area is not None else "Unavailable"
    lines += ["", f"Perimeter: {solution.perimeter:,.4f} {unit}", f"Area: {area}"]
    lines += ["", *[f"WARNING: {warning}" for warning in solution.warnings]]
    if len(result.solutions) > 1:
        lines += ["Multiple equation solutions exist. Review geometry and engineering constraints."]
    lines += ["", "Numerical residuals are in Calculation Details.",
              "Closure is imposed; it does not measure survey accuracy."]
    return "\n".join(lines)


def export_report(path, model, result, selected=0):
    csv_output = str(path).lower().endswith(".csv")
    with open(path, "w", encoding="utf-8-sig" if csv_output else "utf-8", newline="") as stream:
        writer = csv.writer(stream) if csv_output else None

        def section(title, headers, rows):
            rows = [[str(value) for value in row] for row in rows]
            if writer:
                writer.writerow([title])
                writer.writerow(headers)
                writer.writerows(rows)
                writer.writerow([])
            else:
                stream.write(title + "\n" + "=" * len(title) + "\n")
                if headers:
                    widths = [max(len(row[i]) for row in [list(headers)] + rows)
                              for i in range(len(headers))]
                    for row in [list(headers)] + rows:
                        stream.write("  ".join(v.ljust(w) for v, w in zip(row, widths)) + "\n")
                else:
                    stream.write("\n".join(row[0] for row in rows) + "\n")
                stream.write("\n")

        section("ENGINEERING TRAVERSE SOLVER", (), [[CASES[result.case]],
                [f"Units: {model.unit}; selected solution: {selected + 1}"]])
        section("INPUT DATA", ("Side", "Length", "Bearing"),
                [[i, "Unknown" if r.length_unknown else r.length,
                  "Unknown" if r.bearing_unknown else r.bearing]
                 for i, r in enumerate(model.rows, 1)])
        for index, solution in enumerate(result.solutions):
            section(f"SOLUTION {index + 1}", (),
                    [[line] for line in solution_summary(model, result, index).splitlines()])
            section("TRAVERSE CALCULATIONS", TABLE_HEADERS,
                    [[f"{v:.8f}" if isinstance(v, float) else v for v in row]
                     for row in calculation_rows(solution)])
            section("COORDINATES", ("Point", "Northing", "Easting"),
                    [[f"P{i + 1}" if i < len(solution.sides) else "Return to P1",
                      f"{n:.8f}", f"{e:.8f}"]
                     for i, (e, n) in enumerate(solution.coordinates)])
        section("CALCULATION DETAILS", (), [[line] for line in result.details.splitlines()])


# -------------------------- Desktop views -------------------------------

try:
    import tkinter as tk
    from tkinter import filedialog, messagebox, ttk
except ImportError:
    tk = None


def screen_geometry(width, height):
    """Leave room for desktop panels and window decorations, including on 1366×768."""
    return min(1200, max(320, width - 60)), min(820, max(240, height - 110))


class PlotView:
    """One plotting backend. An unavailable/broken Matplotlib cannot stop the solver."""
    def __init__(self, parent):
        self.figure = self.canvas = self.axes = None
        frame = ttk.Frame(parent)
        frame.pack(fill="both", expand=True)
        self.message = ttk.Label(frame, wraplength=350, padding=8)
        self.message.pack(side="bottom", fill="x")
        try:
            from matplotlib.figure import Figure
            from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
            # Avoid layout= so releases before Matplotlib 3.5 remain compatible.
            self.figure = Figure(figsize=(4, 2.5), dpi=100, tight_layout=True)
            self.axes = self.figure.add_subplot(111)
            self.canvas = FigureCanvasTkAgg(self.figure, master=frame)
            toolbar = NavigationToolbar2Tk(self.canvas, frame, pack_toolbar=False)
            toolbar.pack(side="bottom", fill="x")
            self.canvas.get_tk_widget().pack(fill="both", expand=True)
        except Exception as error:
            for widget in frame.winfo_children():
                if widget is not self.message:
                    widget.destroy()
            self.figure = self.canvas = self.axes = None
            self.message.configure(text=f"Plot unavailable ({type(error).__name__}). "
                                   "Install or update: python -m pip install -U matplotlib. "
                                   "Numerical solving and export remain available.")
        self.render(None, None, 0)

    def render(self, model, result, selected):
        if self.axes is None:
            return
        axes = self.axes
        axes.clear()
        if result:
            solution = result.solutions[selected]
            # Plot local coordinates; global offsets must not obscure small residuals.
            for i, (a, b) in enumerate(zip(solution.local_coordinates,
                                         solution.local_coordinates[1:])):
                row = model.rows[i]
                color = "#e11d48" if row.length_unknown or row.bearing_unknown else "#2563eb"
                axes.plot([a[0], b[0]], [a[1], b[1]], "o-", color=color, linewidth=2)
                axes.annotate("", xy=b, xytext=a, arrowprops={"arrowstyle": "->", "color": color})
                axes.annotate(f"P{i + 1}", a, xytext=(5, 5), textcoords="offset points")
            axes.plot(0, 0, "o", color="#d18b09", markersize=8)
            axes.set_xlabel(f"Local Easting ({model.unit})")
            axes.set_ylabel(f"Local Northing ({model.unit})")
            axes.set_title(f"Solution {selected + 1}: {solution.geometry_label}", fontsize=10)
            self.message.configure(text="Blue: known · Red: solved · Orange: P1\n"
                                   f"Origin E={model.initial_e}, N={model.initial_n}")
        else:
            axes.set_title("Solve a traverse to plot its geometry", fontsize=10)
            self.message.configure(text="Plotting uses solved geometry only.")
        axes.set_aspect("equal", adjustable="datalim")
        axes.grid(True, alpha=0.3)
        self.canvas.draw_idle()


class InputTable:
    def __init__(self, parent, model, on_change):
        self.frame = ttk.Frame(parent)
        self.model, self.on_change, self.variables = model, on_change, []
        self.canvas = tk.Canvas(self.frame, height=180, highlightthickness=0)
        scrollbar = ttk.Scrollbar(self.frame, command=self.canvas.yview)
        self.canvas.configure(yscrollcommand=scrollbar.set)
        scrollbar.pack(side="right", fill="y")
        self.canvas.pack(side="left", fill="both", expand=True)
        self.body = ttk.Frame(self.canvas, padding=5)
        window = self.canvas.create_window((0, 0), window=self.body, anchor="nw")
        self.body.bind("<Configure>", lambda _: self.canvas.configure(
            scrollregion=self.canvas.bbox("all")))
        self.canvas.bind("<Configure>", lambda e: self.canvas.itemconfigure(window, width=e.width))
        for event in ("<MouseWheel>", "<Button-4>", "<Button-5>"):
            parent.winfo_toplevel().bind(event, self.scroll, add="+")
        self.rebuild()

    def rebuild(self):
        for widget in self.body.winfo_children():
            widget.destroy()
        self.variables = []
        headers = ("Side", "Length", "Unknown length", "Bearing / Azimuth", "Unknown bearing")
        for column, text in enumerate(headers):
            ttk.Label(self.body, text=text).grid(row=0, column=column, padx=8, sticky="w")
        self.body.columnconfigure(3, weight=1)
        for i, row in enumerate(self.model.rows):
            values = [tk.StringVar(value=row.length), tk.BooleanVar(value=row.length_unknown),
                      tk.StringVar(value=row.bearing), tk.BooleanVar(value=row.bearing_unknown)]
            ttk.Label(self.body, text=str(i + 1)).grid(row=i + 1, column=0, padx=8)
            widgets = []
            for column, variable in enumerate(values, 1):
                command = lambda i=i: self.commit(i)
                if column in (2, 4):
                    widget = ttk.Checkbutton(self.body, variable=variable, command=command)
                else:
                    widget = ttk.Entry(self.body, textvariable=variable, width=20)
                    widget.bind("<KeyRelease>", lambda _, i=i: self.commit(i))
                    widget.bind("<FocusOut>", lambda _, i=i: self.commit(i))
                widget.grid(row=i + 1, column=column, sticky="ew", padx=8, pady=3)
                widgets.append(widget)
            self.variables.append((values, widgets))
            self.set_states(i)

    def set_states(self, index):
        values, widgets = self.variables[index]
        widgets[0].configure(state="disabled" if values[1].get() else "normal")
        widgets[2].configure(state="disabled" if values[3].get() else "normal")

    def commit(self, index):
        if index >= len(self.model.rows):
            return
        values, _ = self.variables[index]
        row = DraftSide(values[0].get(), values[2].get(), values[1].get(), values[3].get())
        if row != self.model.rows[index]:
            self.model.rows[index] = row
            self.set_states(index)
            self.on_change(data_changed=True)

    def scroll(self, event):
        x, y = event.x_root - self.canvas.winfo_rootx(), event.y_root - self.canvas.winfo_rooty()
        if 0 <= x < self.canvas.winfo_width() and 0 <= y < self.canvas.winfo_height():
            steps = -1 if getattr(event, "num", None) == 4 else 1
            if getattr(event, "delta", 0):
                steps = -max(1, abs(event.delta) // 120) if event.delta > 0 else max(
                    1, abs(event.delta) // 120)
            self.canvas.yview_scroll(steps, "units")
            return "break"


def text_panel(parent):
    frame = ttk.Frame(parent)
    text = tk.Text(frame, wrap="word", height=6, width=40, padx=10, pady=8, state="disabled")
    scrollbar = ttk.Scrollbar(frame, command=text.yview)
    text.configure(yscrollcommand=scrollbar.set)
    scrollbar.pack(side="right", fill="y")
    text.pack(fill="both", expand=True)
    return frame, text


def set_text(widget, content):
    widget.configure(state="normal")
    widget.delete("1.0", "end")
    widget.insert("1.0", content)
    widget.configure(state="disabled")


class TraverseApp:
    def __init__(self, root):
        self.root, self.model = root, InputModel([DraftSide() for _ in range(4)])
        self.result, self.selected = None, 0
        self.plot_window = self.popup = None
        width, height = screen_geometry(root.winfo_screenwidth(), root.winfo_screenheight())
        root.geometry(f"{width}x{height}+20+20")
        root.minsize(min(850, width), min(500, height))
        root.title("Engineering 1035 · Traverse Solver")
        style = ttk.Style(root)
        style.theme_use("clam")
        root.configure(bg="#edf2f7")
        style.configure("TFrame", background="#edf2f7")
        style.configure("TLabel", background="#edf2f7", foreground="#18314b")
        style.configure("Header.TLabel", background="#18314b", foreground="white")
        style.configure("TButton", padding=(8, 5))
        style.configure("Treeview", rowheight=25)
        ttk.Label(root, text="ENGINEERING TRAVERSE SOLVER", font=("Helvetica", 17, "bold"),
                  padding=10, style="Header.TLabel").pack(fill="x")
        self.status = tk.StringVar(value="Ready · Load an example or enter your traverse.")
        ttk.Label(root, textvariable=self.status, padding=(10, 5)).pack(side="bottom", fill="x")
        main = ttk.Frame(root, padding=(10, 0, 10, 5))
        main.pack(fill="both", expand=True)
        settings = ttk.Frame(main)
        settings.pack(fill="x", pady=3)
        self.count = tk.StringVar(value="4")
        ttk.Label(settings, text="Sides").pack(side="left")
        count = ttk.Spinbox(settings, from_=3, to=20, width=4,
                            textvariable=self.count, command=self.resize)
        count.pack(side="left", padx=(4, 12))
        count.bind("<Return>", self.resize)
        count.bind("<FocusOut>", self.resize)
        self.origin_vars = []
        for name, text in (("initial_n", "Initial N"), ("initial_e", "Initial E")):
            ttk.Label(settings, text=text).pack(side="left")
            variable = tk.StringVar(value=getattr(self.model, name))
            self.origin_vars.append(variable)
            entry = ttk.Entry(settings, textvariable=variable, width=12)
            entry.pack(side="left", padx=(4, 12))
            entry.bind("<KeyRelease>", lambda _, n=name, v=variable: self.change_origin(n, v.get()))
            entry.bind("<FocusOut>", lambda _, n=name, v=variable: self.change_origin(n, v.get()))
        self.unit = tk.StringVar(value="m")
        unit = ttk.Combobox(settings, textvariable=self.unit, values=("m", "ft"),
                            state="readonly", width=4)
        unit.pack(side="left")
        unit.bind("<<ComboboxSelected>>", self.change_unit)
        ttk.Button(settings, text="Save Inputs", command=lambda: self.input_file(False)).pack(
            side="right", padx=4)
        ttk.Button(settings, text="Load Inputs", command=lambda: self.input_file(True)).pack(
            side="right")
        ttk.Label(main, text="Bearings: 125.5 · 125°30'00\" · N 30° E. "
                  "0°=North, 90°=East. Units label values; they do not convert them.").pack(
                      anchor="w", pady=3)
        actions = ttk.Frame(main)
        actions.pack(fill="x", pady=3)
        ttk.Button(actions, text="Solve Traverse", command=self.solve).pack(side="left")
        self.example_case = ttk.Combobox(actions,
                                         values=("Case I", "Case II", "Case III", "Case IV"),
                                         state="readonly", width=9)
        self.example_case.current(0)
        self.example_case.pack(side="left", padx=(12, 4))
        ttk.Button(actions, text="Load Example", command=self.load_example).pack(side="left")
        ttk.Button(actions, text="Clear", command=self.clear).pack(side="left", padx=4)
        self.export_button = ttk.Button(actions, text="Export Results", command=self.export)
        self.export_button.pack(side="right")
        self.plot_button = ttk.Button(actions, text="Enlarge Plot", command=self.enlarge_plot)
        self.plot_button.pack(side="right", padx=4)
        self.case_label = ttk.Label(main, wraplength=1100)
        self.case_label.pack(fill="x", pady=3)
        self.panes = ttk.Panedwindow(main, orient="vertical")
        self.panes.pack(fill="both", expand=True)
        self.table = InputTable(self.panes, self.model, self.refresh)
        self.panes.add(self.table.frame, weight=1)
        results = ttk.Frame(self.panes)
        self.panes.add(results, weight=2)
        self.choice = ttk.Combobox(results, state="disabled")
        self.choice.pack(fill="x", pady=(4, 2))
        self.choice.bind("<<ComboboxSelected>>", self.select_solution)
        self.result_panes = ttk.Panedwindow(results, orient="horizontal")
        self.result_panes.pack(fill="both", expand=True)
        tabs = ttk.Notebook(self.result_panes)
        self.result_panes.add(tabs, weight=1)
        frame, self.summary = text_panel(tabs)
        tabs.add(frame, text="Results")
        frame, self.details = text_panel(tabs)
        tabs.add(frame, text="Calculation Details")
        frame = ttk.Frame(tabs)
        tabs.add(frame, text="Side Calculations & Coordinates")
        self.calculations = ttk.Treeview(frame, columns=TABLE_HEADERS, show="headings", height=4)
        for column in TABLE_HEADERS:
            self.calculations.heading(column, text=column)
            self.calculations.column(column, width=120, stretch=False)
        horizontal = ttk.Scrollbar(frame, orient="horizontal", command=self.calculations.xview)
        vertical = ttk.Scrollbar(frame, command=self.calculations.yview)
        self.calculations.configure(xscrollcommand=horizontal.set, yscrollcommand=vertical.set)
        horizontal.pack(side="bottom", fill="x")
        vertical.pack(side="right", fill="y")
        self.calculations.pack(fill="both", expand=True)
        plot_frame = ttk.Frame(self.result_panes)
        self.result_panes.add(plot_frame, weight=1)
        self.plot = PlotView(plot_frame)
        self.panes.bind("<Map>", lambda _: root.after_idle(self.balance_panes))
        self.refresh()

    def balance_panes(self):
        self.panes.sashpos(0, int(self.panes.winfo_height() * 0.33))
        self.result_panes.sashpos(0, int(self.result_panes.winfo_width() * 0.52))

    def refresh(self, data_changed=False):
        if data_changed:
            self.result = None
            self.status.set("Inputs changed · Solve to update results.")
        try:
            self.case_label.configure(text="Detected: " + CASES[self.model.case()])
        except TraverseError:
            self.case_label.configure(
                text="Mark two unknown quantities to select Case I–IV automatically.")
        self.calculations.delete(*self.calculations.get_children())
        self.export_button.configure(state="normal" if self.result else "disabled")
        can_plot = self.result is not None and self.plot.figure is not None
        self.plot_button.configure(state="normal" if can_plot else "disabled")
        if self.result:
            labels = [f"Solution {i + 1} — {s.geometry_label}"
                      for i, s in enumerate(self.result.solutions)]
            self.choice.configure(state="readonly", values=labels)
            self.selected = min(self.selected, len(labels) - 1)
            self.choice.current(self.selected)
            set_text(self.summary, solution_summary(self.model, self.result, self.selected))
            set_text(self.details, self.result.details)
            for row in calculation_rows(self.result.solutions[self.selected]):
                values = [f"{v:.4f}" if isinstance(v, float) else v for v in row]
                values[3] = display_azimuth(row[3])
                self.calculations.insert("", "end", values=values)
        else:
            self.choice.configure(state="disabled", values=())
            self.choice.set("")
            set_text(self.summary, "Enter known values, mark two Unknown fields, then solve.")
            set_text(self.details, "")
        self.plot.render(self.model, self.result, self.selected)
        if self.popup:
            self.popup.render(self.model, self.result, self.selected)

    def resize(self, *_):
        try:
            if not re.fullmatch(r"[0-9]+", self.count.get()):
                raise TraverseError("Enter an integer side count from 3 to 20.")
            count = int(self.count.get())
            if count != len(self.model.rows):
                self.model.resize(count)
                self.table.rebuild()
                self.refresh(data_changed=True)
            return True
        except TraverseError as error:
            self.count.set(str(len(self.model.rows)))
            messagebox.showerror("Invalid Side Count", str(error), parent=self.root)
            return False

    def change_origin(self, name, value):
        if value != getattr(self.model, name):
            setattr(self.model, name, value)
            self.refresh(data_changed=True)

    def change_unit(self, *_):
        self.model.unit = self.unit.get()
        self.refresh()

    def select_solution(self, *_):
        self.selected = self.choice.current()
        self.refresh()

    def replace_model(self, model):
        self.model = self.table.model = model
        self.count.set(str(len(model.rows)))
        self.origin_vars[0].set(model.initial_n)
        self.origin_vars[1].set(model.initial_e)
        self.unit.set(model.unit)
        self.table.rebuild()
        self.selected = 0
        self.refresh(data_changed=True)

    def load_example(self):
        self.replace_model(InputModel.from_example(self.example_case.current()))

    def clear(self):
        self.replace_model(InputModel([DraftSide() for _ in self.model.rows], unit=self.model.unit))

    def solve(self):
        if not self.resize():
            return
        try:
            result = solve_traverse(self.model.sides(), self.model.initial_e, self.model.initial_n)
        except (TraverseError, ArithmeticError) as error:
            message = str(error) if isinstance(error, TraverseError) else (
                "Values exceed the supported numeric range. Use smaller lengths or "
                "coordinates, or choose units that keep the entered values manageable.")
            messagebox.showerror("Cannot Solve", message, parent=self.root)
            return
        self.result, self.selected = result, 0
        self.refresh()
        self.status.set(f"Solved · {len(result.solutions)} equation solution(s). "
                        "Review geometry warnings.")

    def input_file(self, loading):
        path = (filedialog.askopenfilename if loading else filedialog.asksaveasfilename)(
            parent=self.root, defaultextension=".json", filetypes=(("Traverse inputs", "*.json"),))
        if not path:
            return
        try:
            if loading:
                self.replace_model(InputModel.load(path))
            else:
                self.model.save(path)
            self.status.set("Inputs loaded." if loading else
                            "Inputs saved (including unknown flags).")
        except (OSError, ValueError, TypeError) as error:
            messagebox.showerror("Input File Error", str(error), parent=self.root)

    def enlarge_plot(self):
        if self.result is None or self.plot.figure is None:
            return
        if self.plot_window is None:
            self.plot_window = tk.Toplevel(self.root)
            width, height = screen_geometry(self.root.winfo_screenwidth(),
                                             self.root.winfo_screenheight())
            self.plot_window.geometry(f"{min(850, width)}x{height}+30+20")
            self.plot_window.title("Solved Traverse Plot")
            self.popup = PlotView(self.plot_window)
            self.plot_window.protocol("WM_DELETE_WINDOW", self.close_plot)
        self.popup.render(self.model, self.result, self.selected)
        self.plot_window.lift()

    def close_plot(self):
        self.plot_window.destroy()
        self.plot_window = self.popup = None

    def export(self):
        if self.result is None:
            return
        path = filedialog.asksaveasfilename(parent=self.root, defaultextension=".txt",
                    filetypes=(("Text report", "*.txt"), ("CSV report", "*.csv")))
        if path:
            try:
                export_report(path, self.model, self.result, self.selected)
                self.status.set("Exported all solutions and numerical details.")
            except OSError as error:
                messagebox.showerror("Export Failed", str(error), parent=self.root)


def main():
    if tk is None:
        print("Tkinter is required (Ubuntu/Debian: install python3-tk).", file=sys.stderr)
        return 1
    try:
        root = tk.Tk()
    except tk.TclError as error:
        print(f"A desktop display is required: {error}", file=sys.stderr)
        return 1
    TraverseApp(root)
    root.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
