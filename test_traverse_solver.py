"""Run with ``python -m unittest -v test_traverse_solver``."""
import math
from dataclasses import is_dataclass, replace
import csv
import json
from pathlib import Path
import tempfile
import unittest

from traverse_solver import (
    TraverseError,
    TraverseSide,
    InputModel,
    display_azimuth,
    example,
    export_report,
    number,
    parse_bearing,
    quadrantal,
    solve_traverse,
)


def bearing(east, north):
    return math.degrees(math.atan2(east, north)) % 360


def complete_pentagon():
    points = [(0, 0), (3, 5), (9, 6), (12, 2), (5, -3)]
    return [
        TraverseSide(math.hypot(b[0] - a[0], b[1] - a[1]),
                     bearing(b[0] - a[0], b[1] - a[1]))
        for a, b in zip(points, points[1:] + points[:1])
    ]


class SolverTests(unittest.TestCase):
    def assert_closed(self, solution):
        self.assertLessEqual(solution.closure, solution.tolerance)
        self.assertAlmostEqual(solution.local_coordinates[-1][0], 0, places=8)
        self.assertAlmostEqual(solution.local_coordinates[-1][1], 0, places=8)
        self.assertTrue(all(side.length > 0 for side in solution.sides))

    def test_four_examples_and_named_result(self):
        for case in range(4):
            with self.subTest(case=case):
                result = solve_traverse(example(case), 1500, 2200)
                self.assertTrue(is_dataclass(result))
                self.assertEqual(result.case, case)
                self.assertIn("Closure", result.details)
                self.assertTrue(result.solutions)
                for solution in result.solutions:
                    self.assert_closed(solution)
                    self.assertAlmostEqual(solution.perimeter, 360)
                    self.assertAlmostEqual(solution.coordinates[-1][0], 1500)
                    self.assertAlmostEqual(solution.coordinates[-1][1], 2200)
                self.assertAlmostEqual(result.solutions[0].area, 8000)

    def test_nonrectangular_polygon_all_cases_and_unknown_positions(self):
        complete = complete_pentagon()
        for case in range(4):
            for i in range(len(complete)):
                j = (i + 1) % len(complete)
                sides = list(complete)
                if case == 0:
                    sides[i] = TraverseSide(None, None)
                elif case == 1:
                    sides[i] = replace(sides[i], length=None)
                    sides[j] = replace(sides[j], azimuth=None)
                elif case == 2:
                    sides[i] = replace(sides[i], length=None)
                    sides[j] = replace(sides[j], length=None)
                else:
                    sides[i] = replace(sides[i], azimuth=None)
                    sides[j] = replace(sides[j], azimuth=None)
                with self.subTest(case=case, omitted=i):
                    result = solve_traverse(sides)
                    # Compare vectors to avoid a special case at North (0 == 360).
                    expected = [
                        (s.length * math.sin(math.radians(s.azimuth)),
                         s.length * math.cos(math.radians(s.azimuth)))
                        for s in complete
                    ]
                    self.assertTrue(any(
                        all(math.dist(actual, wanted) < 1e-7
                            for actual, wanted in zip(solution.components, expected))
                        for solution in result.solutions
                    ))
                    for solution in result.solutions:
                        self.assert_closed(solution)

    def test_nearly_parallel_lengths_report_conditioning(self):
        angle = math.radians(0.01)
        east = 100 * math.sin(angle)
        north = 80 + 100 * math.cos(angle)
        sides = [TraverseSide(None, 0), TraverseSide(None, 0.01),
                 TraverseSide(math.hypot(east, north), bearing(-east, -north))]
        result = solve_traverse(sides)
        self.assertAlmostEqual(result.condition_number, 11459.1559, places=2)
        solution = result.solutions[0]
        self.assertAlmostEqual(solution.sides[0].length, 80, places=6)
        self.assertAlmostEqual(solution.sides[1].length, 100, places=6)
        self.assertTrue(solution.warnings)
        self.assertIn("condition", " ".join(solution.warnings).lower())

    def test_negative_solution_has_actionable_diagnostic(self):
        sides = [TraverseSide(None, 0), TraverseSide(None, 90),
                 TraverseSide(5, 45)]
        with self.assertRaises(TraverseError) as caught:
            solve_traverse(sides)
        diagnostic = str(caught.exception).lower()
        self.assertIn("negative", diagnostic)
        self.assertIn("revers", diagnostic)
        self.assertIn("180", diagnostic)
        self.assertIn("quadrant", diagnostic)

    def test_case_ii_two_positive_solutions_and_tangency(self):
        result = solve_traverse([
            TraverseSide(None, 0), TraverseSide(2, None), TraverseSide(5, 180)
        ])
        self.assertEqual(len(result.solutions), 2)
        self.assertEqual({round(s.sides[0].length) for s in result.solutions}, {3, 7})
        tangent = solve_traverse([
            TraverseSide(None, 0), TraverseSide(5, None),
            TraverseSide(math.sqrt(34), bearing(-5, -3)),
        ])
        self.assertEqual(len(tangent.solutions), 1)
        self.assertAlmostEqual(tangent.solutions[0].sides[0].length, 3, places=6)

    def test_rotated_tangencies_do_not_split_into_false_solutions(self):
        # A ray/circle tangent rotated off North must still have one root.
        angle = math.radians(1)
        east = 3 * math.sin(angle) - 5 * math.cos(angle)
        north = 3 * math.cos(angle) + 5 * math.sin(angle)
        case_ii = solve_traverse([
            TraverseSide(None, 1), TraverseSide(5, None),
            TraverseSide(math.hypot(east, north), bearing(-east, -north)),
        ])
        self.assertEqual(len(case_ii.solutions), 1)
        self.assertAlmostEqual(case_ii.solutions[0].sides[0].length, 3)
        # Two unit vectors adding to a length-two vector are coincident.
        case_iv = solve_traverse([
            TraverseSide(1, None), TraverseSide(1, None), TraverseSide(2, 120),
        ])
        self.assertEqual(len(case_iv.solutions), 1)
        self.assertAlmostEqual(case_iv.solutions[0].sides[0].azimuth, 300)

    def test_case_iv_crossed_solution_is_flagged_and_sorted_last(self):
        result = solve_traverse(example(3))
        self.assertEqual(len(result.solutions), 2)
        first, second = result.solutions
        self.assertEqual(first.geometry_label, "Simple polygon")
        self.assertAlmostEqual(first.area, 8000)
        self.assertEqual(second.geometry_label, "Crossed or degenerate")
        self.assertIsNone(second.area)
        self.assertTrue(second.warnings)
        tangent = solve_traverse([
            TraverseSide(3, None), TraverseSide(2, None), TraverseSide(5, 180)
        ])
        self.assertEqual(len(tangent.solutions), 1)

    def test_impossible_indeterminate_and_singular_inputs(self):
        datasets = {
            "singular": [TraverseSide(None, 0), TraverseSide(None, 180),
                         TraverseSide(3, 90)],
            "case ii no intersection": [TraverseSide(None, 0), TraverseSide(2, None),
                                         TraverseSide(5, 270)],
            "case iv impossible": [TraverseSide(1, None), TraverseSide(1, None),
                                   TraverseSide(5, 90)],
            "case iv indeterminate": [TraverseSide(1, None), TraverseSide(1, None),
                                      TraverseSide(3, 0), TraverseSide(3, 180)],
            "case i indeterminate": [TraverseSide(None, None), TraverseSide(3, 0),
                                     TraverseSide(3, 180)],
            "bad length": [TraverseSide(None, None), TraverseSide(-2, 90),
                           TraverseSide(3, 180)],
            "nonfinite length": [TraverseSide(None, None), TraverseSide(math.nan, 90),
                                 TraverseSide(3, 180)],
            "too few sides": example(0)[:2],
            "no unknowns": [TraverseSide(1, 0)] * 4,
            "too many unknowns": [TraverseSide(None, None)] * 4,
        }
        for label, sides in datasets.items():
            with self.subTest(label=label), self.assertRaises(TraverseError):
                solve_traverse(sides)

    def test_area_and_local_geometry_survive_coordinate_offsets(self):
        result = solve_traverse(example(0), 1e12, -1e12)
        solution = result.solutions[0]
        self.assertAlmostEqual(solution.area, 8000)
        self.assert_closed(solution)
        self.assertLess(math.dist(solution.local_coordinates[0],
                                  solution.local_coordinates[-1]), solution.tolerance)


class ParsingTests(unittest.TestCase):
    def test_ascii_numbers_and_finite_values(self):
        for text, expected in (("1000", 1000), (" -12.5 ", -12.5),
                               (".25", .25), ("+1.2e3", 1200)):
            with self.subTest(text=text):
                self.assertEqual(number(text), expected)
        for text in ("1_000", "１２", "١٢", "nan", "inf", "1e999", "", "1,000"):
            with self.subTest(text=text), self.assertRaises(TraverseError):
                number(text)

    def test_bearing_formats(self):
        examples = (("N 30° E", 30), ("S 45 E", 135), ("S 20 W", 200),
                    ("N 70 W", 290), ('125°30\'00"', 125.5),
                    ("125°30′00″", 125.5), ("125:30:00", 125.5),
                    ("125 30 00", 125.5), ("360°", 0), ("N 0 W", 0))
        for text, expected in examples:
            with self.subTest(text=text):
                self.assertAlmostEqual(parse_bearing(text), expected)
        for text in ("-1", "361", "nan", "inf", "N 91 E", "N 30 N",
                     '20°60\'00"', "20:10:60", "", "1_0", "１２", "N ١٠ E"):
            with self.subTest(text=text), self.assertRaises(TraverseError):
                parse_bearing(text)

    def test_cardinal_rounding_happens_before_quadrant_selection(self):
        expected = {
            89.9999999: 'N 90°00\'00.00" E',
            179.9999999: 'S 0°00\'00.00" E',
            269.9999999: 'N 90°00\'00.00" W',
            359.9999999: 'N 0°00\'00.00" E',
        }
        for angle, text in expected.items():
            with self.subTest(angle=angle):
                self.assertEqual(quadrantal(angle), text)
        self.assertEqual(display_azimuth(359.999999999), "0.0000")
        for angle in range(360):
            with self.subTest(angle=angle):
                self.assertAlmostEqual(parse_bearing(quadrantal(angle)), angle)


class InputModelTests(unittest.TestCase):
    def test_unknown_masks_detect_case_without_dummy_numeric_values(self):
        for case in range(4):
            with self.subTest(case=case):
                model = InputModel.from_example(case)
                self.assertEqual(model.case(), case)
                self.assertEqual(model.sides(), example(case))
                model.rows[2].length = "unfinished input"
                self.assertEqual(model.case(), case)
                with self.assertRaises(TraverseError):
                    model.sides()

    def test_resize_preserves_entered_rows(self):
        model = InputModel.from_example(0)
        model.rows[0].length = "123.45"
        model.resize(7)
        self.assertEqual(len(model.rows), 7)
        self.assertEqual(model.rows[0].length, "123.45")
        model.resize(3)
        self.assertEqual(len(model.rows), 3)
        self.assertEqual(model.rows[0].length, "123.45")

    def test_input_file_round_trip_preserves_draft_text_and_units(self):
        model = InputModel.from_example(3)
        model.initial_e, model.initial_n, model.unit = "1000.5", "-300", "ft"
        model.rows[0].bearing = "unfinished bearing"
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "traverse.json"
            model.save(path)
            restored = InputModel.load(path)
        self.assertEqual(restored, model)
        self.assertEqual(restored.rows[0].bearing, "unfinished bearing")
        # A label change keeps physical values in their original numerical units.
        original = restored.sides()
        restored.unit = "m"
        self.assertEqual(restored.sides(), original)

    def test_input_loader_rejects_invalid_shape(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "traverse.json"
            InputModel.from_example(0).save(path)
            valid = json.loads(path.read_text(encoding="utf-8"))
            mutations = (
                lambda data: data.update(version=99),
                lambda data: data.update(rows=data["rows"][:2]),
                lambda data: data.update(unit="yards"),
                lambda data: data["rows"][0].update(length_unknown="false"),
                lambda data: data["rows"][0].update(length=80),
            )
            for mutate in mutations:
                data = json.loads(json.dumps(valid))
                mutate(data)
                path.write_text(json.dumps(data), encoding="utf-8")
                with self.subTest(data=data), self.assertRaises(TraverseError):
                    InputModel.load(path)


class ExportTests(unittest.TestCase):
    def test_text_report_is_formatted_and_csv_is_parseable(self):
        model = InputModel.from_example(0)
        model.unit = "ft"
        model.initial_e, model.initial_n = "1500", "2200"
        result = solve_traverse(model.sides(), 1500, 2200)
        with tempfile.TemporaryDirectory() as directory:
            text_path = Path(directory) / "report.txt"
            csv_path = Path(directory) / "report.csv"
            export_report(text_path, model, result)
            export_report(csv_path, model, result)
            report = text_path.read_text(encoding="utf-8-sig")
            with csv_path.open(encoding="utf-8-sig", newline="") as file:
                rows = list(csv.reader(file))
        self.assertNotIn("\t", report)
        for expected in ("Perimeter", "Area", "Side", "ft"):
            self.assertIn(expected, report)
            self.assertTrue(any(expected in cell for row in rows for cell in row))
        self.assertTrue(any(len(row) >= 5 for row in rows))


if __name__ == "__main__":
    unittest.main()
