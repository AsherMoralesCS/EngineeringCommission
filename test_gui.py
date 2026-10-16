"""Display-dependent regression tests: ``python -m unittest -v test_gui``."""
import math
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import traverse_solver as app_module
from traverse_solver import DraftSide, InputModel, TraverseApp


class GuiRegressionTests(unittest.TestCase):
    def setUp(self):
        if app_module.tk is None:
            self.skipTest("Tkinter is unavailable")
        try:
            self.root = app_module.tk.Tk()
        except app_module.tk.TclError as error:
            self.skipTest(f"A desktop display is required: {error}")
        self.width = patch.object(self.root, "winfo_screenwidth", return_value=1366)
        self.height = patch.object(self.root, "winfo_screenheight", return_value=768)
        self.errors = patch.object(app_module.messagebox, "showerror")
        self.width.start()
        self.height.start()
        self.showerror = self.errors.start()

    def tearDown(self):
        if hasattr(self, "root"):
            self.root.update_idletasks()
            self.root.destroy()
            self.width.stop()
            self.height.stop()
            self.errors.stop()
            self.showerror.assert_not_called()

    def make_app(self):
        app = TraverseApp(self.root)
        self.root.update()
        return app

    def load_and_solve(self, app, case=0):
        app.example_case.current(case)
        app.load_example()
        app.solve()
        self.root.update()
        self.assertIsNotNone(app.result)

    def test_laptop_window_footer_and_twenty_row_mouse_scrolling(self):
        app = self.make_app()
        self.assertLessEqual(self.root.winfo_y() + self.root.winfo_height(), 768)
        self.assertLessEqual(self.root.winfo_x() + self.root.winfo_width(), 1366)
        footer = next(widget for widget in self.root.winfo_children()
                      if isinstance(widget, app_module.ttk.Label)
                      and str(widget.cget("textvariable")) == str(app.status))
        app.table.variables[0][0][0].set("123.45")
        app.table.commit(0)
        app.count.set("20")
        app.resize()
        self.root.update()
        self.assertEqual(len(app.model.rows), 20)
        self.assertEqual(app.model.rows[0].length, "123.45")
        self.assertGreater(footer.winfo_height(), 10)
        self.assertLessEqual(footer.winfo_y() + footer.winfo_height(), self.root.winfo_height())
        canvas = app.table.canvas
        canvas.yview_moveto(0)
        event = SimpleNamespace(x_root=canvas.winfo_rootx() + 10,
                                y_root=canvas.winfo_rooty() + 10, num=5, delta=0)
        self.assertEqual(app.table.scroll(event), "break")
        self.assertGreater(canvas.yview()[0], 0)

    def test_units_preserve_results_and_data_edits_invalidate_them(self):
        app = self.make_app()
        self.load_and_solve(app)
        result = app.result
        original_sides = list(result.solutions[0].sides)
        app.unit.set("ft")
        app.change_unit()
        self.assertIs(app.result, result)
        self.assertEqual(result.solutions[0].sides, original_sides)
        self.assertIn("ft", app.summary.get("1.0", "end"))
        app.table.variables[0][0][0].set("81")
        app.table.commit(0)
        self.assertIsNone(app.result)
        self.assertTrue(app.export_button.instate(["disabled"]))

    def assert_plot_construction_failure_is_recoverable(self, error):
        try:
            import matplotlib.figure
        except ImportError:
            self.skipTest("Matplotlib is unavailable")
        with patch("matplotlib.figure.Figure", side_effect=error):
            app = self.make_app()
        self.assertIsNone(app.plot.figure)
        self.assertIsNone(app.plot.canvas)
        self.assertIn("Plot unavailable", app.plot.message.cget("text"))
        self.load_and_solve(app)
        self.assertTrue(app.plot_button.instate(["disabled"]))
        self.assertFalse(app.export_button.instate(["disabled"]))

    def test_older_matplotlib_typeerror_keeps_solver_usable(self):
        self.assert_plot_construction_failure_is_recoverable(TypeError("unsupported option"))

    def test_matplotlib_runtimeerror_keeps_solver_usable(self):
        self.assert_plot_construction_failure_is_recoverable(RuntimeError("broken backend"))

    def test_large_offset_plot_uses_local_geometry_and_flags_crossed_choice(self):
        app = self.make_app()
        if app.plot.axes is None:
            self.skipTest("Matplotlib is unavailable")
        model = InputModel.from_example(3)
        model.initial_e, model.initial_n = "1000000000000", "-1000000000000"
        app.replace_model(model)
        app.solve()
        self.root.update()
        self.assertIsNotNone(app.result)
        choices = app.choice.cget("values")
        self.assertIn("Simple polygon", choices[0])
        self.assertIn("Crossed or degenerate", choices[1])
        labels = [text.get_text() for text in app.plot.axes.texts if text.get_text()]
        self.assertEqual(labels, ["P1", "P2", "P3", "P4"])
        plotted_x = [abs(float(x)) for line in app.plot.axes.lines for x in line.get_xdata()]
        self.assertLessEqual(max(plotted_x), 101)
        app.choice.current(1)
        app.select_solution()
        self.assertEqual(app.selected, 1)
        self.assertIn("Crossed", app.summary.get("1.0", "end"))

    def test_nearly_parallel_warning_is_visible_in_results(self):
        app = self.make_app()
        angle = math.radians(.01)
        east, north = 100 * math.sin(angle), 80 + 100 * math.cos(angle)
        closing_bearing = math.degrees(math.atan2(-east, -north)) % 360
        app.replace_model(InputModel([
            DraftSide("", "0", True, False),
            DraftSide("", ".01", True, False),
            DraftSide(str(math.hypot(east, north)), str(closing_bearing)),
        ]))
        app.solve()
        self.assertIsNotNone(app.result)
        self.assertIn("condition", app.summary.get("1.0", "end").lower())
        self.assertGreater(app.result.condition_number, 10000)

    def test_save_clear_and_load_inputs_through_dialog_actions(self):
        app = self.make_app()
        app.replace_model(InputModel.from_example(3))
        saved_model = app.model
        with tempfile.TemporaryDirectory() as directory:
            path = str(Path(directory) / "inputs.json")
            with patch.object(app_module.filedialog, "asksaveasfilename", return_value=path):
                app.input_file(False)
            app.clear()
            self.assertNotEqual(app.model, saved_model)
            with patch.object(app_module.filedialog, "askopenfilename", return_value=path):
                app.input_file(True)
        self.assertEqual(app.model, saved_model)
        self.assertEqual(app.table.model, saved_model)
        self.assertIsNone(app.result)


if __name__ == "__main__":
    unittest.main()
