Engineering Traverse Solver — source-based flowcharts

Open Engineering_Traverse_Solver_Flowcharts.drawio in diagrams.net to edit.
The PDF contains vector pages and a bookmarked contents panel.
Each SVG is a standalone scalable page. Page identifiers are continuation references.

Source commit: 0cf64fb29a4615db97badf16442ca089ce13d9aa
Source SHA-256: e18d3a23955f429268eb674285ee16d754702207062ae6a13b631dfc1c5613dc
Only traverse_solver.py was used to document behavior; the application was not modified.

Symbols: oval=start/return/raise; rectangle=process; diamond=decision;
parallelogram=I/O; double vertical bars=predefined subprocess;
gray reference panels=non-executable source notes. Every decision branch is labeled.

PAGE INDEX
01. 01 — Program entry and Tkinter event loop
02. 02a — Input-table edits and scrolling
03. 02b — Side count changes and validation
04. 02c — Central view refresh
05. 02d — Origin, unit, and solution-selection callbacks
06. 02e — Example, clear, and model replacement
07. 03a — Parse input rows and decimal numbers
08. 03b — Angle parsing: decimal or DMS
09. 03c — Bearing parsing: quadrantal or whole-circle
10. 03d — Detect the omitted-measurement case
11. 03e — Format quadrantal and azimuth bearings
12. 04a — Solve callback, validation, and solver dispatch
13. 04b — Candidate acceptance, conditioning, and result return
14. 05 — Cases I and III — omitted lengths
15. 06 — Case II — omitted length and another bearing
16. 07 — Case IV — two omitted bearings
17. 08a — Evaluate one candidate
18. 08b — Polygon boundary checks and shoelace area
19. 09a — Input files and result export
20. 09b — JSON input validation
21. 10 — Plot backend initialization and recovery
22. 10a — Solved plot rendering
23. 10b — Enlarged plot window and close callback
24. 11a — Export report sections and all-solution loop
25. 11b — User-facing solution summary
26. 11c — CSV or aligned-text section output
27. 11d — Calculation explanation and solved-value details
