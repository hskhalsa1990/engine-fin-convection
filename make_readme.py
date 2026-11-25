"""Regenerates README.md from the figures and tables in this folder:  python make_readme.py"""
import glob, os
rd = lambda f: open(f).read()
figs = sorted(glob.glob("fig*.png")); g2 = [f for f in figs if "3d" not in f]; g3 = [f for f in figs if "3d" in f]
img = lambda l: "\n".join(f"![{f[:-4]}]({f})" for f in l)
tabs = "\n".join(f"### {os.path.basename(f)[6:-3].replace('_', ' ').title()}\n{rd(f)}" for f in sorted(glob.glob("table_*.md")))
open("README.md", "w").write(f"""# Engine Fin Convection (4-Stroke Engine)

Flat repository (no sub-folders): heat convection from 4-stroke engine cooling fins. Fin theory, finite-difference solvers in four languages, channel flow between fins, pressure drop, fin shapes, spacing and geometry optimisation, annular fins, Monte-Carlo uncertainty, transient warm-up (with animation), ML, graphs, comparison tables and 3D diagrams.

## Files
| File | Purpose |
|---|---|
| `Engine_Fin_Convection_Fluid_Dynamics.ipynb` | Jupyter notebook (Python, C, C++, Rust, Go, Bash, HTML, ML, DL, CV, upgrade sections) |
| `fin_core.py` | Physics: straight, triangular, pin and annular fins; finite difference; Graetz channel flow; Elenbaas natural convection; channel pressure drop |
| `engine_fin.py` | Basic figures, comparison tables, ML, 3D diagrams |
| `engine_fin_advanced.py` | Channel flow, natural-convection spacing, annular fin, Monte-Carlo, transient |
| `engine_fin_optimization.py` | Fin shapes, thermal-hydraulic trade-off, differential-evolution optimisation, animated 2D transient |
| `streamlit_app.py` | Interactive dashboard (`streamlit run streamlit_app.py`) |
| `test_engine_fin.py` | Tests (`pytest -q` or `python test_engine_fin.py`) |
| `fin.c`, `fin.cpp`, `fin.rs`, `fin.go` | Same 1D fin solver in four languages |
| `run_all.sh`, `make_readme.py` | Rebuild everything / regenerate this README |
| `table_*`, `fig*.png`, `anim_fin_transient.gif`, `fin_data_clean.csv` | Generated results |

## Quick start
```bash
pip install -r requirements.txt
./run_all.sh
python make_readme.py
```

## Models and assumptions
- Straight fin: `θ'' − m²θ = 0`, `Q = √(hPkA)·θb·tanh(m·Lc)`; h from flat-plate Nu–Re–Pr correlations. Triangular, pin and annular (Bessel I/K) fins use their standard closed-form efficiencies.
- Channel flow between fins: laminar parabolic flow, isothermal walls, implicit marching energy equation (Graetz). Fully developed Nu converges to 7.54.
- Natural-convection spacing: Elenbaas correlation plus radiation reduced by the fin-to-fin view factor (vertical plate fins). Fin surface assumed at 85 % of the base excess temperature.
- Forced-flow array: laminar (Nu = 7.54 + entrance term, f = 96/Re) or turbulent (Dittus-Boelter, Blasius) in the gap; about 1.5 velocity heads of inlet/outlet loss; bypass flow around the array is ignored, so heat rates are optimistic.
- Optimisation: differential evolution maximises Q for limits on pressure drop and mass, within manufacturing bounds (S 3–20 mm, t 1.5–6 mm, height 10–60 mm). Several optima sit on these bounds or on the mass limit, so the Pareto front flattens at high pressure-drop limits.
- Transient: lumped cylinder head (3 kJ/K, 1 kW, 8 annular fins); animation uses 2D explicit conduction with a hot spot at the base.
- All dimensions, materials and heat inputs are illustrative, not measured engine data. The ML data is synthetic, generated from the same formulas.

## Animation
![fin warm-up](anim_fin_transient.gif)

## Graphs and charts
{img(g2)}

## 3D diagrams
{img(g3)}

Fig13 and fig27 use illustrative h or m values so the temperature gradient is visible.

## Comparison tables
{tabs}
Rust and Go rows in the language table show "not installed" if those compilers were missing when generated. Run times include process startup, except the in-process Python row.

## License
MIT (H.S. Khalsa)
""")
