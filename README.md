# Engine Fin Convection (4-Stroke Engine)

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
![fig01_h_Q_eta_curves](fig01_h_Q_eta_curves.png)
![fig02_temperature_profile](fig02_temperature_profile.png)
![fig03_material_bar](fig03_material_bar.png)
![fig04_efficiency_heatmap](fig04_efficiency_heatmap.png)
![fig05_cooling_pie](fig05_cooling_pie.png)
![fig06_ml_comparison](fig06_ml_comparison.png)
![fig07_correlation](fig07_correlation.png)
![fig08_feature_importance](fig08_feature_importance.png)
![fig09_language_runtime](fig09_language_runtime.png)
![fig14_channel_flow_temperature](fig14_channel_flow_temperature.png)
![fig15_graetz_nusselt](fig15_graetz_nusselt.png)
![fig16_fin_spacing_optimum](fig16_fin_spacing_optimum.png)
![fig17_annular_vs_straight](fig17_annular_vs_straight.png)
![fig18_monte_carlo](fig18_monte_carlo.png)
![fig19_transient_warmup](fig19_transient_warmup.png)
![fig23_fin_shapes_comparison](fig23_fin_shapes_comparison.png)
![fig24_thermal_hydraulic](fig24_thermal_hydraulic.png)
![fig25_pareto_optimisation](fig25_pareto_optimisation.png)
![fig28_transient_2d_snapshots](fig28_transient_2d_snapshots.png)

## 3D diagrams
![fig10_3d_temperature_surface](fig10_3d_temperature_surface.png)
![fig11_3d_Q_surface](fig11_3d_Q_surface.png)
![fig12_3d_bar_materials](fig12_3d_bar_materials.png)
![fig13_3d_engine_fins](fig13_3d_engine_fins.png)
![fig20_3d_channel_temperature](fig20_3d_channel_temperature.png)
![fig21_3d_annular_efficiency](fig21_3d_annular_efficiency.png)
![fig22_3d_spacing_temperature](fig22_3d_spacing_temperature.png)
![fig26_3d_pareto](fig26_3d_pareto.png)
![fig27_3d_fin_shapes](fig27_3d_fin_shapes.png)

Fig13 and fig27 use illustrative h or m values so the temperature gradient is visible.

## Comparison tables
### Annular
| h (W/m2K) | Annular eta | Straight eta | Annular Q per fin (W), theta=160 |
|---|---|---|---|
| 10.0 | 0.984 | 0.988 | 34.7 |
| 30.0 | 0.953 | 0.964 | 101.0 |
| 60.0 | 0.911 | 0.931 | 193.1 |
| 100.0 | 0.862 | 0.89 | 304.2 |

### Channel Nusselt
| U (m/s) | Re (Dh) | Nu local at exit | Nu mean | h channel (W/m2K) | h flat plate (W/m2K) |
|---|---|---|---|---|---|
| 0.3 | 142.0 | 7.54 | 7.82 | 21.7 | 6.8 |
| 1.0 | 474.0 | 7.54 | 8.48 | 23.5 | 12.4 |
| 3.0 | 1421.0 | 7.88 | 10.23 | 28.4 | 21.4 |

### Fin Shapes
| Shape | Height (mm) | Q per fin (W) | Efficiency | Mass (g) | Q per kg (W/kg) |
|---|---|---|---|---|---|
| Rectangular plate | 20.0 | 27.7 | 0.977 | 16.2 | 1710.0 |
| Rectangular plate | 40.0 | 50.4 | 0.921 | 32.4 | 1555.0 |
| Rectangular plate | 60.0 | 68.4 | 0.844 | 48.6 | 1407.0 |
| Triangular plate | 20.0 | 24.9 | 0.971 | 8.1 | 3079.0 |
| Triangular plate | 40.0 | 46.0 | 0.898 | 16.2 | 2839.0 |
| Triangular plate | 60.0 | 61.7 | 0.803 | 24.3 | 2538.0 |
| Cylindrical pin (D=6 mm) | 20.0 | 2.5 | 0.978 | 1.5 | 1661.0 |
| Cylindrical pin (D=6 mm) | 40.0 | 4.6 | 0.923 | 3.1 | 1513.0 |
| Cylindrical pin (D=6 mm) | 60.0 | 6.3 | 0.847 | 4.6 | 1372.0 |

### Geometry
| Height (mm) | Thickness (mm) | Q (W) | Efficiency | Q per area (W/m2) |
|---|---|---|---|---|
| 20.0 | 2.0 | 31.6 | 0.962 | 7515.1 |
| 20.0 | 4.0 | 34.3 | 0.978 | 7792.5 |
| 20.0 | 6.0 | 36.7 | 0.984 | 7986.9 |
| 40.0 | 2.0 | 55.8 | 0.871 | 6807.8 |
| 40.0 | 4.0 | 62.0 | 0.926 | 7375.4 |
| 40.0 | 6.0 | 66.0 | 0.946 | 7679.1 |
| 60.0 | 2.0 | 72.4 | 0.76 | 5937.4 |
| 60.0 | 4.0 | 84.3 | 0.854 | 6801.1 |
| 60.0 | 6.0 | 91.2 | 0.892 | 7240.3 |
| 80.0 | 2.0 | 82.6 | 0.652 | 5096.2 |
| 80.0 | 4.0 | 101.1 | 0.774 | 6163.9 |
| 80.0 | 6.0 | 111.7 | 0.828 | 6726.4 |

### Languages
| Language | Q (W) | Run time incl. startup (ms) | Status |
|---|---|---|---|
| Python (in-process, no process startup) | 85.2048 | 0.24 | ran |
| C | 85.2048 | 2.02 | ran |
| C++ | 85.2048 | 3.68 | ran |
| Rust | - | - | rustc not installed |
| Go | - | - | go not installed |

### Materials
| Material | k (W/mK) | V (m/s) | h (W/m2K) | Q (W) | Efficiency |
|---|---|---|---|---|---|
| Aluminium | 180 | 5 | 27.6 | 43.0 | 0.916 |
| Aluminium | 180 | 15 | 47.9 | 70.3 | 0.865 |
| Aluminium | 180 | 30 | 67.7 | 94.3 | 0.821 |
| Copper | 390 | 5 | 27.6 | 45.0 | 0.959 |
| Copper | 390 | 15 | 47.9 | 75.7 | 0.932 |
| Copper | 390 | 30 | 67.7 | 104.2 | 0.906 |
| Cast iron | 52 | 5 | 27.6 | 36.0 | 0.767 |
| Cast iron | 52 | 15 | 47.9 | 54.0 | 0.664 |
| Cast iron | 52 | 30 | 67.7 | 68.1 | 0.592 |
| Steel | 45 | 5 | 27.6 | 34.8 | 0.742 |
| Steel | 45 | 15 | 47.9 | 51.6 | 0.635 |
| Steel | 45 | 30 | 67.7 | 64.5 | 0.562 |
| Magnesium | 156 | 5 | 27.6 | 42.4 | 0.905 |
| Magnesium | 156 | 15 | 47.9 | 68.9 | 0.848 |
| Magnesium | 156 | 30 | 67.7 | 91.9 | 0.8 |

### Ml Models
| Model | R2 | MAE (W) |
|---|---|---|
| Linear | 0.924 | 6.63 |
| Random Forest | 0.9837 | 2.81 |
| Gradient Boosting | 0.9881 | 2.49 |
| Neural Net (MLP) | 0.9917 | 2.13 |

### Montecarlo
| Mean Q (W) | Std (W) | 2.5 % | 97.5 % | Most sensitive input |
|---|---|---|---|---|
| 69.8 | 13.1 | 44.0 | 95.4 | h |

### Pareto
| dP limit (Pa) | S (mm) | t (mm) | Height (mm) | Q (W) | dP (Pa) | Mass (kg) |
|---|---|---|---|---|---|---|
| 30.0 | 7.26 | 1.5 | 60.0 | 1200.0 | 30.0 | 0.55 |
| 50.0 | 3.56 | 1.5 | 60.0 | 2634.0 | 50.0 | 0.96 |
| 80.0 | 3.0 | 1.5 | 55.6 | 3109.0 | 61.3 | 1.0 |
| 120.0 | 3.0 | 1.5 | 55.6 | 3109.0 | 61.3 | 1.0 |
| 200.0 | 3.0 | 1.5 | 55.6 | 3109.0 | 61.3 | 1.0 |
| 400.0 | 3.0 | 1.5 | 55.6 | 3109.0 | 61.3 | 1.0 |
| 800.0 | 3.0 | 1.5 | 55.6 | 3109.0 | 61.3 | 1.0 |

### Spacing
| Base temp (C) | Optimal S (mm) | Number of fins | h total (W/m2K) | Q total (W) |
|---|---|---|---|---|
| 150.0 | 7.0 | 20.1 | 7.8 | 156.0 |
| 200.0 | 6.6 | 20.8 | 8.3 | 251.0 |
| 250.0 | 6.6 | 20.8 | 8.9 | 352.0 |

### Thermal Hydraulic
| V (m/s) | S (mm) | Re (gap) | h (W/m2K) | Q (W) | dP (Pa) | Pump power (mW) |
|---|---|---|---|---|---|---|
| 2.0 | 4.0 | 1167.0 | 36.9 | 1405.0 | 14.3 | 229.5 |
| 2.0 | 8.0 | 1833.0 | 23.4 | 608.0 | 6.4 | 102.4 |
| 2.0 | 12.0 | 2500.0 | 14.2 | 287.0 | 4.9 | 78.1 |
| 5.0 | 4.0 | 2917.0 | 48.0 | 1789.0 | 72.2 | 2889.9 |
| 5.0 | 8.0 | 4583.0 | 34.5 | 878.0 | 38.1 | 1523.6 |
| 5.0 | 12.0 | 6250.0 | 29.5 | 580.0 | 29.8 | 1192.6 |
| 10.0 | 4.0 | 5833.0 | 83.6 | 2930.0 | 276.9 | 22148.8 |
| 10.0 | 8.0 | 9167.0 | 60.0 | 1463.0 | 149.0 | 11921.3 |
| 10.0 | 12.0 | 12500.0 | 51.3 | 974.0 | 117.6 | 9404.7 |
| 20.0 | 4.0 | 11667.0 | 145.6 | 4642.0 | 1066.6 | 170662.5 |
| 20.0 | 8.0 | 18333.0 | 104.5 | 2375.0 | 584.8 | 93570.8 |
| 20.0 | 12.0 | 25000.0 | 89.3 | 1599.0 | 464.5 | 74320.0 |

### Transient
| Case | Max T (C) | T at 900 s (C) |
|---|---|---|
| V = 1 m/s | 261.5 | 261.5 |
| V = 3 m/s | 222.2 | 222.2 |
| V = 8 m/s | 181.3 | 181.3 |
| V = 20 m/s | 143.6 | 143.6 |
| idle → 20 m/s → idle | 221.8 | 221.8 |

Rust and Go rows in the language table show "not installed" if those compilers were missing when generated. Run times include process startup, except the in-process Python row.

## License
MIT (H.S. Khalsa)
