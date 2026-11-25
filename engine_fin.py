"""Engine fin convection: analysis, ML, comparison tables, 2D/3D figures.
Run:  python engine_fin.py   (all outputs are written to this same folder)"""
import os, time, subprocess, shutil
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.linear_model import LinearRegression
from sklearn.neural_network import MLPRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score, mean_absolute_error

np.random.seed(42)
from fin_core import *
plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": .3})
w, t0 = w0, t0

def save(name): plt.tight_layout(); plt.savefig(name, dpi=130); plt.close()
def md_table(df, fn):
    df.to_csv(fn.replace(".md", ".csv"), index=False)
    cols = list(df.columns); s = "| " + " | ".join(cols) + " |\n|" + "---|"*len(cols) + "\n"
    for _, r in df.iterrows(): s += "| " + " | ".join(str(v) for v in r) + " |\n"
    open(fn, "w").write(s); return s

# ---------- 1. basic curves ----------
Vs = np.linspace(0.5, 35, 70)
hs = np.array([h_forced(v)[0] for v in Vs])
fig, ax = plt.subplots(1, 3, figsize=(16, 4.5))
ax[0].plot(Vs, hs, lw=2); ax[0].set(title="Heat transfer coefficient h", xlabel="Air velocity (m/s)", ylabel="h (W/m²K)")
for kk, nm in [(120, "k=120"), (180, "k=180 (Al)"), (230, "k=230")]:
    ax[1].plot(Vs, [fin(h_forced(v)[0], kk)[0] for v in Vs], label=nm)
ax[1].set(title="Fin heat rate Q", xlabel="Air velocity (m/s)", ylabel="Q (W)"); ax[1].legend()
ax[2].plot(Vs, [fin(h_forced(v)[0])[1] for v in Vs], color="C3"); ax[2].set(title="Fin efficiency η", xlabel="Air velocity (m/s)", ylabel="η")
save("fig01_h_Q_eta_curves.png")

# ---------- 2. temperature profile + numeric vs analytic ----------
h60 = 60.0; Qfd, th = fd_fin(h60); xs = np.linspace(0, L0, len(th))*1000
m = fin(h60)[2]; Lc = L0+t0/2
ana = (Tb-Tinf)*np.cosh(m*(Lc-xs/1000))/np.cosh(m*Lc)+Tinf
plt.figure(figsize=(7, 4.5)); plt.plot(xs, th+Tinf, lw=2, label="Finite difference"); plt.plot(xs, ana, "--", label="Analytical")
plt.xlabel("Fin length (mm)"); plt.ylabel("Temperature (°C)"); plt.title("Fin temperature profile (h=60)"); plt.legend(); save("fig02_temperature_profile.png")

# ---------- 3. comparison tables ----------
mats = {"Aluminium": 180, "Copper": 390, "Cast iron": 52, "Steel": 45, "Magnesium": 156}
rows = []
for n, kk in mats.items():
    for V in (5, 15, 30):
        h, Re = h_forced(V); Q, e, _ = fin(h, kk); rows.append([n, kk, V, round(h, 1), round(Q, 1), round(e, 3)])
mat_df = pd.DataFrame(rows, columns=["Material", "k (W/mK)", "V (m/s)", "h (W/m2K)", "Q (W)", "Efficiency"])
md_table(mat_df, "table_materials.md")
geo = []
for L in (0.02, 0.04, 0.06, 0.08):
    for t in (0.002, 0.004, 0.006):
        h = h_forced(15)[0]; Q, e, _ = fin(h, L=L, t=t); geo.append([L*1000, t*1000, round(Q, 1), round(e, 3), round(Q/(2*w*L+w*t)/1, 1)])
geo_df = pd.DataFrame(geo, columns=["Height (mm)", "Thickness (mm)", "Q (W)", "Efficiency", "Q per area (W/m2)"])
md_table(geo_df, "table_geometry.md")
piv = mat_df.pivot(index="Material", columns="V (m/s)", values="Q (W)")
piv.plot(kind="bar", figsize=(8, 4.5)); plt.ylabel("Q (W)"); plt.title("Material comparison (bar chart)"); plt.xticks(rotation=0); save("fig03_material_bar.png")
hm = geo_df.pivot(index="Height (mm)", columns="Thickness (mm)", values="Efficiency")
plt.figure(figsize=(6, 4.5)); plt.imshow(hm.values, cmap="viridis", origin="lower", aspect="auto"); plt.colorbar(label="η")
plt.xticks(range(3), hm.columns); plt.yticks(range(4), hm.index); plt.xlabel("Thickness (mm)"); plt.ylabel("Height (mm)")
for i in range(4):
    for j in range(3): plt.text(j, i, f"{hm.values[i,j]:.2f}", ha="center", color="w")
plt.title("Efficiency heat-map"); plt.grid(False); save("fig04_efficiency_heatmap.png")
labels = ["Air-cooled\n(natural)", "Forced 5 m/s", "Forced 15 m/s", "Forced 30 m/s"]
hh = [8, h_forced(5)[0], h_forced(15)[0], h_forced(30)[0]]
plt.figure(figsize=(7, 4.5)); plt.pie([fin(x)[0] for x in hh], labels=labels, autopct="%1.0f%%"); plt.title("Share of Q by cooling condition (pie)"); save("fig05_cooling_pie.png")

# ---------- 4. data science + ML ----------
N = 1500; V = np.random.uniform(.5, 35, N); Lf = np.random.uniform(.02, .08, N)
tf = np.random.uniform(.002, .006, N); kf = np.random.choice([120, 180, 230], N)
rec = []
for v, l, tt, kk in zip(V, Lf, tf, kf):
    h, Re = h_forced(v); Q, e, _ = fin(h, kk, l, tt); rec.append((v, l, tt, kk, Re, h, Q*(1+np.random.normal(0, .03)), e))
df = pd.DataFrame(rec, columns=["V", "L", "t", "k", "Re", "h", "Q", "eta"]); df.to_csv("fin_data_clean.csv", index=False)
X = df[["V", "L", "t", "k"]]; y = df["Q"]; Xtr, Xte, ytr, yte = train_test_split(X, y, test_size=.2, random_state=0)
sx, sy = StandardScaler().fit(Xtr), StandardScaler().fit(ytr.values.reshape(-1, 1))
models = {"Linear": LinearRegression(), "Random Forest": RandomForestRegressor(200, random_state=0),
          "Gradient Boosting": GradientBoostingRegressor(random_state=0)}
res, preds = [], {}
for n, mdl in models.items():
    mdl.fit(Xtr, ytr); p = mdl.predict(Xte); preds[n] = p
    res.append([n, round(r2_score(yte, p), 4), round(mean_absolute_error(yte, p), 2)])
nn = MLPRegressor(hidden_layer_sizes=(64, 64), max_iter=2000, random_state=0).fit(sx.transform(Xtr), sy.transform(ytr.values.reshape(-1, 1)).ravel())
p = sy.inverse_transform(nn.predict(sx.transform(Xte)).reshape(-1, 1)).ravel(); preds["Neural Net (MLP)"] = p
res.append(["Neural Net (MLP)", round(r2_score(yte, p), 4), round(mean_absolute_error(yte, p), 2)])
ml_df = pd.DataFrame(res, columns=["Model", "R2", "MAE (W)"]); md_table(ml_df, "table_ml_models.md")
fig, ax = plt.subplots(1, 2, figsize=(13, 4.5))
ax[0].bar(ml_df["Model"], ml_df["R2"], color=["C0", "C1", "C2", "C3"]); ax[0].set(title="ML model comparison (R²)", ylim=(0, 1)); ax[0].tick_params(axis="x", rotation=15)
for n, p in preds.items(): ax[1].scatter(yte, p, s=6, alpha=.5, label=n)
lim = [yte.min(), yte.max()]; ax[1].plot(lim, lim, "k--"); ax[1].set(title="Predicted vs actual Q", xlabel="Actual (W)", ylabel="Predicted (W)"); ax[1].legend(markerscale=3)
save("fig06_ml_comparison.png")
cm = df[["V", "L", "t", "k", "h", "Q", "eta"]].corr()
plt.figure(figsize=(6, 5)); plt.imshow(cm, cmap="coolwarm", vmin=-1, vmax=1); plt.colorbar(); plt.xticks(range(7), cm.columns, rotation=45); plt.yticks(range(7), cm.columns)
plt.grid(False); plt.title("Correlation matrix"); save("fig07_correlation.png")
fi = RandomForestRegressor(200, random_state=0).fit(Xtr, ytr).feature_importances_
plt.figure(figsize=(6, 4)); plt.barh(X.columns, fi, color="C4"); plt.title("Feature importance (Random Forest)"); save("fig08_feature_importance.png")

# ---------- 5. language comparison (C, C++, Rust, Go, Python) ----------
def timeit(cmd):
    t = time.perf_counter(); o = subprocess.run(cmd, capture_output=True, text=True, shell=True); return o.stdout.strip(), (time.perf_counter()-t)*1000
lang = []
t = time.perf_counter(); q, _ = fd_fin(60); lang.append(["Python (in-process, no process startup)", round(q, 4), round((time.perf_counter()-t)*1000, 2), "ran"])
build = {"C": ("gcc", "gcc -O2 fin.c -o fin_c -lm", "./fin_c 60"), "C++": ("g++", "g++ -O2 fin.cpp -o fin_cpp", "./fin_cpp 60"),
         "Rust": ("rustc", "rustc -O fin.rs -o fin_rs", "./fin_rs 60"), "Go": ("go", "go build -o fin_go fin.go", "./fin_go 60")}
for n, (tool, b, r) in build.items():
    if shutil.which(tool):
        subprocess.run(b, shell=True, capture_output=True); o, ms = timeit(r)
        try: lang.append([n, float(o), round(ms, 2), "ran"])
        except ValueError: lang.append([n, "-", "-", "build failed"])
    else: lang.append([n, "-", "-", f"{tool} not installed"])
lang_df = pd.DataFrame(lang, columns=["Language", "Q (W)", "Run time incl. startup (ms)", "Status"]); md_table(lang_df, "table_languages.md")
ok = lang_df[lang_df.Status == "ran"]
plt.figure(figsize=(7, 4)); plt.bar(ok["Language"], ok["Run time incl. startup (ms)"].astype(float), color="C2"); plt.ylabel("ms"); plt.title("Language run-time comparison (solver, N=200)"); save("fig09_language_runtime.png")
for f in ("fin_c", "fin_cpp", "fin_rs", "fin_go"):
    if os.path.exists(f): os.remove(f)

# ---------- 6. 2D field + 3D diagrams ----------
nx, ny = 120, 60; T = np.full((ny, nx), Tinf); T[:, 0] = Tb; dx = L0/nx*2; Bi = 60*dx/k0
for _ in range(6000):
    Tn = T.copy(); Tn[1:-1, 1:-1] = (T[2:, 1:-1]+T[:-2, 1:-1]+T[1:-1, 2:]+T[1:-1, :-2]+Bi*Tinf*.3)/(4+Bi*.3)
    Tn[:, 0] = Tb; Tn[0, :] = Tn[1, :]; Tn[-1, :] = Tn[-2, :]; Tn[:, -1] = Tn[:, -2]; T = Tn
X2, Y2 = np.meshgrid(np.linspace(0, L0*1000, nx), np.linspace(0, w*1000, ny))
fig = plt.figure(figsize=(8, 5.5)); a = fig.add_subplot(111, projection="3d")
s = a.plot_surface(X2, Y2, T, cmap="inferno", linewidth=0); fig.colorbar(s, shrink=.6, label="T (°C)")
a.set(xlabel="Length (mm)", ylabel="Width (mm)", zlabel="T (°C)", title="3D temperature surface of fin"); a.view_init(30, -60); save("fig10_3d_temperature_surface.png")

Vg, Lg = np.meshgrid(np.linspace(1, 35, 40), np.linspace(.02, .08, 40))
Qg = np.array([[fin(h_forced(v)[0], L=l)[0] for v, l in zip(rv, rl)] for rv, rl in zip(Vg, Lg)])
fig = plt.figure(figsize=(8, 5.5)); a = fig.add_subplot(111, projection="3d")
s = a.plot_surface(Vg, Lg*1000, Qg, cmap="viridis"); fig.colorbar(s, shrink=.6, label="Q (W)")
a.set(xlabel="Velocity (m/s)", ylabel="Fin height (mm)", zlabel="Q (W)", title="3D surface: Q vs velocity and height"); a.view_init(28, -125); save("fig11_3d_Q_surface.png")

fig = plt.figure(figsize=(8, 5.5)); a = fig.add_subplot(111, projection="3d")
Vb = [5, 15, 30]; xs_ = np.arange(len(mats))
for j, v in enumerate(Vb):
    qs = [fin(h_forced(v)[0], kk)[0] for kk in mats.values()]
    a.bar3d(xs_, [j]*len(xs_), np.zeros(len(xs_)), .6, .6, qs, alpha=.85, color=f"C{j}")
a.set_xticks(xs_+.3); a.set_xticklabels(list(mats), rotation=-20, fontsize=7); a.set_yticks(np.arange(3)+.3); a.set_yticklabels(Vb)
a.set(ylabel="Velocity (m/s)", zlabel="Q (W)", title="3D bar chart: material × velocity"); save("fig12_3d_bar_materials.png")

fig = plt.figure(figsize=(8, 6)); a = fig.add_subplot(111, projection="3d")
th_, zz = np.meshgrid(np.linspace(0, 2*np.pi, 60), np.linspace(0, .12, 2)); r0 = .04
a.plot_surface(r0*np.cos(th_), r0*np.sin(th_), zz, color="gray", alpha=.7)
mm = 60.0; Lc = .03  # illustrative m (1/m) to make the gradient visible
rr, tt2 = np.meshgrid(np.linspace(r0, r0+Lc, 20), np.linspace(0, 2*np.pi, 60))
Tfin = Tinf+(Tb-Tinf)*np.cosh(mm*(Lc-(rr-r0)))/np.cosh(mm*Lc)
cmap = plt.cm.inferno; nrm = plt.Normalize(Tinf, Tb)
for z in np.linspace(.01, .11, 6):
    a.plot_surface(rr*np.cos(tt2), rr*np.sin(tt2), np.full_like(rr, z), facecolors=cmap(nrm(Tfin)), shade=False)
a.set(title="3D engine cylinder with circumferential fins (colour = temperature)", xlabel="x (m)", ylabel="y (m)", zlabel="z (m)")
a.set_box_aspect((1, 1, .8)); a.view_init(25, -50)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=nrm); fig.colorbar(sm, ax=a, shrink=.6, label="T (°C)"); save("fig13_3d_engine_fins.png")
print("done"); print(ml_df); print(lang_df)
