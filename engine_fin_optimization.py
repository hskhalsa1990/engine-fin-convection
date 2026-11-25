"""Fin shapes, thermal-hydraulic trade-off, multi-objective optimisation, animated 2D transient.
Run: python engine_fin_optimization.py"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from mpl_toolkits.mplot3d import Axes3D  # noqa
from scipy.optimize import differential_evolution
from scipy.special import iv
from fin_core import *

plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": .3})
def save(n): plt.tight_layout(); plt.savefig(n, dpi=130); plt.close()
def md_table(df, fn):
    df.to_csv(fn.replace(".md", ".csv"), index=False)
    s = "| " + " | ".join(df.columns) + " |\n|" + "---|"*len(df.columns) + "\n"
    for _, r in df.iterrows(): s += "| " + " | ".join(str(v) for v in r) + " |\n"
    open(fn, "w").write(s)

# ---- 1. Fin shape comparison: rectangular, triangular, pin (same base thickness/diameter scale)
h = 40.0; t = 0.003; D = 0.006; Ls = np.linspace(0.01, 0.08, 50)
def rect(L): Q, e, _ = fin(h, k0, L, t); return Q, e, RHO_AL*w0*t*L
shapes = {"Rectangular plate": rect, "Triangular plate": lambda L: triangular_fin(h, k0, L, t), "Cylindrical pin (D=6 mm)": lambda L: pin_fin(h, k0, L, D)}
fig, ax = plt.subplots(1, 3, figsize=(16, 4.5)); rows = []
for n, f in shapes.items():
    Q, e, m = np.array([f(L) for L in Ls]).T
    ax[0].plot(Ls*1000, Q, label=n); ax[1].plot(Ls*1000, e, label=n); ax[2].plot(Ls*1000, Q/m, label=n)
    for L in (0.02, 0.04, 0.06):
        q, ee, mm = f(L); rows.append([n, L*1000, round(q, 1), round(ee, 3), round(mm*1000, 1), round(q/mm, 0)])
ax[0].set(xlabel="Fin height (mm)", ylabel="Q per fin (W)", title="Heat rate"); ax[1].set(xlabel="Fin height (mm)", ylabel="η", title="Efficiency")
ax[2].set(xlabel="Fin height (mm)", ylabel="Q / mass (W/kg)", title="Specific heat rate (h = 40, k = 180)"); ax[0].legend(fontsize=8)
save("fig23_fin_shapes_comparison.png")
md_table(pd.DataFrame(rows, columns=["Shape", "Height (mm)", "Q per fin (W)", "Efficiency", "Mass (g)", "Q per kg (W/kg)"]), "table_fin_shapes.md")

# ---- 2. Thermal-hydraulic trade-off (forced flow between fins)
Sg = np.linspace(0.002, 0.02, 80); V = 5.0
R = [array_forced(V, s) for s in Sg]
fig, ax = plt.subplots(1, 3, figsize=(16, 4.5))
ax[0].plot(Sg*1000, [r["Q"] for r in R]); ax[0].set(xlabel="Spacing S (mm)", ylabel="Q (W)", title=f"Total heat (V = {V:.0f} m/s)")
ax[1].semilogy(Sg*1000, [r["dp"] for r in R], color="C3"); ax[1].set(xlabel="Spacing S (mm)", ylabel="Pressure drop (Pa)", title="Pressure drop")
ax[2].plot(Sg*1000, [r["Q"]/max(r["pump"], 1e-9) for r in R], color="C2"); ax[2].set(xlabel="Spacing S (mm)", ylabel="Q / pumping power (-)", title="Thermal-hydraulic performance")
save("fig24_thermal_hydraulic.png")
rows = []
for Vv in (2, 5, 10, 20):
    for s in (0.004, 0.008, 0.012):
        r = array_forced(Vv, s); rows.append([Vv, s*1000, round(r["Re"]), round(r["h"], 1), round(r["Q"]), round(r["dp"], 1), round(r["pump"]*1000, 1)])
md_table(pd.DataFrame(rows, columns=["V (m/s)", "S (mm)", "Re (gap)", "h (W/m2K)", "Q (W)", "dP (Pa)", "Pump power (mW)"]), "table_thermal_hydraulic.md")

# ---- 3. Multi-objective optimisation: max Q subject to pressure-drop and mass limits (V = 10 m/s)
def optimise(dp_max, m_max=1.0, V=5.0):
    def obj(x):
        r = array_forced(V, x[0], x[1], x[2])
        return -r["Q"]+1e5*max(0, r["dp"]/dp_max-1)+1e5*max(0, r["mass"]/m_max-1)
    sol = differential_evolution(obj, [(0.003, 0.02), (0.0015, 0.006), (0.01, 0.06)], seed=1, maxiter=80, popsize=15, tol=1e-8, polish=False)
    return sol.x, array_forced(V, *sol.x)
rows, P = [], []
for dpm in (30, 50, 80, 120, 200, 400, 800):
    x, r = optimise(dpm); rows.append([dpm, round(x[0]*1000, 2), round(x[1]*1000, 2), round(x[2]*1000, 1), round(r["Q"]), round(r["dp"], 1), round(r["mass"], 2)])
    P.append((r["Q"], r["dp"], r["mass"], x[0]))
par = pd.DataFrame(rows, columns=["dP limit (Pa)", "S (mm)", "t (mm)", "Height (mm)", "Q (W)", "dP (Pa)", "Mass (kg)"]); md_table(par, "table_pareto.md")
fig, ax = plt.subplots(1, 2, figsize=(13, 4.5))
ax[0].plot(par["dP (Pa)"], par["Q (W)"], "o-"); ax[0].set(xlabel="Pressure drop (Pa)", ylabel="Maximum Q (W)", title="Pareto front: heat vs pressure drop (V = 5 m/s, mass ≤ 1 kg)")
for c in ("S (mm)", "t (mm)"): ax[1].plot(par["dP limit (Pa)"], par[c], "o-", label=c)
ax[1].plot(par["dP limit (Pa)"], par["Height (mm)"]/10, "o-", label="Height / 10 (mm)"); ax[1].set_xscale("log"); ax[1].legend()
ax[1].set(xlabel="Allowed pressure drop (Pa)", ylabel="mm", title="Optimal fin geometry"); save("fig25_pareto_optimisation.png")
P = np.array(P); fig = plt.figure(figsize=(8, 5.5)); a = fig.add_subplot(111, projection="3d")
s = a.scatter(P[:, 1], P[:, 2], P[:, 0], c=P[:, 3]*1000, s=70, cmap="viridis"); fig.colorbar(s, shrink=.6, label="optimal S (mm)")
a.plot(P[:, 1], P[:, 2], P[:, 0], "k-", alpha=.4); a.set(xlabel="dP (Pa)", ylabel="Mass (kg)", zlabel="Q (W)", title="3D Pareto points: Q, pressure drop and mass"); a.view_init(25, -125); save("fig26_3d_pareto.png")

# ---- 4. 3D fin shapes (voxels) coloured by temperature (illustrative h = 150 so the gradient is visible)
nz, nc = 30, 20; fig = plt.figure(figsize=(14, 5)); hh = 150.0; L = 0.04
cx = (nc-1)/2; X_, Y_ = np.meshgrid(np.arange(nc), np.arange(nc), indexing="ij")
for i, (name, kind) in enumerate([("Rectangular", "r"), ("Triangular", "t"), ("Pin", "p")]):
    vox = np.zeros((nc, nc, nz), bool); th = np.zeros((nc, nc, nz))
    for z in range(nz):
        zz = (z+.5)/nz*L
        if kind == "r":
            m = np.sqrt(2*hh/(k0*t)); T = np.cosh(m*(L-zz))/np.cosh(m*L); vox[:, :, z] = np.abs(X_-cx) <= 2.5
        elif kind == "t":
            m = np.sqrt(2*hh/(k0*t)); T = iv(0, 2*m*np.sqrt(L*(L-zz)))/iv(0, 2*m*L); vox[:, :, z] = np.abs(X_-cx) <= 2.5*(1-zz/L)+0.5
        else:
            m = np.sqrt(4*hh/(k0*D)); T = np.cosh(m*(L-zz))/np.cosh(m*L); vox[:, :, z] = (X_-cx)**2+(Y_-cx)**2 <= 3.0**2
        th[:, :, z] = T
    cols = plt.cm.inferno(np.clip((th-0.2)/0.8, 0, 1))
    a = fig.add_subplot(1, 3, i+1, projection="3d"); a.voxels(vox, facecolors=cols, edgecolor=None)
    a.set(title=name+" fin"); a.set_xticks([]); a.set_yticks([]); a.set_zticks([]); a.view_init(22, -60); a.set_box_aspect((1, 1, 1.5))
fig.suptitle("3D fin shapes coloured by temperature (bright = hot, base at bottom)"); save("fig27_3d_fin_shapes.png")

# ---- 5. 2D transient conduction in a fin plate with a hot spot at the base (animated)
Lf, Wf, tf, hc = 0.05, 0.10, 0.003, 40.0; nxx, nyy = 50, 100; dx = Lf/nxx; dy = Wf/nyy
al = k0/(RHO_AL*900.0); dt = 0.2*min(dx, dy)**2/al; loss = 2*hc/(RHO_AL*900.0*tf)
yy = np.linspace(0, Wf, nyy); Tbase = Tinf+(Tb-Tinf)*(0.75+0.25*np.exp(-((yy-Wf/2)/0.015)**2))
T = np.full((nxx, nyy), Tinf); T[0, :] = Tbase; frames = []; tmax = 40.0; nst = int(tmax/dt); every = nst//40
for s in range(nst+1):
    Tn = T.copy()
    Tn[1:-1, 1:-1] = T[1:-1, 1:-1]+dt*(al*((T[2:, 1:-1]-2*T[1:-1, 1:-1]+T[:-2, 1:-1])/dx**2+(T[1:-1, 2:]-2*T[1:-1, 1:-1]+T[1:-1, :-2])/dy**2)-loss*(T[1:-1, 1:-1]-Tinf))
    Tn[-1, 1:-1] = Tn[-2, 1:-1]; Tn[:, 0] = Tn[:, 1]; Tn[:, -1] = Tn[:, -2]; Tn[0, :] = Tbase; T = Tn
    if s % every == 0: frames.append((s*dt, T.copy()))
fig, ax = plt.subplots(figsize=(6, 3.4)); im = ax.imshow(frames[0][1], origin="lower", extent=[0, Wf*1000, 0, Lf*1000], vmin=Tinf, vmax=Tb, cmap="inferno", aspect="auto")
fig.colorbar(im, label="T (°C)"); ax.set(xlabel="Fin width (mm)", ylabel="Height from base (mm)"); ax.grid(False); ttl = ax.set_title("")
def upd(i): im.set_data(frames[i][1]); ttl.set_text(f"Fin warm-up with a hot spot at the base, t = {frames[i][0]:.1f} s"); return im, ttl
FuncAnimation(fig, upd, frames=len(frames)).save("anim_fin_transient.gif", writer=PillowWriter(fps=8)); plt.close(fig)
fig, ax = plt.subplots(1, 4, figsize=(16, 3.4)); idx = [1, len(frames)//5, len(frames)//2, len(frames)-1]
for a, i in zip(ax, idx):
    c = a.imshow(frames[i][1], origin="lower", extent=[0, Wf*1000, 0, Lf*1000], vmin=Tinf, vmax=Tb, cmap="inferno", aspect="auto"); a.set_title(f"t = {frames[i][0]:.1f} s"); a.grid(False)
fig.colorbar(c, ax=ax, label="T (°C)"); fig.savefig("fig28_transient_2d_snapshots.png", dpi=130, bbox_inches="tight"); plt.close(fig)
print("optimisation done")
