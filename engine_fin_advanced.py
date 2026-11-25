"""Advanced analyses: channel flow between fins (Graetz), fin-spacing optimisation,
annular fin (Bessel), Monte-Carlo uncertainty, transient warm-up, extra 3D figures.
Run: python engine_fin_advanced.py"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa
from fin_core import *

plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": .3})
np.random.seed(7)

def save(n): plt.tight_layout(); plt.savefig(n, dpi=130); plt.close()
def md_table(df, fn):
    df.to_csv(fn.replace(".md", ".csv"), index=False)
    s = "| " + " | ".join(df.columns) + " |\n|" + "---|"*len(df.columns) + "\n"
    for _, r in df.iterrows(): s += "| " + " | ".join(str(v) for v in r) + " |\n"
    open(fn, "w").write(s)

# ---- 1. Laminar flow + heat transfer in the gap between two fins (Graetz problem)
b, Lw = 0.006, 0.10
res = {U: graetz_channel(U, b, Lw) for U in (0.3, 1.0, 3.0)}
r = res[0.3]
fig, ax = plt.subplots(1, 2, figsize=(13, 4.5), gridspec_kw={"width_ratios": [1, 3]})
ax[0].plot(r["u"], r["y"]*1000, lw=2); ax[0].set(xlabel="u (m/s)", ylabel="y across gap (mm)", title="Velocity profile (parabolic)")
cs = ax[1].contourf(r["x"]*1000, r["y"]*1000, r["theta"].T, 20, cmap="inferno"); fig.colorbar(cs, ax=ax[1], label="θ = (T−T∞)/(Tw−T∞)")
ax[1].set(xlabel="Flow direction x (mm)", ylabel="y (mm)", title="Temperature development between two fins (U=0.3 m/s)"); ax[1].grid(False)
save("fig14_channel_flow_temperature.png")

plt.figure(figsize=(7, 4.5))
for U, d in res.items():
    ok = (d["x_plus"] > 2e-3) & np.isfinite(d["Nu"])
    plt.semilogx(d["x_plus"][ok], d["Nu"][ok], label=f"U={U} m/s (Re={d['Re']:.0f})")
plt.axhline(7.54, color="k", ls="--", label="Fully developed Nu = 7.54"); plt.ylim(0, 25)
plt.xlabel("x⁺ = x / (Dh·Re·Pr)"); plt.ylabel("Local Nusselt number"); plt.title("Graetz problem: Nu development"); plt.legend(); save("fig15_graetz_nusselt.png")
rows = []
for U, d in res.items():
    hch = d["Nu_mean"]*d["kf"]/(2*b); hfp = h_forced(U, Lc=Lw)[0]
    rows.append([U, round(d["Re"]), round(d["Nu"][-1], 2) if np.isfinite(d["Nu"][-1]) else "-", round(d["Nu_mean"], 2), round(hch, 1), round(hfp, 1)])
md_table(pd.DataFrame(rows, columns=["U (m/s)", "Re (Dh)", "Nu local at exit", "Nu mean", "h channel (W/m2K)", "h flat plate (W/m2K)"]), "table_channel_nusselt.md")

fig = plt.figure(figsize=(8, 5.5)); a = fig.add_subplot(111, projection="3d")
sub = r["theta"][::15, ::2]; X, Y = np.meshgrid(r["y"][::2]*1000, r["x"][::15]*1000)
s = a.plot_surface(Y, X, sub, cmap="inferno"); fig.colorbar(s, shrink=.6, label="θ")
a.set(xlabel="x (mm)", ylabel="y (mm)", zlabel="θ", title="3D temperature field in the fin gap"); a.view_init(28, -55); save("fig20_3d_channel_temperature.png")

# ---- 2. Fin spacing optimisation (natural convection + radiation, Elenbaas)
Sg = np.linspace(0.002, 0.03, 80)
Tbs = (150, 200, 250)
curves = {T_: [array_heat(S, Tw=T_)[0] for S in Sg] for T_ in Tbs}
fig, ax = plt.subplots(1, 2, figsize=(13, 4.5))
for T_, q in curves.items():
    ax[0].plot(Sg*1000, q, label=f"base {T_} °C"); i = int(np.argmax(q)); ax[0].plot(Sg[i]*1000, q[i], "ko")
ax[0].set(xlabel="Fin spacing S (mm)", ylabel="Total heat Q (W)", title="Optimal fin spacing (natural convection + radiation)"); ax[0].legend()
Ts_ = Tinf+.85*(Tb-Tinf)
ax[1].plot(Sg*1000, [elenbaas_h(S, Ts=Ts_) for S in Sg], label="natural convection h")
ax[1].plot(Sg*1000, [h_rad(S, Ts=Ts_) for S in Sg], label="radiation h (reduced by view factor)"); ax[1].set(xlabel="S (mm)", ylabel="h (W/m²K)", title="Heat-transfer coefficients vs spacing (base 200 °C)"); ax[1].legend()
save("fig16_fin_spacing_optimum.png")
rows = []
for T_, q in curves.items():
    i = int(np.argmax(q)); Q, n, h = array_heat(Sg[i], Tw=T_)
    rows.append([T_, round(Sg[i]*1000, 1), round(n, 1), round(h, 1), round(Q, 0)])
md_table(pd.DataFrame(rows, columns=["Base temp (C)", "Optimal S (mm)", "Number of fins", "h total (W/m2K)", "Q total (W)"]), "table_spacing.md")

Tg, Sg2 = np.meshgrid(np.linspace(100, 300, 30), np.linspace(.003, .03, 30))
Qg = np.array([[array_heat(s, Tw=tb)[0] for tb, s in zip(rt, rs)] for rt, rs in zip(Tg, Sg2)])
fig = plt.figure(figsize=(8, 5.5)); a = fig.add_subplot(111, projection="3d")
s = a.plot_surface(Tg, Sg2*1000, Qg, cmap="plasma"); fig.colorbar(s, shrink=.6, label="Q (W)")
a.set(xlabel="Base temperature (°C)", ylabel="Spacing (mm)", zlabel="Q (W)", title="3D: array heat vs base temperature and spacing"); a.view_init(28, -130); save("fig22_3d_spacing_temperature.png")

# ---- 3. Annular (engine-cylinder) fin vs straight fin
hs = np.linspace(5, 150, 60)
fig, ax = plt.subplots(1, 2, figsize=(13, 4.5))
for kk in (52, 180, 390): ax[0].plot(hs, [annular_eta(h, kk) for h in hs], label=f"annular, k={kk}")
ax[0].plot(hs, [fin(h, 180, 0.03)[1] for h in hs], "k--", label="straight, k=180"); ax[0].set(xlabel="h (W/m²K)", ylabel="η", title="Fin efficiency"); ax[0].legend()
for r1 in (0.03, 0.04, 0.06): ax[1].plot(hs, [annular_eta(h, 180, 0.003, r1, 0.03) for h in hs], label=f"r1={r1*1000:.0f} mm")
ax[1].set(xlabel="h (W/m²K)", ylabel="η", title="Annular fin, effect of cylinder radius (k=180)"); ax[1].legend(); save("fig17_annular_vs_straight.png")
md_table(pd.DataFrame([[h, round(annular_eta(h), 3), round(fin(h, 180, 0.03)[1], 3), round(annular_Q(h, 160), 1)] for h in (10, 30, 60, 100)],
         columns=["h (W/m2K)", "Annular eta", "Straight eta", "Annular Q per fin (W), theta=160"]), "table_annular.md")
Hg, Lg = np.meshgrid(np.linspace(10, 150, 40), np.linspace(.01, .06, 40))
Eg = np.array([[annular_eta(h, 180, .003, .04, l) for h, l in zip(rh, rl)] for rh, rl in zip(Hg, Lg)])
fig = plt.figure(figsize=(8, 5.5)); a = fig.add_subplot(111, projection="3d")
s = a.plot_surface(Hg, Lg*1000, Eg, cmap="viridis"); fig.colorbar(s, shrink=.6, label="η")
a.set(xlabel="h (W/m²K)", ylabel="Fin height (mm)", zlabel="η", title="3D: annular fin efficiency"); a.view_init(25, -120); save("fig21_3d_annular_efficiency.png")

# ---- 4. Monte-Carlo uncertainty (V = 15 m/s)
n = 20000; h = h_forced(15)[0]*np.random.normal(1, .2, n).clip(.3)
P = pd.DataFrame({"h": h, "k": np.random.normal(180, 9, n), "L": np.random.normal(.05, .001, n), "t": np.random.normal(.003, .00015, n),
                  "Tb": np.random.normal(200, 10, n), "Tinf": np.random.normal(40, 5, n)})
P["Q"] = fin(P.h, P.k, P.L, P.t, P.Tb, P.Tinf)[0]
lo, hi = np.percentile(P.Q, [2.5, 97.5])
fig, ax = plt.subplots(1, 2, figsize=(13, 4.5))
ax[0].hist(P.Q, 60, color="C0", alpha=.8); ax[0].axvline(lo, color="r", ls="--"); ax[0].axvline(hi, color="r", ls="--", label="95 % interval"); ax[0].axvline(P.Q.mean(), color="k")
ax[0].set(xlabel="Q (W)", ylabel="count", title=f"Monte-Carlo: Q = {P.Q.mean():.1f} ± {P.Q.std():.1f} W"); ax[0].legend()
sp = P.corr(method="spearman")["Q"].drop("Q").sort_values(); ax[1].barh(sp.index, sp.values, color="C1"); ax[1].set(title="Sensitivity (Spearman rank correlation with Q)")
save("fig18_monte_carlo.png")
md_table(pd.DataFrame([[round(P.Q.mean(), 1), round(P.Q.std(), 1), round(lo, 1), round(hi, 1), sp.index[-1]]], columns=["Mean Q (W)", "Std (W)", "2.5 %", "97.5 %", "Most sensitive input"]), "table_montecarlo.md")

# ---- 5. Transient engine warm-up (lumped head + annular fins)
Cth, Qin, Nf = 3000.0, 1000.0, 8
def run(Vfun, tend=900):
    T = 40.0; out = []
    for tt in range(tend):
        V = Vfun(tt); hh = h_forced(max(V, .5), Tw=max(T, 41))[0]
        T += (Qin-Nf*annular_Q(hh, T-Tinf))/Cth; out.append(T)
    return np.array(out)
cases = {"V = 1 m/s": lambda t: 1, "V = 3 m/s": lambda t: 3, "V = 8 m/s": lambda t: 8, "V = 20 m/s": lambda t: 20,
         "idle → 20 m/s → idle": lambda t: 20 if 300 <= t < 600 else .5}
plt.figure(figsize=(8, 4.8)); rows = []
for n_, f in cases.items():
    T = run(f); plt.plot(T, label=n_); rows.append([n_, round(T.max(), 1), round(T[-1], 1)])
plt.axhline(250, color="r", ls=":", label="assumed limit 250 °C"); plt.ylim(30, 420)
plt.xlabel("Time (s)"); plt.ylabel("Cylinder-head temperature (°C)"); plt.title("Warm-up with 8 annular fins, 1 kW heat input"); plt.legend(fontsize=8); save("fig19_transient_warmup.png")
md_table(pd.DataFrame(rows, columns=["Case", "Max T (C)", "T at 900 s (C)"]), "table_transient.md")
print("advanced done")
