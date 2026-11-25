"""Core physics for engine-fin convection (no plotting, safe to import)."""
import numpy as np
from scipy.special import iv, kv
from scipy.linalg import solve_banded

k0, Tb, Tinf = 180.0, 200.0, 40.0        # W/mK, deg C, deg C
L0, w0, t0 = 0.05, 0.10, 0.003           # fin height, width, thickness (m)
_trapz = getattr(np, "trapezoid", None) or np.trapz
SIGMA = 5.670374419e-8

def air_props(T):
    """Approximate air properties at temperature T (deg C): nu, k, Pr."""
    Tk = T + 273.15
    return 1.516e-5*(Tk/293.15)**1.75, 0.0263*(Tk/293.15)**0.8, 0.707

def h_forced(V, Lc=0.10, Tw=Tb, Ti=Tinf):
    """Flat-plate forced convection h (laminar/turbulent) and Re."""
    nu, kf, Pr = air_props((Tw+Ti)/2); Re = max(V, 1e-6)*Lc/nu
    Nu = 0.664*Re**.5*Pr**(1/3) if Re < 5e5 else (0.037*Re**.8-871)*Pr**(1/3)
    return Nu*kf/Lc, Re

def fin(h, k=k0, L=L0, t=t0, Tw=Tb, Ti=Tinf, w=w0):
    """Straight rectangular fin with corrected length: returns Q, efficiency, m."""
    P = 2*(w+t); Ac = w*t; m = np.sqrt(h*P/(k*Ac)); Lc = L+t/2
    return np.sqrt(h*P*k*Ac)*(Tw-Ti)*np.tanh(m*Lc), np.tanh(m*Lc)/(m*Lc), m

def fd_fin(h, N=200, k=k0, L=L0, t=t0, w=w0, Tw=Tb, Ti=Tinf):
    """1D finite-difference fin (convective tip). Returns Q and theta(x)."""
    P = 2*(w+t); Ac = w*t; m2 = h*P/(k*Ac); dx = L/N
    a = np.ones(N+1); b = np.full(N+1, -(2+m2*dx*dx)); c = np.ones(N+1); d = np.zeros(N+1)
    b[0] = 1; c[0] = 0; d[0] = Tw-Ti; a[0] = 0
    a[N] = 2; b[N] = -(2+m2*dx*dx+2*dx*h/k); c[N] = 0
    for i in range(1, N+1):
        f = a[i]/b[i-1]; b[i] -= f*c[i-1]; d[i] -= f*d[i-1]
    x = np.zeros(N+1); x[N] = d[N]/b[N]
    for i in range(N-1, -1, -1): x[i] = (d[i]-c[i]*x[i+1])/b[i]
    return -k*Ac*(-3*x[0]+4*x[1]-x[2])/(2*dx), x

def annular_eta(h, k=k0, t=t0, r1=0.04, L=0.03):
    """Efficiency of an annular (circumferential) fin, Bessel-function solution."""
    m = np.sqrt(2*h/(k*t)); r2 = r1+L+t/2
    num = 2*r1/(m*(r2**2-r1**2))*(kv(1, m*r1)*iv(1, m*r2)-iv(1, m*r1)*kv(1, m*r2))
    den = iv(0, m*r1)*kv(1, m*r2)+kv(0, m*r1)*iv(1, m*r2)
    return num/den

def annular_Q(h, theta, k=k0, t=t0, r1=0.04, L=0.03):
    r2 = r1+L+t/2
    return annular_eta(h, k, t, r1, L)*h*2*np.pi*(r2**2-r1**2)*theta

def graetz_channel(U, b=0.006, L=0.10, Ts=Tb, Ti=Tinf, ny=81, nx=1500):
    """Laminar parallel-plate channel (fluid between two fins), isothermal walls.
    Parabolic velocity, implicit marching in x. Returns dict with x, y, theta(x,y),
    Nu_local, x_plus, Nu_mean, u(y)."""
    nu, kf, Pr = air_props((Ts+Ti)/2); alpha = nu/Pr; Dh = 2*b
    y = np.linspace(0, b, ny); dy = y[1]-y[0]
    u = 1.5*U*(1-(2*y/b-1)**2); ui = u[1:-1]
    th = np.zeros(ny); th[0] = th[-1] = 1.0
    dx = L/nx; r = alpha*dx/(ui*dy**2); n = ny-2
    ab = np.zeros((3, n)); ab[0, 1:] = -r[:-1]; ab[1] = 1+2*r; ab[2, :-1] = -r[1:]
    TH = np.zeros((nx+1, ny)); TH[0] = th
    for s in range(1, nx+1):
        rhs = th[1:-1].copy(); rhs[0] += r[0]; rhs[-1] += r[-1]
        th[1:-1] = solve_banded((1, 1), ab, rhs); TH[s] = th
    xs = np.linspace(0, L, nx+1)
    thm = _trapz(u*TH, y, axis=1)/_trapz(u, y)
    grad = -(-3*TH[:, 0]+4*TH[:, 1]-TH[:, 2])/(2*dy)
    with np.errstate(divide="ignore", invalid="ignore"):
        Nu = grad*Dh/(1-thm)
    Re = U*Dh/nu
    Nu_mean = (kf/alpha*U*b*np.log(1/(1-thm[-1]))/(2*L))*Dh/kf
    return dict(x=xs, y=y, theta=TH, Nu=Nu, x_plus=xs/(Dh*Re*Pr), Nu_mean=Nu_mean,
                u=u, Re=Re, thm=thm, kf=kf)

def elenbaas_h(S, H=w0, Ts=Tb, Ti=Tinf):
    """Natural-convection h between vertical isothermal parallel plates (Elenbaas)."""
    nu, kf, Pr = air_props((Ts+Ti)/2); alpha = nu/Pr
    Ra = 9.81/((Ts+Ti)/2+273.15)*(Ts-Ti)*S**4/(nu*alpha*H)
    Nu = Ra/24*(1-np.exp(-35/Ra))**0.75
    return Nu*kf/S

def h_rad(S, H=w0, Ts=Tb, Ti=Tinf, eps=0.9):
    """Linearised radiation h to ambient, reduced by fin-to-fin view factor."""
    F = np.sqrt(1+(S/H)**2)-S/H
    Tk, Ta = Ts+273.15, Ti+273.15
    return eps*SIGMA*(Tk**2+Ta**2)*(Tk+Ta)*(1-F)

def array_heat(S, t=t0, L=0.04, Hb=0.20, k=k0, Tw=Tb, Ti=Tinf):
    """Total natural-convection + radiation heat from a vertical fin array on a base
    of height Hb with spacing S (Elenbaas). Fin surface assumed at Ti + 0.85*(Tw-Ti)
    for the h correlations. Forced flow is not included (it would not be a function of S
    in this simple model). Returns Q_total, number of fins, h_total."""
    Ts = Ti+0.85*(Tw-Ti); n = Hb/(S+t)
    h = elenbaas_h(S, Ts=Ts, Ti=Ti)+h_rad(S, Ts=Ts, Ti=Ti)
    Qf = fin(h, k, L, t, Tw, Ti)[0]
    return n*(Qf+h*S*w0*(Tw-Ti)), n, h

# ---------------- v3: other fin shapes, hydraulics, thermal-hydraulic array ----------------
RHO_AL = 2700.0

def triangular_fin(h, k=k0, L=0.04, t=t0, w=w0, Tw=Tb, Ti=Tinf):
    """Triangular (tapered) plate fin. Returns Q, efficiency, mass of one fin (aluminium)."""
    m = np.sqrt(2*h/(k*t)); x = 2*m*L
    eta = iv(1, x)/(m*L*iv(0, x)); As = 2*w*np.sqrt(L**2+(t/2)**2)
    return eta*h*As*(Tw-Ti), eta, RHO_AL*w*t*L/2

def pin_fin(h, k=k0, L=0.04, D=0.006, Tw=Tb, Ti=Tinf):
    """Cylindrical pin fin with corrected length. Returns Q, efficiency, mass."""
    m = np.sqrt(4*h/(k*D)); Lc = L+D/4; Q = np.sqrt(h*np.pi*D*k*np.pi*D**2/4)*(Tw-Ti)*np.tanh(m*Lc)
    return Q, np.tanh(m*Lc)/(m*Lc), RHO_AL*np.pi*D**2/4*L

def channel_forced(V, S, t=t0, H=w0, Ts=Tb, Ti=Tinf):
    """Forced flow in the gap S between two parallel fins (flow along width H).
    V = approach velocity; channel velocity from continuity. Laminar: Nu = 7.54 + entrance
    term (Edwards), f = 96/Re. Turbulent (Re > 2300): Dittus-Boelter, f = 0.316 Re^-0.25.
    Pressure drop includes ~1.5 velocity heads of inlet/outlet losses (approximate)."""
    U = V*(S+t)/S; Tf = (Ts+Ti)/2
    nu, kf, Pr = air_props(Tf); rho = 101325/(287.0*(Tf+273.15)); Dh = 2*S; Re = U*Dh/nu
    if Re < 2300:
        Gz = Dh/H*Re*Pr; Nu = 7.54+0.03*Gz/(1+0.016*Gz**(2/3)); f = 96/Re
    else:
        Nu = 0.023*Re**0.8*Pr**0.4; f = 0.316*Re**-0.25
    dp = (1.5+f*H/Dh)*rho*U**2/2
    return Nu*kf/Dh, dp, Re, U

def array_forced(V, S, t=t0, L=0.04, Hb=0.20, k=k0, Tw=Tb, Ti=Tinf):
    """Fin array in forced flow, ignoring bypass around the array.
    Returns dict: Q (W), dp (Pa), mass (kg), n fins, h, Re, pump power (W)."""
    n = Hb/(S+t); h, dp, Re, U = channel_forced(V, S, t, Ts=Ti+0.85*(Tw-Ti), Ti=Ti)
    Qf = fin(h, k, L, t, Tw, Ti)[0]; Q = n*(Qf+h*S*w0*(Tw-Ti))
    return dict(Q=Q, dp=dp, mass=n*RHO_AL*w0*t*L, n=n, h=h, Re=Re, pump=dp*V*Hb*L)
