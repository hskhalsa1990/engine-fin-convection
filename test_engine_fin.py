"""Run with: pytest -q   (or: python test_engine_fin.py)"""
import numpy as np
from fin_core import fin, fd_fin, h_forced, annular_eta, graetz_channel, array_heat, triangular_fin, pin_fin, channel_forced, array_forced

def test_zero_h_gives_zero_heat():
    assert fin(1e-9)[0] < 1e-3

def test_high_conductivity_efficiency_near_one():
    assert fin(30, k=1e5)[1] > 0.999

def test_fd_matches_analytical():
    h = 60; assert abs(fd_fin(h)[0]-fin(h)[0])/fin(h)[0] < 0.01

def test_h_increases_with_velocity():
    assert h_forced(20)[0] > h_forced(5)[0] > h_forced(1)[0]

def test_annular_efficiency_bounds():
    for h in (5, 50, 150): assert 0 < annular_eta(h) < 1

def test_graetz_fully_developed_nusselt():
    r = graetz_channel(0.3); assert abs(r["Nu"][-1]-7.54)/7.54 < 0.03

def test_spacing_has_interior_optimum():
    S = np.linspace(0.002, 0.03, 60); q = [array_heat(s)[0] for s in S]
    assert 0 < int(np.argmax(q)) < len(S)-1


def test_triangular_and_pin_bounds():
    for h in (10, 60, 120):
        assert 0 < triangular_fin(h)[1] < 1 and 0 < pin_fin(h)[1] < 1

def test_triangular_lighter_than_rectangular_for_same_height():
    from fin_core import RHO_AL, w0, t0
    assert triangular_fin(40)[2] < RHO_AL*w0*t0*0.04

def test_pressure_drop_grows_with_velocity():
    assert channel_forced(10, 0.005)[1] > channel_forced(5, 0.005)[1] > channel_forced(2, 0.005)[1]

def test_narrow_gap_has_higher_pressure_drop():
    assert array_forced(5, 0.003)["dp"] > array_forced(5, 0.01)["dp"]

if __name__ == "__main__":
    for n, f in list(globals().items()):
        if n.startswith("test_"): f(); print("PASS", n)
