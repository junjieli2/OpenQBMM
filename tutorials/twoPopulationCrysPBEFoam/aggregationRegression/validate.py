#!/usr/bin/env python3
import math
import pathlib
import re

case = pathlib.Path(__file__).resolve().parent
end = case / "0.1"

def value(path):
    text = path.read_text()
    match = re.search(r"internalField\s+uniform\s+([^;]+);", text)
    if not match:
        raise RuntimeError(f"Cannot read uniform internalField from {path}")
    return float(match.group(1))

s0 = value(end / "moment.0.singleCrystal")
a0 = value(end / "moment.0.agglomerate")
s3 = value(end / "moment.3.singleCrystal")
a3 = value(end / "moment.3.agglomerate")

Kss, Ksa, Kaa = 2e-7, 1e-7, 5e-8
s, a = 1e6, 0.0
dt = 1e-6
steps = round(0.1/dt)

def rhs(state):
    ss, aa = state
    return (
        -Kss*ss*ss - Ksa*ss*aa,
        0.5*Kss*ss*ss - 0.5*Kaa*aa*aa,
    )

for _ in range(steps):
    k1 = rhs((s, a))
    k2 = rhs((s + 0.5*dt*k1[0], a + 0.5*dt*k1[1]))
    k3 = rhs((s + 0.5*dt*k2[0], a + 0.5*dt*k2[1]))
    k4 = rhs((s + dt*k3[0], a + dt*k3[1]))
    s += dt*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0])/6
    a += dt*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1])/6

rel_s = abs(s0 - s)/max(abs(s), 1)
rel_a = abs(a0 - a)/max(abs(a), 1)
initial_m3 = 4.125e-6
m3_drift = abs((s3 + a3) - initial_m3)/initial_m3

all_values = []
for population in ("singleCrystal", "agglomerate"):
    for order in range(6):
        all_values.append(value(end / f"moment.{order}.{population}"))

assert a0 > 0, "ss collisions did not create the initially empty agglomerate population"
assert rel_s <= 2e-6, f"single m0 differs from RK4: {rel_s}"
assert rel_a <= 2e-6, f"agglomerate m0 differs from RK4: {rel_a}"
assert m3_drift <= 1e-9, f"total third-moment drift: {m3_drift}"
assert all(math.isfinite(v) and v >= -1e-14 for v in all_values)
print(f"PASS aggregationRegression rel_s={rel_s:.3e} rel_a={rel_a:.3e} m3_drift={m3_drift:.3e}")
