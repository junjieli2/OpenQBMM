#!/usr/bin/env python3
import math
import pathlib
import re

case = pathlib.Path(__file__).resolve().parent
end = case / "0.2"
log = (case / "log.twoPopulationCrysPBEFoam").read_text()

assert "Temperature mode: prescribed" in log
assert log.count("Selecting aggregationKernel linnikovAggregation") == 3
assert log.count("K1D3 =") == 3

def value(path):
    text = path.read_text()
    match = re.search(r"internalField\s+uniform\s+([^;]+);", text)
    if not match:
        raise RuntimeError(path)
    return float(match.group(1))

rhop = 2666.0
shape_factor = 0.5235987755982988
initial_m3 = 4.125e-4 + 8.25e-5
initial_total = 150.0 + rhop*shape_factor*initial_m3
final_m3 = sum(
    value(end / f"moment.3.{population}")
    for population in ("singleCrystal", "agglomerate")
)
final_c = value(end / "C")
final_total = final_c + rhop*shape_factor*final_m3
conservation_error = abs(final_total - initial_total)/abs(initial_total)

values = [
    value(end / f"moment.{order}.{population}")
    for population in ("singleCrystal", "agglomerate")
    for order in range(6)
]
gs = value(end / "growthRate.singleCrystal")
ga = value(end / "growthRate.agglomerate")

assert conservation_error <= 1e-6, conservation_error
assert final_c < 150.0, "cooling/growth did not consume solute"
assert gs > ga > 0.0, (gs, ga)
assert all(math.isfinite(v) and v >= -1e-14 for v in values)
print(
    "PASS K2SO4Batch "
    f"conservation_error={conservation_error:.3e} Gs={gs:.3e} Ga={ga:.3e}"
)
