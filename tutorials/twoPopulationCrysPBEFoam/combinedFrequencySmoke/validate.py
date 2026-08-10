#!/usr/bin/env python3
import math
import pathlib
import re

case = pathlib.Path(__file__).resolve().parent
end = case / "0.02"
log = (case / "log.twoPopulationCrysPBEFoam").read_text()

for selected in (
    "Selecting crystalCollisionFrequency sum",
    "Selecting crystalCollisionFrequency turbulent",
    "Selecting crystalCollisionFrequency shear",
    "Selecting crystalAggregationEfficiency Ilievski",
):
    assert selected in log, f"runtime-selection path not exercised: {selected}"

def value(path):
    text = path.read_text()
    match = re.search(r"internalField\s+uniform\s+([^;]+);", text)
    if not match:
        raise RuntimeError(path)
    return float(match.group(1))

initial_m3 = 4.125e-4 + 8.25e-5
final_m3 = sum(value(end / f"moment.3.{p}") for p in ("singleCrystal", "agglomerate"))
drift = abs(final_m3 - initial_m3)/initial_m3
values = [
    value(end / f"moment.{k}.{p}")
    for p in ("singleCrystal", "agglomerate")
    for k in range(6)
]
assert drift <= 1e-9, f"total third-moment drift: {drift}"
assert all(math.isfinite(v) and v >= -1e-14 for v in values)
print(f"PASS combinedFrequencySmoke m3_drift={drift:.3e}")
