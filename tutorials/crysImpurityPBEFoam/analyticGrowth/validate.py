#!/usr/bin/env python3
"""Strict, dependency-free regression check for the one-cell tutorial."""

from __future__ import annotations

import math
import re
import sys
from pathlib import Path


CASE = Path(__file__).resolve().parent
FLOAT = r"[+-]?(?:(?:\d+(?:\.\d*)?)|(?:\.\d+))(?:[eE][+-]?\d+)?"


class ValidationError(RuntimeError):
    pass


def internal_scalar(path: Path) -> float:
    if not path.is_file():
        raise ValidationError(f"missing field: {path}")

    contents = path.read_text(encoding="utf-8")
    uniform = re.search(
        rf"\binternalField\s+uniform\s+({FLOAT})\s*;", contents
    )
    if uniform:
        value = float(uniform.group(1))
    else:
        nonuniform = re.search(
            r"\binternalField\s+nonuniform\s+List<scalar>\s+"
            r"(\d+)\s*\((.*?)\)\s*;",
            contents,
            flags=re.DOTALL,
        )
        if not nonuniform:
            raise ValidationError(f"cannot parse scalar internalField: {path}")

        count = int(nonuniform.group(1))
        values = [float(token) for token in re.findall(FLOAT, nonuniform.group(2))]
        if count != 1 or len(values) != 1:
            raise ValidationError(
                f"expected one cell in {path}, found declared={count}, parsed={len(values)}"
            )
        value = values[0]

    if not math.isfinite(value):
        raise ValidationError(f"non-finite internal value in {path}: {value}")
    return value


def latest_time() -> tuple[float, Path]:
    times: list[tuple[float, Path]] = []
    for entry in CASE.iterdir():
        if entry.is_dir() and re.fullmatch(FLOAT, entry.name):
            time_value = float(entry.name)
            if time_value > 0.0:
                times.append((time_value, entry))

    if not times:
        raise ValidationError("no written time directory greater than zero")
    return max(times, key=lambda item: item[0])


def close(label: str, actual: float, expected: float, rtol: float, atol: float) -> None:
    error = abs(actual - expected)
    tolerance = atol + rtol * abs(expected)
    if error > tolerance:
        raise ValidationError(
            f"{label}: actual={actual:.12g}, expected={expected:.12g}, "
            f"error={error:.3g}, tolerance={tolerance:.3g}"
        )


def moment(nodes: tuple[tuple[float, float, float], ...], i: int, j: int) -> float:
    return sum(weight * length**i * impurity**j for weight, length, impurity in nodes)


def validate_source_contracts() -> None:
    """Check the algebraic invariants used by aggregation and breakup."""
    length1, length2 = 1.2, 2.4
    impurity1, impurity2 = 0.2, 0.7
    length_new = (length1**3 + length2**3) ** (1.0 / 3.0)
    impurity_new = impurity1 + impurity2

    for order in ((3, 0), (0, 1)):
        i, j = order
        aggregation_bracket = (
            length_new**i * impurity_new**j
            - length1**i * impurity1**j
            - length2**i * impurity2**j
        )
        close(
            f"aggregation conservation for moment {order}",
            aggregation_bracket,
            0.0,
            rtol=0.0,
            atol=2.0e-12,
        )

    parent_length = 1.7
    for order in ((3, 0), (0, 1)):
        i, j = order
        daughter_order = i + 3 * j
        daughter_moment = 2.0 * (
            parent_length / 2.0 ** (1.0 / 3.0)
        ) ** daughter_order
        if j > 0:
            daughter_moment /= parent_length ** (3 * j)
        breakup_bracket = daughter_moment - parent_length**i
        close(
            f"breakup conservation for moment {order}",
            breakup_bracket,
            0.0,
            rtol=0.0,
            atol=2.0e-12,
        )


def validate() -> None:
    validate_source_contracts()
    time_value, time_dir = latest_time()
    close("final time", time_value, 5.0, rtol=0.0, atol=1.0e-10)

    # These constants mirror populationBalanceProperties and the initial fields.
    temperature0 = 298.15
    solute0 = 150.0
    impurity0 = 2.0
    ki = 0.5
    ks = 0.01
    rho_particle = 2000.0
    shape_factor = 0.523598775598
    size_ref = 1.0e-4
    eta = 0.01
    rho_impurity = 1200.0
    surface_factor = 3.14159265359
    growth_coefficient = 1.0e-6

    initial_nodes = (
        (1.0e9, 1.0e-4, 0.01),
        (5.0e8, 2.0e-4, 0.02),
    )

    csat = (
        -249.369
        + 0.3761 * temperature0
        + 0.0029415 * temperature0**2
    )
    sigma = (solute0 - csat) / csat
    theta = ki * impurity0 / (1.0 + ki * impurity0 + ks * solute0)
    pure_growth = growth_coefficient * sigma
    length_growth = (1.0 - theta) * pure_growth
    m_ref = rho_particle * shape_factor * size_ref**3
    incorporation = eta * rho_impurity * surface_factor / m_ref

    final_nodes = tuple(
        (
            weight,
            length + length_growth * time_value,
            impurity
            + incorporation
            * theta
            / 3.0
            * ((length + length_growth * time_value) ** 3 - length**3),
        )
        for weight, length, impurity in initial_nodes
    )

    orders = ((0, 0), (1, 0), (2, 0), (3, 0), (0, 1), (1, 1))
    actual_moments: dict[tuple[int, int], float] = {}

    for order in orders:
        suffix = f"{order[0]}{order[1]}"
        field_name = f"moment.{suffix}.populationBalance"
        actual = internal_scalar(time_dir / field_name)
        expected = moment(final_nodes, *order)
        initial = moment(initial_nodes, *order)

        if actual < -1.0e-12 * max(1.0, abs(expected)):
            raise ValidationError(f"{field_name} is negative: {actual:.12g}")

        close(field_name, actual, expected, rtol=5.0e-6, atol=1.0e-12)
        actual_moments[order] = actual

        if order != (0, 0) and not actual > initial * (1.0 + 1.0e-7):
            raise ValidationError(
                f"{field_name} did not show the expected positive growth trend: "
                f"initial={initial:.12g}, final={actual:.12g}"
            )

    close(
        "number conservation",
        actual_moments[(0, 0)],
        moment(initial_nodes, 0, 0),
        rtol=1.0e-10,
        atol=1.0e-6,
    )

    temperature = internal_scalar(time_dir / "T")
    solute = internal_scalar(time_dir / "C")
    impurity = internal_scalar(time_dir / "Ci")
    theta_actual = internal_scalar(time_dir / "thetaImpurity")

    if solute < 0.0 or impurity < 0.0:
        raise ValidationError(f"negative dissolved concentration: C={solute}, Ci={impurity}")

    close("T conservation", temperature, temperature0, rtol=1.0e-11, atol=1.0e-10)
    close("C conservation", solute, solute0, rtol=1.0e-11, atol=1.0e-10)
    close("Ci conservation", impurity, impurity0, rtol=1.0e-11, atol=1.0e-12)

    theta_formula = ki * impurity / (1.0 + ki * impurity + ks * solute)
    if not 0.0 < theta_actual < 1.0:
        raise ValidationError(f"thetaImpurity is outside (0,1): {theta_actual}")
    close("competitive adsorption theta", theta_actual, theta_formula, 1.0e-10, 1.0e-12)

    print(
        "PASS: source contracts hold and six bivariate moments match analytic growth; "
        f"theta={theta_actual:.9g}, C={solute:.9g}, Ci={impurity:.9g} at t={time_value:g}"
    )


if __name__ == "__main__":
    try:
        validate()
    except (OSError, ValueError, ValidationError) as error:
        print(f"FAIL: {error}", file=sys.stderr)
        sys.exit(1)
