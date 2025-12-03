"""Command line helper to generate a light curve for an eclipsing binary."""
from __future__ import annotations

import argparse
import csv
import importlib.util
import math
from pathlib import Path

if importlib.util.find_spec("matplotlib"):
    import matplotlib

    # Use a non-interactive backend so the script works in headless environments while still
    # saving the plot alongside the CSV output.
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
else:  # pragma: no cover - optional dependency
    plt = None

from binary_lightcurve import BinarySystem, Star, period_from_parameters


def write_svg_light_curve(path: Path, x_values: list[float], fluxes: list[float], x_label: str) -> None:
    """Write a minimal SVG polyline plot for the light curve when matplotlib is missing."""

    width, height = 800, 450
    margin = 50
    x_min, x_max = min(x_values), max(x_values)
    y_min, y_max = min(fluxes), max(fluxes)
    x_range = x_max - x_min or 1.0
    y_range = y_max - y_min or 1.0

    def sx(x: float) -> float:
        return margin + (x - x_min) / x_range * (width - 2 * margin)

    def sy(y: float) -> float:
        # Flip y because SVG origin is at the top-left corner
        return height - margin - (y - y_min) / y_range * (height - 2 * margin)

    polyline_points = " ".join(f"{sx(x):.2f},{sy(y):.2f}" for x, y in zip(x_values, fluxes))
    axes = (
        f'<line x1="{margin}" y1="{margin}" x2="{margin}" y2="{height-margin}" stroke="black" stroke-width="1" />'
        f'<line x1="{margin}" y1="{height-margin}" x2="{width-margin}" y2="{height-margin}" stroke="black" stroke-width="1" />'
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        f"""
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}">
  <rect width="100%" height="100%" fill="white" />
  {axes}
  <polyline fill="none" stroke="steelblue" stroke-width="2" points="{polyline_points}" />
  <text x="{width / 2}" y="{height - margin}" font-size="14" text-anchor="middle" dy="28">{x_label}</text>
  <text x="{margin}" y="{margin}" font-size="14" text-anchor="start" dy="-10">Relative flux</text>
</svg>
"""
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Eclipsing binary light curve simulator")
    parser.add_argument("--mass1", type=float, default=1.989e30, help="Primary mass in kg")
    parser.add_argument("--mass2", type=float, default=1.5e30, help="Secondary mass in kg")
    parser.add_argument("--radius1", type=float, default=6.957e8, help="Primary radius in m")
    parser.add_argument("--radius2", type=float, default=5e8, help="Secondary radius in m")
    parser.add_argument("--temp1", type=float, default=6000, help="Primary temperature in K")
    parser.add_argument("--temp2", type=float, default=5000, help="Secondary temperature in K")
    parser.add_argument("--semi-major-axis", type=float, default=2.5e10, help="Semi-major axis in m")
    parser.add_argument("--eccentricity", type=float, default=0.1, help="Orbital eccentricity (0-1)")
    parser.add_argument("--inclination", type=float, default=math.radians(87), help="Inclination in radians")
    parser.add_argument("--phases", type=int, default=400, help="Number of phase samples")
    parser.add_argument("--output", type=Path, default=Path("light_curve.csv"), help="CSV file for output")
    parser.add_argument(
        "--figure",
        type=Path,
        default=None,
        help="Optional path for saving the light-curve plot (PNG). Defaults to output stem with .png",
    )
    parser.add_argument(
        "--x-axis",
        choices=["phase", "day", "minute"],
        default="phase",
        help="Use orbital phase or convert to time in days/minutes for plots and CSV",
    )
    return parser.parse_args()


def compute_light_curve(
    *,
    mass1: float,
    mass2: float,
    radius1: float,
    radius2: float,
    temp1: float,
    temp2: float,
    semi_major_axis: float,
    eccentricity: float,
    inclination: float,
    phases: int,
) -> tuple[list[float], list[float], list[float], list[float], float]:
    """Create a binary system and return phases, fluxes, RVs, and the orbital period."""

    primary = Star(mass1, radius1, temp1)
    secondary = Star(mass2, radius2, temp2)

    period = period_from_parameters(semi_major_axis, mass1, mass2)
    system = BinarySystem(
        primary=primary,
        secondary=secondary,
        semi_major_axis=semi_major_axis,
        eccentricity=eccentricity,
        inclination=inclination,
        period=period,
    )

    phases_values, fluxes, rv_primary, rv_secondary = system.light_curve(phases)
    return phases_values, fluxes, rv_primary, rv_secondary, period


def main() -> None:
    args = parse_args()

    phases, fluxes, rv1, rv2, period = compute_light_curve(
        mass1=args.mass1,
        mass2=args.mass2,
        radius1=args.radius1,
        radius2=args.radius2,
        temp1=args.temp1,
        temp2=args.temp2,
        semi_major_axis=args.semi_major_axis,
        eccentricity=args.eccentricity,
        inclination=args.inclination,
        phases=args.phases,
    )
    time_seconds = [phase * period for phase in phases]
    if args.x_axis == "day":
        x_values = [t / 86400 for t in time_seconds]
        x_label = "Time (days)"
        time_header = "time_days"
    elif args.x_axis == "minute":
        x_values = [t / 60 for t in time_seconds]
        x_label = "Time (minutes)"
        time_header = "time_minutes"
    else:
        x_values = phases
        x_label = "Orbital phase"
        time_header = None

    with args.output.open("w", newline="") as f:
        writer = csv.writer(f)
        header = ["phase"]
        if time_header:
            header.append(time_header)
        header += ["relative_flux", "rv_primary", "rv_secondary"]
        writer.writerow(header)

        rows = []
        for idx in range(len(phases)):
            row = [phases[idx]]
            if time_header:
                row.append(x_values[idx])
            row += [fluxes[idx], rv1[idx], rv2[idx]]
            rows.append(row)
        writer.writerows(rows)

    print(f"Saved {len(phases)} samples to {args.output.resolve()}")
    print(f"Orbital period: {period / 86400:.2f} days")

    fig_path = args.figure if args.figure is not None else args.output.with_suffix(".png")
    if plt is not None:
        fig, (ax_flux, ax_rv) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

        ax_flux.plot(x_values, fluxes, color="tab:blue", linewidth=1.5)
        ax_flux.set_ylabel("Relative flux")
        ax_flux.set_title("Eclipsing binary light curve")
        ax_flux.grid(True, alpha=0.3)

        ax_rv.plot(x_values, rv1, color="tab:red", linewidth=1.2, label="Primary RV")
        ax_rv.plot(x_values, rv2, color="tab:green", linewidth=1.2, label="Secondary RV")
        ax_rv.set_xlabel(x_label)
        ax_rv.set_ylabel("Line-of-sight velocity (m/s)")
        ax_rv.set_title("Radial velocity curves")
        ax_rv.legend()
        ax_rv.grid(True, alpha=0.3)

        fig.tight_layout()
        fig.savefig(fig_path)
        print(f"Saved plot to {fig_path.resolve()}")
    else:
        svg_path = fig_path if fig_path.suffix.lower() == ".svg" else fig_path.with_suffix(".svg")
        write_svg_light_curve(svg_path, x_values, fluxes, x_label)
        print(
            "matplotlib is not installed; saved a fallback SVG plot to " f"{svg_path.resolve()}"
        )


if __name__ == "__main__":
    main()
