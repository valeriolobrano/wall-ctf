"""CLI entry point for CATI.

Usage:
    python -m cati wall.json
    python -m cati wall.json --roots 100 --poles 20 --coefficients 30
    python -m cati wall.json --sampling-time 1 --periods 10 --output results.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from cati.wall import Wall
from cati.ctf import compute_ctf


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        prog="cati",
        description=(
            "Compute Conduction Transfer Function (CTF) coefficients "
            "for a multilayer wall (Mitalas-Stephenson method)."
        ),
    )
    parser.add_argument(
        "wall_file",
        help="Path to a JSON file describing the wall assembly.",
    )
    parser.add_argument(
        "--roots", type=int, default=50,
        help="Number of roots of B(s) to find (default: 50).",
    )
    parser.add_argument(
        "--poles", type=int, default=None,
        help="Number of poles to use for coefficients (default: same as roots).",
    )
    parser.add_argument(
        "--coefficients", type=int, default=49,
        help="Maximum number of z-domain coefficients (default: 49).",
    )
    parser.add_argument(
        "--sampling-time", type=float, default=1.0,
        help="Sampling period in hours (default: 1.0).",
    )
    parser.add_argument(
        "--periods", type=int, default=10,
        help="Number of 24h periods for transient simulation (default: 10).",
    )
    parser.add_argument(
        "--harmonics", type=int, default=120,
        help="Number of harmonics for Fourier validation (default: 120).",
    )
    parser.add_argument(
        "--t-int", type=float, default=24.0,
        help="Constant internal temperature in °C (default: 24.0).",
    )
    parser.add_argument(
        "--no-mitalas", action="store_true",
        help="Disable the Mitalas instruction.",
    )
    parser.add_argument(
        "--no-fourier", action="store_true",
        help="Skip Fourier validation.",
    )
    parser.add_argument(
        "--output", "-o", type=str, default=None,
        help="Output JSON file for results (default: print to stdout).",
    )

    args = parser.parse_args(argv)

    # Load wall
    wall_path = Path(args.wall_file)
    if not wall_path.exists():
        print(f"Error: file not found: {wall_path}", file=sys.stderr)
        sys.exit(1)

    with open(wall_path) as f:
        wall_data = json.load(f)

    wall = Wall.from_dict(wall_data)

    # Load temperature profile if present
    temperature_profile = None
    if "temperature_profile" in wall_data:
        tp = wall_data["temperature_profile"]
        temperature_profile = np.array(tp, dtype=float)
        if len(temperature_profile) == 24:
            # Add wrap-around: hour 24 = hour 0
            temperature_profile = np.append(temperature_profile, temperature_profile[0])

    # Compute CTF
    result = compute_ctf(
        wall,
        n_roots=args.roots,
        n_poles=args.poles,
        n_coefficients=args.coefficients,
        sampling_time=args.sampling_time,
        n_periods=args.periods,
        n_harmonics=args.harmonics,
        T_int=args.t_int,
        temperature_profile=temperature_profile,
        use_mitalas=not args.no_mitalas,
        validate_fourier=not args.no_fourier,
    )

    # Format output
    output = {
        "wall_name": wall.name,
        "thermal_transmittance_W_m2K": round(result.thermal_transmittance, 6),
        "sampling_time_hours": result.sampling_time,
        "n_poles": result.n_poles,
        "n_coefficients": result.n_coefficients,
        "b_coefficients (1/B, external T)": [
            round(float(v), 10) for v in result.b_coeffs[:result.n_coefficients + 2]
        ],
        "c_coefficients (A/B, internal T)": [
            round(float(v), 10) for v in result.c_coeffs[:result.n_coefficients + 2]
        ],
        "d_coefficients (denominator)": [
            round(float(v), 10) for v in result.d_coeffs[:result.n_coefficients + 1]
        ],
        "poles": [round(float(v), 10) for v in result.poles],
    }
    if result.fourier_error is not None:
        output["fourier_validation_error_percent"] = round(result.fourier_error, 4)

    json_str = json.dumps(output, indent=2)

    if args.output:
        out_path = Path(args.output)
        with open(out_path, "w") as f:
            f.write(json_str + "\n")
        print(f"Results written to {out_path}")
    else:
        print(json_str)


if __name__ == "__main__":
    main()
