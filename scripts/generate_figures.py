"""Generate figures for the README documentation.

Produces a comparison plot of CTF vs Fourier heat flux for a heavy concrete
wall, demonstrating the accuracy of the Z-transform method.

Usage:
    uv run python scripts/generate_figures.py
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use("Agg")

from cati import Wall, Layer, compute_ctf
from cati.ctf import (
    _fourier_validate, _find_roots, _compute_heaviside,
    _compute_denominator, _compute_numerators, _compute_signal_output,
)


def generate_ctf_vs_fourier():
    """Generate CTF vs Fourier comparison for a heavy concrete wall."""

    wall = Wall(
        name="Heavy concrete wall",
        layers=[
            Layer(name="External surface", resistance=0.04),
            Layer(name="External plaster", thickness=0.02, density=1800,
                  specific_heat=1000, conductivity=0.9),
            Layer(name="Concrete block", thickness=0.25, density=2400,
                  specific_heat=1000, conductivity=1.4),
            Layer(name="Internal plaster", thickness=0.02, density=1400,
                  specific_heat=1000, conductivity=0.7),
            Layer(name="Internal surface", resistance=0.13),
        ],
    )

    # Summer sol-air temperature profile for Palermo (typical July day)
    profile = np.array([
        25.0, 24.0, 23.5, 23.0, 22.5, 23.0,
        24.0, 26.0, 28.0, 30.0, 32.0, 34.0,
        36.0, 37.0, 37.5, 37.0, 36.0, 34.0,
        32.0, 30.0, 28.0, 27.0, 26.0, 25.5, 25.0,
    ])
    T_int = 24.0

    # --- Compute CTF coefficients ---
    result = compute_ctf(
        wall, n_roots=30, n_coefficients=20,
        temperature_profile=profile, T_int=T_int,
        n_periods=20, validate_fourier=True,
    )
    nc = result.n_coefficients

    # --- Reconstruct CTF flux for the last period ---
    Tc_s = 3600.0
    tau_max = 24 * 20
    T_ext = np.zeros(tau_max + 1)
    T_ext[0] = profile[0]
    for p in range(20):
        for i in range(1, 25):
            T_ext[p * 24 + i] = profile[i]

    q_ctf_full = _compute_signal_output(
        result.b_coeffs, result.c_coeffs, result.d_coeffs,
        T_ext, T_int, nc, tau_max, 1.0,
    )
    flux_ctf = q_ctf_full[tau_max - 24:tau_max + 1]

    # --- Fourier reference ---
    flux_fourier = _fourier_validate(wall, profile, T_int, 120, 1.0)

    hours = np.arange(25)

    # ---- FIGURE 1: CTF vs Fourier comparison ----
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7.5), height_ratios=[3, 1])
    fig.suptitle(
        f"CTF vs Fourier validation - {wall.name}\n"
        f"U = {wall.thermal_transmittance:.3f} W/(m$^2$K), "
        f"$\\Delta$ = 1 h, {result.n_poles} poles, PME = {result.fourier_error:.2f}%",
        fontsize=12,
    )

    # Top: fluxes
    ax1.plot(hours, flux_fourier[:25], "o-", color="#2563eb", markersize=5,
             linewidth=2, label="Fourier (harmonic analysis)")
    ax1.plot(hours, flux_ctf[:25], "s--", color="#dc2626", markersize=5,
             linewidth=1.5, label="Z-transform (CTF)")
    ax1.set_ylabel("Heat flux $q_i$ [W/m$^2$]", fontsize=11)
    ax1.legend(fontsize=10, loc="upper left")
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, 24)

    # Bottom: input temperature
    ax2.fill_between(hours, profile[:25], alpha=0.15, color="#f59e0b")
    ax2.plot(hours, profile[:25], "^-", color="#f59e0b", markersize=4,
             linewidth=1.5, label="Sol-air temperature")
    ax2.axhline(y=T_int, color="#6b7280", linestyle=":", linewidth=1,
                label=f"Internal air T = {T_int} C")
    ax2.set_xlabel("Time [h]", fontsize=11)
    ax2.set_ylabel("Temperature [C]", fontsize=11)
    ax2.legend(fontsize=9, loc="upper left")
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(0, 24)

    plt.tight_layout()
    fig.savefig("docs/ctf_vs_fourier.png", dpi=150, bbox_inches="tight")
    print("Saved docs/ctf_vs_fourier.png")
    plt.close(fig)

    # ---- FIGURE 2: Wall cross-section schematic (text-based) ----
    print(f"\nResults:")
    print(f"  U-value: {wall.thermal_transmittance:.4f} W/(m2K)")
    print(f"  Poles used: {result.n_poles}")
    print(f"  Coefficients: {nc}")
    print(f"  Fourier PME: {result.fourier_error:.2f}%")
    print(f"  b coefficients: {result.b_coeffs[:nc+1]}")
    print(f"  d coefficients: {result.d_coeffs[:nc+1]}")


if __name__ == "__main__":
    generate_ctf_vs_fourier()
