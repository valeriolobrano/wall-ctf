"""Conduction Transfer Function (CTF) coefficient computation.

This module implements the Transfer Function Method (TFM) for computing
the Z-transform coefficients of the heat conduction transfer functions
through multilayer wall assemblies. The method is based on the
Mitalas-Stephenson approach, widely used in TRNSYS, DOE-2, and ASHRAE
procedures for dynamic building energy simulation.

Algorithm overview:

1. **Transmission matrix**: For each homogeneous layer, the thermal
   transmission matrix in the Laplace domain is (Ref. [1], Eq. 8):

       | a  b |     | cosh(L*sqrt(s/alpha))        sinh(L*sqrt(s/alpha)) / (k*sqrt(s/alpha)) |
       | c  d |  =  | k*sqrt(s/alpha)*sinh(...)     cosh(L*sqrt(s/alpha))                     |

   The overall wall matrix is the product of all layer matrices including
   the surface resistance matrices (Ref. [1], Eq. 9-10).

2. **Root finding**: The poles of the transfer function are the values
   of s that make DEN(s) = B(s) = 0, where B is the (0,1) element of
   the overall matrix. Using the substitution sqrt(s) = j*xi, the search
   is reduced to the real axis (Ref. [1], Sec. 3; Ref. [2], Appendix A).

3. **Heaviside expansion**: The ramp response is expanded as (Ref. [1], Eq. 4):
       O(s) = 1/s^2 * NUM(s)/DEN(s) = C0/s^2 + C1/s + sum resn/(s - sn)

   where resn = NUM(sn) / (sn^2 * DEN'(sn))  (Ref. [1], Eq. 26; Ref. [2], Eq. A.4).

4. **Procedure I** (optimal pole/residue selection): Residues are sorted by
   magnitude and only significant ones (|resn| > 1e-10) are retained.
   This dramatically improves reliability for massive walls (Ref. [2], Sec. 5).

5. **Z-domain coefficients**: The denominator is den(z) = prod(1 - e^{sn*Delta}*z^{-1})
   (Ref. [1], Eq. 6) and numerators are computed from the ramp response
   convolution (Ref. [1], Eq. 5; Ref. [2], Eq. A.5-A.8).

6. **Fourier validation**: An independent solution using the complex thermal
   quadrupole (harmonic analysis) provides a reference for error assessment
   (Ref. [1], Sec. 4; Ref. [2], Sec. 3).

References:
    [1] G. Beccali, M. Cellura, V. Lo Brano, A. Orioli,
        "Single thermal zone balance solved by Transfer Function Method",
        Energy and Buildings, 37 (2005) 1268-1277.
        DOI: 10.1016/j.enbuild.2005.02.010

    [2] G. Beccali, M. Cellura, V. Lo Brano, A. Orioli,
        "Is the transfer function method reliable in a European building
        context? A theoretical analysis and a case study in the south
        of Italy", Applied Thermal Engineering, 25 (2005) 341-357.
        DOI: 10.1016/j.applthermaleng.2004.06.010

    [3] Mitalas, G.P., Stephenson, D.G. (1967). Room Thermal Response
        Factors. ASHRAE Transactions, 73(1).

Copyright (c) 2005-2026 Valerio Lo Brano, Università degli Studi di Palermo.
License: CC BY-NC 4.0 (non-commercial use).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from cati.wall import Wall


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------

@dataclass
class CTFResult:
    """Result of a CTF coefficient computation.

    The heat flux at the internal surface is computed as:
        q_i(n) = sum_j b_j * T_e(n-j) - sum_{j>=1} d_j * q_i(n-j)
                 - T_i * sum_j c_j

    where:
        b_j: coefficients driven by external temperature (numerator of 1/B)
        c_j: coefficients driven by internal temperature (numerator of A/B)
        d_j: denominator (common to both transfer functions)

    Attributes:
        b_coeffs: Numerator coefficients for external temperature [1/B].
        c_coeffs: Numerator coefficients for internal temperature [A/B].
        d_coeffs: Common denominator coefficients.
        n_poles: Number of poles used.
        n_coefficients: Number of z-coefficients.
        sampling_time: Sampling period in hours.
        thermal_transmittance: Steady-state U-value [W/(m²·K)].
        poles: Array of poles s_n in the Laplace domain.
        fourier_error: Percentage error vs Fourier solution (if computed).
    """

    b_coeffs: np.ndarray = field(default_factory=lambda: np.array([]))
    c_coeffs: np.ndarray = field(default_factory=lambda: np.array([]))
    d_coeffs: np.ndarray = field(default_factory=lambda: np.array([]))
    n_poles: int = 0
    n_coefficients: int = 0
    sampling_time: float = 1.0
    thermal_transmittance: float = 0.0
    poles: np.ndarray = field(default_factory=lambda: np.array([]))
    fourier_error: float | None = None


# ---------------------------------------------------------------------------
# 2×2 matrix helpers (pure numpy)
# ---------------------------------------------------------------------------

def _mat_mul(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Multiply two 2×2 matrices."""
    return a @ b


def _mat_zeros() -> np.ndarray:
    return np.zeros((2, 2))


def _mat_identity() -> np.ndarray:
    return np.eye(2)


# ---------------------------------------------------------------------------
# Transmission matrix for a single layer at Laplace variable s = -xi²
# ---------------------------------------------------------------------------

def _layer_matrix(xi: float, thickness: float, conductivity: float,
                  diffusivity: float, resistance: float) -> np.ndarray:
    """Compute the transmission matrix of a single layer evaluated at s = -xi^2.

    Ref. [1] Eq. (8): For a homogeneous isotropic layer of thickness L,
    conductivity k, and thermal diffusivity alpha = k/(rho*Cp), the Laplace-
    domain transmission matrix relates temperatures and heat fluxes on the
    two sides. Evaluated at s = -xi^2 (negative real axis), the hyperbolic
    functions become trigonometric:

        cosh(L*sqrt(s/alpha)) -> cos(xi*L/sqrt(alpha))
        sinh(L*sqrt(s/alpha))/(k*sqrt(s/alpha)) -> sin(xi*L/sqrt(alpha))/(k*xi/sqrt(alpha))

    For air gaps (thickness=0), the matrix is simply [[1, R], [0, 1]].
    """
    M = np.eye(2)
    if thickness == 0.0:
        M[0, 1] = resistance
        return M

    if xi == 0.0:
        M[0, 1] = thickness / conductivity
        return M

    sqrt_alpha = math.sqrt(diffusivity)
    u = xi * thickness / sqrt_alpha
    cos_u = math.cos(u)
    sin_u = math.sin(u)
    gamma = xi * conductivity / sqrt_alpha  # k * xi / sqrt(alpha)

    M[0, 0] = cos_u
    M[0, 1] = sin_u / gamma
    M[1, 0] = -gamma * sin_u
    M[1, 1] = cos_u
    return M


def _layer_matrix_deriv_s0(thickness: float, conductivity: float,
                           diffusivity: float) -> np.ndarray:
    """Derivative of layer transmission matrix w.r.t. s, evaluated at s = 0.

    Needed for the C1 coefficient of the Heaviside expansion (Ref. [1] Eq. 4).
    Obtained from the Taylor expansion of the matrix elements around s = 0.
    For a conductive layer:

        dM/ds|_{s=0} = | L^2/(2*alpha)       L^3/(6*k*alpha) |
                       | k*L/alpha            L^2/(2*alpha)   |
    """
    if thickness == 0.0:
        return _mat_zeros()
    L = thickness
    alpha = diffusivity
    k = conductivity
    d = np.zeros((2, 2))
    d[0, 0] = 0.5 * L ** 2 / alpha
    d[0, 1] = L ** 3 / (6.0 * k * alpha)
    d[1, 0] = k * L / alpha
    d[1, 1] = d[0, 0]
    return d


def _layer_matrix_deriv(xi: float, thickness: float, conductivity: float,
                        diffusivity: float) -> np.ndarray:
    """Derivative of layer transmission matrix w.r.t. s, at s = -xi^2 (xi > 0).

    Needed for the residue computation (Ref. [2] Eq. A.4):
        resn = NUM(sn) / (sn^2 * d/ds{DEN(s)}|_{s=sn})

    The derivative dM/ds is obtained analytically from dM/dxi via the
    chain rule: since s = -xi^2, ds = -2*xi*dxi, hence dM/ds = dM/dxi / (-2*xi).
    """
    if thickness == 0.0:
        return _mat_zeros()

    L = thickness
    alpha = diffusivity
    k = conductivity
    sqrt_alpha = math.sqrt(alpha)
    u = xi * L / sqrt_alpha
    cos_u = math.cos(u)
    sin_u = math.sin(u)

    d = np.zeros((2, 2))
    # dA/ds = 0.5 * L/sqrt(alpha) * sin(u) / xi
    d[0, 0] = 0.5 * L / sqrt_alpha * sin_u / xi
    # dB/ds = sqrt(alpha)/(2k) * [sin(u)/xi³ - (L/sqrt(alpha))*cos(u)/xi²]
    d[0, 1] = 0.5 * sqrt_alpha / k * (sin_u / xi ** 3 - (L / sqrt_alpha) * cos_u / xi ** 2)
    # dC/ds = k/(2*sqrt(alpha)) * [sin(u)/xi + (L/sqrt(alpha))*cos(u)]
    d[1, 0] = 0.5 * k / sqrt_alpha * (sin_u / xi + (L / sqrt_alpha) * cos_u)
    # dD/ds = dA/ds (symmetric matrix)
    d[1, 1] = d[0, 0]
    return d


# ---------------------------------------------------------------------------
# Overall transmission matrix and its derivative
# ---------------------------------------------------------------------------

def _compute_overall_matrix(xi: float, wall: Wall) -> np.ndarray:
    """Compute the overall transmission matrix M = He * M2 * ... * M_{N-1} * Hi."""
    # Start with external surface resistance
    R_ext = wall.external_surface_resistance
    M = np.array([[1.0, R_ext], [0.0, 1.0]])

    # Multiply by each material layer
    for layer in wall.material_layers:
        M_layer = _layer_matrix(xi, layer.thickness, layer.conductivity,
                                layer.thermal_diffusivity, layer.resistance)
        M = M @ M_layer

    # Multiply by internal surface resistance
    R_int = wall.internal_surface_resistance
    Hi = np.array([[1.0, R_int], [0.0, 1.0]])
    M = M @ Hi
    return M


def _compute_overall_matrix_deriv(xi: float, wall: Wall,
                                  layer_matrices: list[np.ndarray]) -> np.ndarray:
    """Compute dM_total/ds using the matrix product rule.

    The overall matrix is M = He * M_2 * ... * M_{N-1} * Hi (Ref. [1] Eq. 10).
    Its derivative w.r.t. s is (Leibniz rule for matrix products):

        dM/ds = sum_i [He * M_2 * ... * dM_i/ds * ... * M_{N-1} * Hi]

    i.e. for each layer i, replace M_i with dM_i/ds and keep all others.
    The (0,1) element of this derivative gives DEN'(s) needed for the
    residue formula (Ref. [2] Eq. A.4).
    """
    mat_layers = wall.material_layers
    n_mat = len(mat_layers)

    R_ext = wall.external_surface_resistance
    He = np.array([[1.0, R_ext], [0.0, 1.0]])
    R_int = wall.internal_surface_resistance
    Hi = np.array([[1.0, R_int], [0.0, 1.0]])

    Z = _mat_zeros()

    for i in range(n_mat):
        layer = mat_layers[i]
        if xi == 0.0:
            dM_i = _layer_matrix_deriv_s0(layer.thickness, layer.conductivity,
                                          layer.thermal_diffusivity)
        else:
            dM_i = _layer_matrix_deriv(xi, layer.thickness, layer.conductivity,
                                       layer.thermal_diffusivity)

        # Build product: He * M_0 * ... * M_{i-1} * dM_i * M_{i+1} * ... * M_{n-1} * Hi
        Q = He.copy()
        for j in range(n_mat):
            if j == i:
                Q = Q @ dM_i
            else:
                Q = Q @ layer_matrices[j]
        Q = Q @ Hi
        Z = Z + Q

    return Z


# ---------------------------------------------------------------------------
# Root finding: find zeros of B(s) along the negative real s-axis
# ---------------------------------------------------------------------------

def _find_roots(wall: Wall, n_roots: int,
                delta_init: float | None = None,
                limit_init: float = 1e-14) -> tuple[np.ndarray, list, list]:
    """Find the first n_roots zeros of B(s) on the negative real axis.

    B(s) is the (0,1) element of the overall transmission matrix.
    We parametrize s = -xi² and scan xi > 0 looking for sign changes.

    Returns:
        poles: array of s_n = -xi_n² values
        mat_at_roots: list of 2×2 matrices M(s_n)
        layer_mats_at_roots: list of lists of individual layer matrices at each root
    """
    if delta_init is None:
        # Adaptive step: estimate root spacing from individual layers.
        # For a single slab, roots are spaced at Δξ = π·√α/L.
        # For a multilayer wall, roots from different layers interleave,
        # so the effective spacing is roughly min_spacing / n_conductive_layers.
        min_spacing = float('inf')
        n_conductive = 0
        for layer in wall.material_layers:
            if layer.thickness > 0 and layer.thermal_diffusivity > 0:
                spacing = math.pi * math.sqrt(layer.thermal_diffusivity) / layer.thickness
                min_spacing = min(min_spacing, spacing)
                n_conductive += 1
        if min_spacing == float('inf') or n_conductive == 0:
            delta_init = 0.01
        else:
            # At least 5 steps per root spacing, divided by number of layers
            delta_init = min_spacing / (5.0 * max(n_conductive, 1))

    poles = np.zeros(n_roots)
    mat_at_roots = []
    layer_mats_at_roots = []

    xi = 0.0
    delta = delta_init
    limit = limit_init
    b_prev = 0.0
    roots_found = 0
    attempt = 0

    while roots_found < n_roots:
        # Compute overall matrix at current xi
        M = _compute_overall_matrix(xi, wall)
        b_curr = M[0, 1]

        if b_curr == b_prev and attempt > 0:
            limit *= 10.0

        if attempt > 0:
            if b_curr * b_prev < 0:
                delta = -delta / 3.0

        if abs(b_curr) < limit or abs(delta) < 1e-30:
            # Root found
            # Store individual layer matrices for derivative computation
            mat_layers = wall.material_layers
            lmats = []
            for layer in mat_layers:
                lm = _layer_matrix(xi, layer.thickness, layer.conductivity,
                                   layer.thermal_diffusivity, layer.resistance)
                lmats.append(lm)

            poles[roots_found] = -xi ** 2
            mat_at_roots.append(M.copy())
            layer_mats_at_roots.append(lmats)
            roots_found += 1

            xi += delta_init
            delta = delta_init
            limit = limit_init
            attempt = 0
            b_prev = 0.0
        else:
            b_prev = b_curr
            attempt += 1
            xi += delta

    return poles, mat_at_roots, layer_mats_at_roots


# ---------------------------------------------------------------------------
# Heaviside partial-fraction expansion
# ---------------------------------------------------------------------------

def _compute_heaviside(wall: Wall, poles: np.ndarray,
                       mat_at_roots: list, layer_mats_at_roots: list,
                       use_mitalas: bool = True) -> tuple:
    """Compute the Heaviside partial-fraction expansion coefficients.

    The heat flux at the internal surface involves two sub-transfer functions
    (Ref. [1] Eq. 11):   Q_i(z) = [1/B]*T_e - [A/B]*T_i

    For each, the ramp response in Laplace domain is (Ref. [1] Eq. 4;
    Ref. [2] Eq. A.3):

        O(s) = 1/s^2 * NUM(s)/DEN(s) = C0/s^2 + C1/s + sum_n resn/(s - sn)

    where:
        C0 = NUM(0)/DEN(0)  (= 1/B(0) = U for the 1/B case)
        C1 = [NUM'(0)*DEN(0) - NUM(0)*DEN'(0)] / DEN(0)^2
        resn = NUM(sn) / (sn^2 * DEN'(sn))   (Ref. [2] Eq. A.4)

    The Mitalas instruction sets C1 = -sum(resn) to ensure the ramp
    response starts exactly at zero for t=0 (Ref. [3]).
    """
    n_roots = len(poles)

    # --- Compute at s=0 (xi=0) ---
    M0 = _compute_overall_matrix(0.0, wall)
    A0_val = M0[0, 0]
    B0_val = M0[0, 1]

    # Layer matrices at xi=0
    mat_layers = wall.material_layers
    lmats_0 = []
    for layer in mat_layers:
        lm = _layer_matrix(0.0, layer.thickness, layer.conductivity,
                           layer.thermal_diffusivity, layer.resistance)
        lmats_0.append(lm)

    # Derivative at s=0
    dM0 = _compute_overall_matrix_deriv(0.0, wall, lmats_0)
    dA0 = dM0[0, 0]
    dB0 = dM0[0, 1]

    C0_1B = 1.0 / B0_val
    C0_AB = A0_val / B0_val

    C1_1B = -dB0 / (B0_val ** 2)
    C1_AB = (B0_val * dA0 - A0_val * dB0) / (B0_val ** 2)

    # --- Compute residues at each pole ---
    residues_1B = np.zeros(n_roots)
    residues_AB = np.zeros(n_roots)

    for k in range(n_roots):
        s_k = poles[k]
        xi_k = math.sqrt(-s_k)
        M_k = mat_at_roots[k]
        A_k = M_k[0, 0]

        # Derivative of overall matrix at this root
        dM_k = _compute_overall_matrix_deriv(xi_k, wall, layer_mats_at_roots[k])
        dB_k = dM_k[0, 1]

        # Residue: dn = 1/(s_n² · dB/ds|_{s_n}) for 1/B
        #          dn = A(s_n)/(s_n² · dB/ds|_{s_n}) for A/B
        denom = s_k ** 2 * dB_k
        residues_1B[k] = 1.0 / denom
        residues_AB[k] = A_k / denom

    # --- Mitalas instruction ---
    if use_mitalas:
        C1_1B = -np.sum(residues_1B)
        C1_AB = -np.sum(residues_AB)

    return C0_1B, C0_AB, C1_1B, C1_AB, residues_1B, residues_AB


# ---------------------------------------------------------------------------
# Denominator coefficients in the z-domain
# ---------------------------------------------------------------------------

def _compute_denominator(Tc: float, poles: np.ndarray,
                         n_poles: int, n_coefficients: int) -> np.ndarray:
    """Compute the common denominator of the Z-transfer functions.

    Ref. [1] Eq. (6):  den(z) = prod_{i=1}^{Np} (1 - e^{sn*Delta} * z^{-1})

    where sn are the Laplace-domain poles and Delta is the sampling period
    in seconds. The z-domain poles are z_i = exp(sn * Delta), all inside the
    unit circle since sn < 0.

    The polynomial is expanded iteratively by multiplying one factor at a time,
    giving O(Np^2) complexity (vs exponential for the recursive elementary
    symmetric polynomial approach used in the original VB code).
    """
    z_vals = np.exp(Tc * poles[:n_poles])

    # Build polynomial by multiplying (1 - z_i * z^{-1}) factors iteratively
    # d[k] is the coefficient of z^{-k}
    n_coeff = min(n_coefficients, n_poles)
    d = np.zeros(n_coefficients + 1)
    d[0] = 1.0

    for i in range(n_poles):
        # Multiply current polynomial by (1 - z_i * z^{-1})
        # Process from high to low to avoid overwriting
        for k in range(min(i + 1, n_coeff), 0, -1):
            d[k] = d[k] - z_vals[i] * d[k - 1]

    return d


# ---------------------------------------------------------------------------
# Numerator coefficients (ramp interpolation method)
# ---------------------------------------------------------------------------

def _compute_numerators(Tc: float, poles: np.ndarray, n_roots: int,
                        n_coefficients: int,
                        C0_1B: float, C0_AB: float,
                        C1_1B: float, C1_AB: float,
                        residues_1B: np.ndarray, residues_AB: np.ndarray,
                        denom: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Compute the numerator coefficients for [1/B] and [A/B] sub-functions.

    Ref. [1] Eq. (5); Ref. [2] Eq. (A.5)-(A.8).

    The numerators are computed from the sampled ramp response O(n) and the
    denominator using a convolution procedure. The ramp response at t = n*Tc is:

        O(n) = C0*n*Tc + C1 + sum_i resn_i * exp(sn_i * n * Tc)

    The M terms are the second differences of the denominator divided by Tc
    (related to the ramp interpolation kernel). The numerator coefficients
    are obtained as the convolution of M and O.

    Returns:
        b_coeffs: numerator of 1/B (driven by external temperature T_e)
        c_coeffs: numerator of A/B (driven by internal temperature T_i)
    """
    # --- Compute O terms (sampled ramp response) ---
    O_1B = np.zeros(n_roots + 2)
    O_AB = np.zeros(n_roots + 2)

    for n in range(1, n_roots + 2):
        t = n * Tc
        sum_1B = np.sum(residues_1B[:n_roots] * np.exp(t * poles[:n_roots]))
        sum_AB = np.sum(residues_AB[:n_roots] * np.exp(t * poles[:n_roots]))
        O_1B[n] = sum_1B + C1_1B + C0_1B * t
        O_AB[n] = sum_AB + C1_AB + C0_AB * t
    # O[0] = 0 (already initialized)

    # --- Compute M terms (second differences of denominator / Tc) ---
    nc = n_coefficients
    M_arr = np.zeros(nc + 2)
    M_minus1 = denom[0] / Tc

    M_arr[0] = (denom[1] - 2.0 * denom[0]) / Tc
    for n in range(1, nc):
        M_arr[n] = (denom[n + 1] - 2.0 * denom[n] + denom[n - 1]) / Tc
    M_arr[nc] = (-2.0 * denom[nc] + denom[nc - 1]) / Tc
    if nc + 1 < len(M_arr):
        M_arr[nc + 1] = denom[nc] / Tc

    # --- Compute numerator coefficients via convolution ---
    c_coeffs = np.zeros(nc + 2)
    d_coeffs = np.zeros(nc + 2)

    for x in range(nc + 2):
        sum_1B = M_minus1 * O_1B[min(x + 1, len(O_1B) - 1)]
        sum_AB = M_minus1 * O_AB[min(x + 1, len(O_AB) - 1)]

        for y in range(x + 1):
            o_idx = x - y
            if o_idx >= len(O_1B):
                o_val_1B = 0.0
                o_val_AB = 0.0
            else:
                o_val_1B = O_1B[o_idx]
                o_val_AB = O_AB[o_idx]

            if y >= len(M_arr):
                m_val = 0.0
            else:
                m_val = M_arr[y]

            sum_1B += m_val * o_val_1B
            sum_AB += m_val * o_val_AB

        c_coeffs[x] = sum_1B
        d_coeffs[x] = sum_AB

    return c_coeffs, d_coeffs


# ---------------------------------------------------------------------------
# Signal output computation (z-domain recurrence)
# ---------------------------------------------------------------------------

def _compute_signal_output(b_coeffs: np.ndarray, c_coeffs: np.ndarray,
                           d_coeffs: np.ndarray, T_ext: np.ndarray,
                           T_int: float, n_coefficients: int,
                           tau_max: int, sampling_time: float) -> np.ndarray:
    """Compute heat flux at the internal surface using the CTF recurrence.

    Ref. [1] Eq. (29); Ref. [2] Eq. (A.8):

        T_{x,i}(n*Delta) = sum_j num_j * I_i[(n-j)*Delta]
                         - sum_j den_j * T_{x,i}[(n-j)*Delta]

    In our notation for the wall heat flux (Ref. [1] Eq. 11):

        q_i(N) = sum_j b_j * T_e(N-j) - sum_{j>=1} d_j * q_i(N-j) - T_i * sum(c_j)

    Args:
        b_coeffs: numerator coefficients for T_e (1/B)
        c_coeffs: numerator coefficients for T_i (A/B)
        d_coeffs: denominator coefficients
        T_ext: external (sol-air) temperature array, indexed 0..tau_max
        T_int: constant internal temperature
        n_coefficients: number of z-coefficients to use
        tau_max: total number of time steps to simulate
        sampling_time: sampling period in hours
    """
    nc = n_coefficients
    q = np.zeros(tau_max + 1)
    sum_c = np.sum(c_coeffs[:nc + 1])
    n_per_day = int(24.0 / sampling_time)

    def get_T_ext(idx: int) -> float:
        """Get external temperature, wrapping to first period for negative indices."""
        if 0 <= idx <= tau_max:
            return T_ext[idx]
        # Wrap to first period
        wrapped = idx % n_per_day
        if wrapped <= 0:
            wrapped += n_per_day
        return T_ext[wrapped]

    for N in range(1, tau_max + 1):
        temp1 = 0.0
        for j in range(nc + 1):
            temp1 += b_coeffs[j] * get_T_ext(N - j)

        temp2 = 0.0
        for j in range(1, nc + 1):
            if N - j > 0:
                temp2 += d_coeffs[j] * q[N - j]

        q[N] = temp1 - temp2 - T_int * sum_c

    return q


# ---------------------------------------------------------------------------
# Fourier (complex quadrupole) validation
# ---------------------------------------------------------------------------

def _fourier_validate(wall: Wall, temperature_profile: np.ndarray,
                      T_int: float, n_harmonics: int,
                      sampling_time: float) -> np.ndarray:
    """Compute heat flux using the Fourier (complex quadrupole) method.

    This is the analytical solution used to validate the CTF results.

    Args:
        wall: Wall definition.
        temperature_profile: 24-point hourly temperature profile.
        T_int: Internal temperature [°C].
        n_harmonics: Number of Fourier harmonics.
        sampling_time: Sampling period [hours].

    Returns:
        flux_fourier: Array of heat flux values [W/m²] at each sampling time.
    """
    # --- Interpolate temperature profile to 240 points ---
    n_points = 240
    input_signal = _interpolate_profile(temperature_profile, n_points)

    # --- Compute Fourier coefficients via vectorized DFT ---
    N = n_points
    sig = input_signal[1:N + 1]  # N samples
    A0 = np.mean(sig)

    # Vectorized: m = 1..n_harmonics, n = 1..N
    m_arr = np.arange(1, n_harmonics + 1)
    n_arr = np.arange(1, N + 1)
    angles = np.outer(m_arr, n_arr) * (2.0 * math.pi / N)  # (M, N)
    Am = np.zeros(n_harmonics + 1)
    Bm = np.zeros(n_harmonics + 1)
    Am[1:] = (2.0 / N) * (np.cos(angles) @ sig)
    Bm[1:] = (2.0 / N) * (np.sin(angles) @ sig)

    amplitude = np.sqrt(Am ** 2 + Bm ** 2)
    phase = np.arctan2(Am, Bm)
    phase[phase < 0] += 2.0 * math.pi

    # --- Compute output for each harmonic using complex quadrupole ---
    # Vectorized: compute all harmonics at once using numpy complex matrices.
    h_ext = 1.0 / wall.external_surface_resistance
    h_int = 1.0 / wall.internal_surface_resistance
    R_total = wall.total_thermal_resistance
    U = 1.0 / R_total

    M_harmonics = n_harmonics
    omegas = (2.0 * math.pi / (24.0 * 3600.0)) * np.arange(1, M_harmonics + 1)

    # Build complex transfer matrices for all harmonics at once.
    # Each harmonic m has a 2x2 complex matrix. We batch as (M, 2, 2).
    M_ext = np.zeros((M_harmonics, 2, 2), dtype=complex)
    M_ext[:, 0, 0] = 1.0
    M_ext[:, 0, 1] = 1.0 / h_ext
    M_ext[:, 1, 1] = 1.0

    Mc = M_ext.copy()  # running product

    for layer in wall.material_layers:
        M_layer = np.zeros((M_harmonics, 2, 2), dtype=complex)
        if layer.thickness == 0.0:
            M_layer[:, 0, 0] = 1.0
            M_layer[:, 0, 1] = layer.resistance
            M_layer[:, 1, 1] = 1.0
        else:
            alpha_p = np.sqrt(
                omegas * layer.specific_heat * layer.density
                * layer.conductivity / 2.0
            )
            beta = layer.thickness * alpha_p / layer.conductivity
            T_v = np.cosh(beta)
            u_v = np.sinh(beta)
            w_v = np.cos(beta)
            v_v = np.sin(beta)

            M_layer[:, 0, 0] = T_v * w_v + 1j * u_v * v_v
            M_layer[:, 0, 1] = ((u_v * w_v + T_v * v_v)
                                + 1j * (T_v * v_v - u_v * w_v)) / (2.0 * alpha_p)
            M_layer[:, 1, 0] = (alpha_p * (u_v * w_v - T_v * v_v)
                                + 1j * alpha_p * (T_v * v_v + u_v * w_v))
            M_layer[:, 1, 1] = M_layer[:, 0, 0]

        # Batched 2x2 complex matrix multiply: Mc = Mc @ M_layer
        Mc_new = np.zeros_like(Mc)
        Mc_new[:, 0, 0] = Mc[:, 0, 0] * M_layer[:, 0, 0] + Mc[:, 0, 1] * M_layer[:, 1, 0]
        Mc_new[:, 0, 1] = Mc[:, 0, 0] * M_layer[:, 0, 1] + Mc[:, 0, 1] * M_layer[:, 1, 1]
        Mc_new[:, 1, 0] = Mc[:, 1, 0] * M_layer[:, 0, 0] + Mc[:, 1, 1] * M_layer[:, 1, 0]
        Mc_new[:, 1, 1] = Mc[:, 1, 0] * M_layer[:, 0, 1] + Mc[:, 1, 1] * M_layer[:, 1, 1]
        Mc = Mc_new

    # Internal surface
    M_int = np.zeros((M_harmonics, 2, 2), dtype=complex)
    M_int[:, 0, 0] = 1.0
    M_int[:, 0, 1] = 1.0 / h_int
    M_int[:, 1, 1] = 1.0
    Mc_new = np.zeros_like(Mc)
    Mc_new[:, 0, 0] = Mc[:, 0, 0] * M_int[:, 0, 0] + Mc[:, 0, 1] * M_int[:, 1, 0]
    Mc_new[:, 0, 1] = Mc[:, 0, 0] * M_int[:, 0, 1] + Mc[:, 0, 1] * M_int[:, 1, 1]
    Mc_new[:, 1, 0] = Mc[:, 1, 0] * M_int[:, 0, 0] + Mc[:, 1, 1] * M_int[:, 1, 0]
    Mc_new[:, 1, 1] = Mc[:, 1, 0] * M_int[:, 0, 1] + Mc[:, 1, 1] * M_int[:, 1, 1]
    Mc = Mc_new

    # T_e / B_complex for each harmonic
    B_complex = Mc[:, 0, 1]  # (M,) complex
    T_complex = amplitude[1:] * np.exp(1j * phase[1:])
    Q_complex = T_complex / B_complex

    amp_qi = np.zeros(n_harmonics + 1)
    phase_qi = np.zeros(n_harmonics + 1)
    amp_qi[1:] = np.abs(Q_complex)
    phase_qi[1:] = np.angle(Q_complex)

    # --- Reconstruct output signal (vectorized) ---
    Qio = (A0 - T_int) * U

    istanti = np.arange(1, n_points + 1)
    omega_arr = (2.0 * math.pi / (n_points * 3600.0)) * m_arr  # (M,)
    # angles_out[m, istante] = omega_m * istante * 3600 + phase_qi_m
    angles_out = np.outer(omega_arr, istanti * 3600.0) + phase_qi[1:, np.newaxis]
    Fl = np.zeros(n_points + 1)
    Fl[1:] = (amp_qi[1:, np.newaxis] * np.sin(angles_out)).sum(axis=0) + Qio

    # Shift so Fl[0] = midnight
    Fl_shifted = np.zeros(n_points + 1)
    Fl_shifted[:n_points] = Fl[1:n_points + 1]
    Fl_shifted[n_points] = Fl_shifted[0]

    # Sample at the desired sampling time
    n_samples = int(24.0 / sampling_time) + 1
    flux_fourier = np.zeros(n_samples)
    points_per_hour = n_points / 24.0
    for h in range(n_samples):
        idx = int(h * sampling_time * points_per_hour)
        if idx < len(Fl_shifted):
            flux_fourier[h] = Fl_shifted[idx]

    return flux_fourier


def _interpolate_profile(profile_24h: np.ndarray, n_points: int) -> np.ndarray:
    """Interpolate a 24-hour profile to n_points using cubic splines.

    Uses the same cubic Hermite interpolation as the VB code.

    Args:
        profile_24h: Temperature values at hours 0..24 (25 values, with [0]=[24]).
        n_points: Number of interpolated points (typically 240).

    Returns:
        signal: Array indexed 0..n_points with interpolated values.
    """
    # Ensure we have 27 points for the interpolation (0..26)
    inp = np.zeros(27)
    inp[0:25] = profile_24h[0:25]
    inp[25] = profile_24h[1]
    inp[26] = profile_24h[2]

    # Compute cubic spline coefficients for each hour-long segment
    A_coeff = np.zeros(25)
    B_coeff = np.zeros(25)
    C_coeff = np.zeros(25)
    D_coeff = np.zeros(25)

    for n in range(1, 25):
        alpha1 = math.atan(inp[n] - inp[n - 1])
        alpha2 = math.atan(inp[n + 1] - inp[n])
        slope1 = math.tan((alpha1 + alpha2) / 2.0)

        alpha3 = math.atan(inp[n + 1] - inp[n])
        alpha4 = math.atan(inp[n + 2] - inp[n + 1])
        slope2 = math.tan((alpha3 + alpha4) / 2.0)

        # Solve 4x4 system for cubic ax³+bx²+cx+d passing through
        # (n, inp[n]) and (n+1, inp[n+1]) with slopes slope1 and slope2
        system = np.array([
            [n ** 3, n ** 2, n, 1],
            [(n + 1) ** 3, (n + 1) ** 2, n + 1, 1],
            [3 * n ** 2, 2 * n, 1, 0],
            [3 * (n + 1) ** 2, 2 * (n + 1), 1, 0]
        ], dtype=float)
        rhs = np.array([inp[n], inp[n + 1], slope1, slope2])

        try:
            sol = np.linalg.solve(system, rhs)
        except np.linalg.LinAlgError:
            sol = np.zeros(4)
            sol[3] = inp[n]

        A_coeff[n] = sol[0]
        B_coeff[n] = sol[1]
        C_coeff[n] = sol[2]
        D_coeff[n] = sol[3]

    # Evaluate the cubic splines at n_points equally spaced points
    step = 0.1  # 24h / 240 = 0.1h
    raw = np.zeros(int(25 / step) + 1)
    sample_idx = int(1.0 / step)  # = 10
    segment = 1
    count = 0

    t = 1.0
    while t <= 25.0 + step / 2:
        if count == int(1.0 / step):
            count = 0
            segment += 1
            if segment == 25:
                segment = 1
        raw[sample_idx] = (A_coeff[segment] * t ** 3 + B_coeff[segment] * t ** 2
                           + C_coeff[segment] * t + D_coeff[segment])
        count += 1
        sample_idx += 1
        t += step

    # Wrap around
    pts_per_hour = int(1.0 / step)
    raw[int(25 / step)] = raw[pts_per_hour]
    for i in range(pts_per_hour + 1):
        raw[i] = raw[int(24 / step) + i]

    # Build output signal (1-indexed to match the Fourier routine)
    signal = np.zeros(n_points + 1)
    for i in range(1, n_points + 1):
        signal[i] = raw[i - 1]

    return signal


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def compute_ctf(wall: Wall,
                n_roots: int = 50,
                n_poles: int | None = None,
                n_coefficients: int = 49,
                sampling_time: float = 1.0,
                n_periods: int = 10,
                n_harmonics: int = 120,
                T_int: float = 24.0,
                temperature_profile: np.ndarray | None = None,
                use_mitalas: bool = True,
                validate_fourier: bool = True) -> CTFResult:
    """Compute CTF coefficients for a multilayer wall.

    Args:
        wall: Wall definition with layers.
        n_roots: Number of roots of B(s) to find.
        n_poles: Number of poles to use for coefficients (default: n_roots).
        n_coefficients: Maximum number of z-domain coefficients.
        sampling_time: Sampling period in hours (must divide 24).
        n_periods: Number of 24h periods for transient simulation.
        n_harmonics: Number of harmonics for Fourier validation.
        T_int: Constant internal temperature [°C].
        temperature_profile: 24-hour sol-air temperature profile (25 values,
            index 0 = midnight, index 24 = midnight again). If None, Fourier
            validation is skipped.
        use_mitalas: Apply Mitalas instruction (recommended).
        validate_fourier: Compute Fourier validation error.

    Returns:
        CTFResult with all computed coefficients and validation metrics.
    """
    wall.validate()

    if n_poles is None:
        n_poles = n_roots

    # 1. Find roots of B(s) = 0
    poles, mat_at_roots, layer_mats_at_roots = _find_roots(wall, n_roots)

    # 2. Compute Heaviside expansion (with all roots)
    C0_1B, C0_AB, C1_1B, C1_AB, residues_1B, residues_AB = _compute_heaviside(
        wall, poles, mat_at_roots, layer_mats_at_roots, use_mitalas=False
    )

    # 3. Procedure I (Beccali et al. 2005): optimal selection of significant
    # poles and residues. Sort by |residue| descending and keep only those
    # with |residue| > sigma (absolute threshold).
    sigma = 1e-10
    res_magnitude = np.maximum(np.abs(residues_1B), np.abs(residues_AB))
    order = np.argsort(-res_magnitude)  # descending

    # Count significant residues
    n_significant = 0
    for i in range(n_roots):
        if res_magnitude[order[i]] > sigma:
            n_significant += 1
        else:
            break
    n_significant = max(n_significant, 1)
    n_significant = min(n_significant, n_poles)

    # Reorder poles and residues by significance
    sig_indices = order[:n_significant]
    poles_sig = poles[sig_indices]
    residues_1B_sig = residues_1B[sig_indices]
    residues_AB_sig = residues_AB[sig_indices]

    # Apply Mitalas instruction on the selected residues
    if use_mitalas:
        C1_1B = -np.sum(residues_1B_sig)
        C1_AB = -np.sum(residues_AB_sig)

    # 4. Compute denominator from the selected poles
    Tc_s = sampling_time * 3600.0  # hours -> seconds
    n_coefficients = min(n_coefficients, n_significant)
    denom = _compute_denominator(Tc_s, poles_sig, n_significant, n_coefficients)

    # Cap n_coefficients to the number of significant denominator terms
    n_eff = n_coefficients
    for k in range(n_coefficients, 0, -1):
        if abs(denom[k]) > 1e-10:
            n_eff = k
            break
    else:
        n_eff = 0
    n_coefficients = min(n_coefficients, max(n_eff, 1))

    # 5. Compute numerator coefficients using selected poles/residues
    b_coeffs, c_coeffs = _compute_numerators(
        Tc_s, poles_sig, n_significant, n_coefficients,
        C0_1B, C0_AB, C1_1B, C1_AB,
        residues_1B_sig, residues_AB_sig, denom
    )

    # 6. Build result
    result = CTFResult(
        b_coeffs=b_coeffs,
        c_coeffs=c_coeffs,
        d_coeffs=denom,
        n_poles=n_significant,
        n_coefficients=n_coefficients,
        sampling_time=sampling_time,
        thermal_transmittance=wall.thermal_transmittance,
        poles=poles_sig,
    )

    # 7. Fourier validation (if a temperature profile is provided)
    if validate_fourier and temperature_profile is not None:
        tau_max = int(24.0 / sampling_time) * n_periods
        n_per_day = int(24.0 / sampling_time)

        # Build sol-air temperature array
        T_ext = np.zeros(tau_max + 1)
        # profile_24h: index 0..24 hourly
        profile = temperature_profile
        T_ext[0] = profile[0]
        for period in range(n_periods):
            for i in range(1, n_per_day + 1):
                hour = int(i * sampling_time)
                idx = period * n_per_day + i
                if idx <= tau_max:
                    T_ext[idx] = profile[hour]

        # Compute CTF output
        q_ctf = _compute_signal_output(
            b_coeffs, c_coeffs, denom, T_ext, T_int, n_coefficients, tau_max,
            sampling_time
        )

        # Extract last period
        flux_ctf = np.zeros(n_per_day + 1)
        t_start = tau_max - n_per_day
        for i in range(n_per_day + 1):
            flux_ctf[i] = q_ctf[t_start + i]

        # Compute Fourier solution
        flux_fourier = _fourier_validate(
            wall, profile, T_int, n_harmonics, sampling_time
        )

        # Compute error
        n_compare = min(len(flux_ctf), len(flux_fourier))
        err = 0.0
        count = 0
        for k in range(1, n_compare):
            if abs(flux_fourier[k]) > 1e-10:
                err += abs((flux_fourier[k] - flux_ctf[k]) / flux_fourier[k])
                count += 1
        if count > 0:
            result.fourier_error = err / count * 100.0

    return result


def _compute_ctf_worker(args: tuple) -> CTFResult:
    """Worker for parallel batch computation (must be top-level for pickling)."""
    wall, kwargs = args
    return compute_ctf(wall, **kwargs)


def compute_ctf_batch(walls: list[Wall],
                      n_workers: int | None = None,
                      **kwargs) -> list[CTFResult]:
    """Compute CTF coefficients for multiple walls in parallel.

    Uses multiprocessing to distribute wall computations across CPU cores.

    Args:
        walls: List of Wall definitions.
        n_workers: Number of parallel workers (default: number of CPU cores).
        **kwargs: Additional arguments passed to compute_ctf for each wall.

    Returns:
        List of CTFResult, one per wall.
    """
    from concurrent.futures import ProcessPoolExecutor
    import os

    if n_workers is None:
        n_workers = os.cpu_count() or 1

    # For single wall or single worker, skip multiprocessing overhead
    if len(walls) <= 1 or n_workers <= 1:
        return [compute_ctf(wall, **kwargs) for wall in walls]

    work_items = [(wall, kwargs) for wall in walls]
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        results = list(executor.map(_compute_ctf_worker, work_items))

    return results
