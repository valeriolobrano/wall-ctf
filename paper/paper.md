---
title: 'wall-ctf: A Python library for computing Conduction Transfer Function coefficients of multilayer building walls'
tags:
  - Python
  - building simulation
  - heat transfer
  - transfer function method
  - Z-transform
  - thermal inertia
authors:
  - name: Valerio Lo Brano
    orcid: 0000-0002-5728-0379
    affiliation: 1
affiliations:
  - name: Department of Engineering, Università degli Studi di Palermo, Italy
    index: 1
date: 9 April 2026
bibliography: paper.bib
---

# Summary

`wall-ctf` is a Python library that computes the Conduction Transfer Function
(CTF) coefficients for multilayer building wall assemblies using the
Z-transform method. The CTF method, also known as the Transfer Function Method
(TFM), is the mathematical backbone of widely used building energy simulation
programs such as TRNSYS, DOE-2, BLAST, and TARP, and is recommended by the
ASHRAE Handbook of Fundamentals for cooling and heating load calculations
[@Mitalas1967; @ASHRAE2021].

The thermal behaviour of a building wall is governed by the Fourier heat
equation, whose solution in the Laplace domain --- following the approach of
Carslaw and Jaeger [@Carslaw1959] --- leads to a compact 2×2 transmission
matrix for each wall layer. For a multilayer wall, the overall matrix is the
ordered product of all layer matrices. Through the Heaviside partial-fraction
expansion and the Z-transform, the continuous-time transfer function is
converted into a ratio of polynomials in $z^{-1}$, whose coefficients
constitute the CTF set.

Given the thermophysical properties of a wall (layer thicknesses, densities,
specific heats, thermal conductivities), `wall-ctf` determines these
coefficients. Once computed, they are a **property of the wall itself** and can
be reused with any arbitrary input signal to compute the heat flux at the wall
surfaces through a simple recursive formula, making the method extremely
efficient for long simulations.

# Statement of need

Despite being developed in the early 1970s [@Mitalas1967] and forming the
basis of the most widely used building energy simulation tools, the Transfer
Function Method has remained inaccessible as a standalone computational tool.
A survey of existing software reveals that:

- **No pip-installable Python library** exists for computing CTF coefficients
  from wall thermophysical properties. The only related open-source project is
  FastCTF [@Khalighi2021], a minimal C++ program (4 GitHub stars) that uses a
  different mathematical approach (Padé approximants) and is not available as
  a library or package.
- **Commercial tools** such as TRNSYS and DOE-2 implement the TFM internally,
  but their CTF computation modules are proprietary, closed-source, and cannot
  be inspected, modified, or used independently of the full simulation
  environment.
- **EnergyPlus** computes CTFs internally using the state-space method, but the
  logic is deeply embedded in the simulation engine and cannot be called as a
  standalone function.
- **ASHRAE Handbook tables** provide pre-computed coefficients only for a
  limited set of North American wall typologies, making them inadequate for
  custom wall definitions or for the diverse construction typologies found in
  European, Mediterranean, and historical buildings.

Beyond the lack of open-source tools, the original Mitalas-Stephenson
algorithm presents well-documented numerical limitations when applied to
massive building walls. As demonstrated in the Ph.D. dissertation of Lo Brano
at the Università degli Studi di Palermo and subsequently published in two
journal articles [@Beccali2005zone; @Beccali2005reliable], a naive application
of the TFM to walls with high thermal inertia (thickness > 0.4 m, typical of
Southern European and Mediterranean historical buildings) can produce
**counterintuitive results**: increasing the number of poles in the transfer
function worsens the accuracy instead of improving it, with Percentage Mean
Errors (PME) escalating from less than 1% to over 1000%.

This behaviour is caused by numerical ill-conditioning that introduces
non-minimum phase zeros (zeros outside the unit circle in the z-plane) and
corrupts the low-frequency response of the transfer function. Using control
system analysis tools --- step response, Bode diagrams, and pole-zero
location plots --- the research identified three fundamental rules
[@Beccali2005reliable]:

1. **The best model is not always the largest.** A model with only 5 poles can
   represent the system as well as one with 15 poles, without any numerical
   problems. Only at high frequencies is the behaviour more approximated.
2. **Increasing the sampling period is a reliable way to eliminate numerical
   problems.** This is however limited by the time resolution of available
   climatic data and the desired output resolution.
3. **Procedure I** --- an automatic algorithm that sorts residues by
   magnitude and retains only the significant ones ($|\hat{d}_n| > 10^{-10}$)
   --- gives a dramatic improvement to the quality of simulations, achieving
   PME < 1% even for massive walls with $\Delta = 1$ h sampling period.

`wall-ctf` fills the gap by providing the first open, transparent, validated,
and pip-installable Python implementation of the complete TFM algorithm,
including the Procedure I improvement.

# Implementation

The library implements the following algorithmic pipeline:

1. **Thermal transmission matrix**: For each homogeneous layer, the 2×2
   Laplace-domain transmission matrix is built from the Carslaw-Jaeger
   analytical solution. The matrix elements involve hyperbolic functions of
   $L\sqrt{s/\alpha}$, where $L$ is the layer thickness and $\alpha$ is the
   thermal diffusivity. For air gaps and surface resistance layers, the matrix
   reduces to $[[1, R], [0, 1]]$.

2. **Adaptive root finding**: The poles (zeros of $B(s)$ on the negative real
   axis) are located using the substitution $\sqrt{s} = j\delta$ with an
   adaptive step size that is automatically computed from the thermophysical
   properties of the layers. This ensures no root is missed even for
   multilayer walls with very different material properties.

3. **Heaviside expansion**: The ramp response is decomposed into partial
   fractions, yielding $C_0$, $C_1$, and the residues at each pole. The
   Mitalas instruction ($C_1 = -\sum \text{res}_n$) ensures the response
   starts at zero for $t = 0$.

4. **Procedure I**: Residues are sorted by absolute magnitude and only those
   above a significance threshold are retained, following
   @Beccali2005reliable.

5. **Z-domain coefficients**: The denominator is computed by iterative
   polynomial multiplication and the numerator from the ramp-response
   convolution. The number of effective coefficients is automatically capped
   based on the significance of the denominator terms.

6. **Fourier validation**: An independent harmonic solution (complex thermal
   quadrupole) provides a quality check --- not an alternative simulation
   method, since it must be recomputed for every input signal.

The library supports parallel computation of multiple walls via Python's
`concurrent.futures`, provides a command-line interface with JSON input/output,
and uses NumPy vectorisation throughout for performance.

# Validation

The CTF output has been validated against the independent Fourier
(harmonic analysis) solution for various wall typologies.
\autoref{fig:validation} shows the comparison for a heavy concrete wall
(plaster 2 cm + concrete 25 cm + plaster 2 cm), yielding PME = 0.15%.

![Comparison between the Z-transform (CTF) output and the Fourier reference
solution for a heavy concrete wall. The two curves are virtually
indistinguishable (PME = 0.15%), confirming that the selected set of poles
and residues accurately represents the wall's thermal dynamics.
\label{fig:validation}](../docs/ctf_vs_fourier.png)

Additional tests with lightweight walls (wood, 10 cm), two-layer walls
(brick + concrete), and various sampling periods confirm PME values
consistently below 5% for standard constructions.

# References
