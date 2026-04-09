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
    corresponding: true
    affiliation: 1
affiliations:
  - name: Department of Engineering, Università degli Studi di Palermo, Viale delle Scienze, 90128 Palermo, Italy
    index: 1
date: 9 April 2026
bibliography: paper.bib
---

# Summary

`wall-ctf` is a Python library that computes the Conduction Transfer Function
(CTF) coefficients for multilayer building wall assemblies using the
Z-transform method. The CTF method, also known as the Transfer Function Method
(TFM), is the mathematical backbone of widely used building energy simulation
programs and is recommended by the ASHRAE Handbook of Fundamentals for cooling
and heating load calculations [@Mitalas1967; @ASHRAE2021].

The thermal behaviour of a building wall is governed by the Fourier heat
equation. The TFM addresses this by working in the Laplace domain, where the
heat conduction equations become algebraic. The analytical solution, derived by
Carslaw and Jaeger [@Carslaw1959], leads to a compact 2×2 transmission matrix
for each wall layer. For a multilayer wall, the overall transmission matrix is
the ordered product of all layer matrices. Through the Heaviside
partial-fraction expansion and the Z-transform, the continuous-time transfer
function is converted into a ratio of polynomials in $z^{-1}$ whose
coefficients constitute the CTF set.

Given the thermophysical properties of a wall (layer thicknesses, densities,
specific heats, thermal conductivities), `wall-ctf` determines these
coefficients. Once computed, they are a property of the wall itself and can be
reused with any arbitrary input signal to compute the heat flux at the wall
surfaces through a simple recursive formula, making the method extremely
efficient for long simulations (e.g., a full year at hourly resolution with
8760 time steps).

# Statement of need

Despite being developed in the early 1970s [@Mitalas1967] and forming the
basis of the most widely used building energy simulation tools, the Transfer
Function Method has remained inaccessible as a standalone computational tool.
Commercial software such as TRNSYS and DOE-2 implement the TFM internally, but
their CTF computation modules are proprietary and cannot be used independently.
EnergyPlus computes CTFs using the state-space method, but the logic is deeply
embedded in the simulation engine. The ASHRAE Handbook provides pre-computed
coefficient tables only for a limited set of North American wall typologies,
inadequate for custom wall definitions or for the diverse construction
typologies found in European and Mediterranean buildings.

Beyond the lack of open-source tools, the original Mitalas-Stephenson
algorithm presents well-documented numerical limitations when applied to
massive building walls with high thermal inertia. As demonstrated in the Ph.D.
dissertation of Lo Brano and subsequently published in two journal articles
[@Beccali2005zone; @Beccali2005reliable], a naive application of the TFM to
walls with thickness above 0.4 m can produce counterintuitive results:
increasing the number of poles in the transfer function worsens the accuracy
instead of improving it, with Percentage Mean Errors (PME) escalating from
less than 1% to over 1000%. This behaviour is caused by numerical
ill-conditioning that introduces non-minimum phase zeros and corrupts the
low-frequency response of the discrete transfer function. The research
identified three fundamental guidelines: (1) the best model is not always the
largest in terms of number of poles; (2) increasing the sampling period
eliminates most numerical problems; (3) an automatic procedure for the optimal
selection of significant poles and residues (Procedure I) dramatically improves
simulation reliability.

`wall-ctf` is aimed at building physicists, energy engineers, and researchers
in building simulation who need a transparent, validated, and extensible tool
for computing CTF coefficients for arbitrary wall compositions.

# State of the field

A survey of existing open-source software reveals that no pip-installable
Python library exists for computing CTF coefficients from wall thermophysical
properties. The only related open-source project is FastCTF [@Khalighi2021], a
minimal C++ program that uses Padé approximants (a different mathematical
approach from the Mitalas-Stephenson Z-transform method) and is not available
as a library or package. General-purpose building simulation wrappers such as
`eppy` or `honeybee-energy` interface with EnergyPlus but do not expose CTF
computation as a standalone function. No equivalent tool exists in R, MATLAB,
or Fortran.

`wall-ctf` differs from these existing approaches by providing: (a) a
standalone, pip-installable library with no dependency on external simulation
engines; (b) the complete Mitalas-Stephenson algorithm including the Heaviside
expansion, the Mitalas instruction, and the ramp-interpolation method for
numerator computation; (c) the Procedure I improvement for optimal pole/residue
selection [@Beccali2005reliable]; (d) built-in Fourier validation for quality
assessment of the computed coefficients; (e) parallel batch computation for
multiple walls.

# Software design

`wall-ctf` is designed around three principles: transparency, efficiency, and
ease of use.

**Transparency.** The implementation follows the mathematical derivation
presented in @Beccali2005zone step by step: transmission matrix construction,
root finding, Heaviside expansion, Z-domain coefficient computation. Each
function is documented with references to the specific equations in the source
publications. This allows researchers to inspect, verify, and extend every
stage of the algorithm.

**Efficiency.** The computationally intensive parts are vectorised with NumPy.
The Fourier validation uses batched complex matrix operations across all
harmonics simultaneously. The denominator polynomial is computed by iterative
multiplication (O($n_p^2$)) rather than the exponential-time recursive
approach used in the original implementation. Multiple walls can be computed in
parallel using Python's `concurrent.futures`.

**Ease of use.** Walls are defined as simple Python objects or JSON files. A
single function call (`compute_ctf`) returns all coefficients. The library
automatically selects significant poles (Procedure I) and caps the number of
effective coefficients based on the denominator's numerical significance. A
command-line interface is provided for batch processing.

A key design decision was the adoption of SI units throughout, which required
solving a subtle numerical issue: the original algorithm (developed in
Imperial units) used a fixed root-finding step size that, when applied to SI
units, was too coarse to resolve closely-spaced roots for multilayer walls.
`wall-ctf` addresses this with an adaptive step size computed from the
thermophysical properties of the individual layers.

# Research impact statement

The mathematical method implemented in `wall-ctf` has been used in published
research since 2005. The two foundational papers [@Beccali2005zone;
@Beccali2005reliable] have been cited in the building simulation literature and
contributed to the understanding of the numerical limitations of the TFM for
massive European building walls. The original software (THELDA/CATI2005,
written in VB.NET) was used to simulate the thermal behaviour of historical
buildings in the Mediterranean area and to build a database of Z-transform
coefficients for Mediterranean building typologies. The release of `wall-ctf`
as an open-source Python library makes this method accessible to the broader
research community for the first time, enabling reproducible research and
integration with modern Python-based building simulation workflows.

# AI usage disclosure

No artificial intelligence tools were used in the development of the
mathematical algorithm, which was entirely derived and validated in the Ph.D.
dissertation of the author and in the referenced publications.

# References
