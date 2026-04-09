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
(CTF) coefficients for multilayer building wall assemblies. The CTF method,
also known as the Z-transform method or Transfer Function Method (TFM), is the
mathematical backbone of widely used building energy simulation programs such
as TRNSYS, DOE-2, BLAST, and TARP, and is recommended by the ASHRAE Handbook
of Fundamentals for cooling and heating load calculations.

Given the thermophysical properties of a wall (layer thicknesses, densities,
specific heats, thermal conductivities), the library determines a compact set
of numerical coefficients that encode the wall's complete thermal dynamic
behaviour. Once computed, these coefficients can be reused with any arbitrary
input signal to compute the heat flux at the wall surfaces through a simple
recursive formula, making the method extremely efficient for long simulations
(e.g., a full year at hourly resolution).

# Statement of need

The Transfer Function Method was developed by Mitalas and Stephenson in the
early 1970s [@Mitalas1967] and has since become one of the most widely adopted
approaches for dynamic thermal simulation of buildings. However, the original
algorithm presents well-documented numerical limitations when applied to
massive building structures with high thermal inertia, typical of the
Mediterranean and Southern European architectural heritage
[@Beccali2005reliable].

These limitations manifest as a counterintuitive degradation of accuracy when
increasing the number of poles in the transfer function: the Percentage Mean
Error (PME) can increase from less than 1% to over 1000% when using 15 poles
instead of 5 for walls with thickness above 0.6 m. This behaviour, caused by
non-minimum phase zeros and numerical ill-conditioning, was first identified
and analysed in the Ph.D. dissertation of Lo Brano at the Università degli
Studi di Palermo, and later published in two peer-reviewed journal articles
[@Beccali2005zone; @Beccali2005reliable].

The research led to the development of **Procedure I**, an automatic algorithm
for the optimal selection of significant poles and residues that dramatically
improves the reliability of the simulation for all wall types. Despite the
importance of this improvement, no open-source implementation has been
available to the scientific community until now.

Existing tools that implement the TFM are either:

- **Proprietary and closed-source** (TRNSYS, DOE-2), preventing researchers
  from inspecting, modifying, or extending the algorithm;
- **Limited to pre-computed coefficient tables** (ASHRAE Handbook), which cover
  only a restricted set of North American wall typologies and do not allow
  custom wall definitions;
- **Not available as standalone libraries**, being embedded in monolithic
  simulation environments.

`wall-ctf` fills this gap by providing an open, transparent, and validated
Python implementation of the complete TFM algorithm, including the Procedure I
improvement, usable both as a library and as a command-line tool.

# Method

The algorithm implemented in `wall-ctf` follows the mathematical framework
described in @Beccali2005zone and @Beccali2005reliable:

1. **Thermal transmission matrix**: For each homogeneous layer, the 2x2
   Laplace-domain transmission matrix is built from the Carslaw-Jaeger
   solution of the heat equation [@Carslaw1959]. The overall wall matrix is
   the ordered product of all layer matrices.

2. **Root finding**: The poles of the transfer function (zeros of the B(s)
   element of the transmission matrix) are located on the negative real axis
   using an adaptive-step search algorithm with the substitution
   $\sqrt{s} = j\delta$.

3. **Heaviside expansion**: The ramp response is decomposed into partial
   fractions, yielding the coefficients $C_0$, $C_1$, and the residues
   $\text{res}_n$ at each pole. The Mitalas instruction is applied to ensure
   numerical stability.

4. **Procedure I**: Residues are sorted by magnitude and only those above a
   significance threshold ($\sigma = 10^{-10}$) are retained, following the
   approach described in @Beccali2005reliable.

5. **Z-domain coefficients**: The denominator polynomial is computed by
   iterative multiplication of $(1 - e^{s_n \Delta} z^{-1})$ factors, and
   the numerator is obtained from the ramp-response convolution.

6. **Fourier validation**: An independent harmonic solution based on the
   complex thermal quadrupole provides a reference for assessing the quality
   of the selected pole/residue set.

The library supports parallel computation of multiple walls via Python's
`concurrent.futures` and provides a command-line interface for batch
processing with JSON input/output.

# Validation

The CTF output has been validated against the independent Fourier
(harmonic analysis) solution for various wall typologies. The Percentage Mean
Error (PME) is consistently below 1% for standard construction walls with
hourly sampling period, confirming that the truncated set of poles and residues
faithfully represents the wall's thermal dynamics (see \autoref{fig:validation}).

![Comparison between the Z-transform (CTF) output and the Fourier reference
solution for a heavy concrete wall (plaster 2 cm + concrete 25 cm + plaster
2 cm). PME = 0.15%.\label{fig:validation}](../docs/ctf_vs_fourier.png)

# Acknowledgements

The algorithm was originally developed as part of the Ph.D. dissertation
of Valerio Lo Brano at the Università degli Studi di Palermo (Department of
Energy and Environmental Research - DREAM), under the supervision of
Prof. Giorgio Beccali and Prof. Aldo Orioli.

# References
