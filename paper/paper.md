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
of Fundamentals for cooling and heating load calculations
[@Mitalas1967; @ASHRAE2021].

The thermal behaviour of a building wall is governed by the Fourier heat
equation, a partial differential equation whose direct solution in the time
domain is complex. The TFM addresses this by working in the Laplace domain,
where the heat conduction equations become algebraic. The analytical solution,
derived by Carslaw and Jaeger [@Carslaw1959], leads to a compact 2×2
transmission matrix for each wall layer. For a multilayer wall, the overall
transmission matrix is simply the ordered product of all layer matrices. The
inverse Laplace transform, performed via the Heaviside partial-fraction
expansion, converts the transfer function into the Z-domain through the
substitution $z = e^{s\Delta}$, yielding a ratio of polynomials in $z^{-1}$
whose coefficients encode the wall's complete thermal dynamic behaviour.

Given the thermophysical properties of a wall (layer thicknesses, densities,
specific heats, thermal conductivities), `wall-ctf` determines these
coefficients. Once computed, they are a **property of the wall itself** and can
be reused with any arbitrary input signal to compute the heat flux at the wall
surfaces through a simple recursive formula, making the method extremely
efficient for long simulations (e.g., a full year at hourly resolution).

# Statement of need

The Transfer Function Method was developed by Mitalas and Stephenson in the
early 1970s [@Mitalas1967] and has since become one of the most widely adopted
approaches for dynamic thermal simulation of buildings. However, the original
algorithm presents well-documented numerical limitations when applied to
massive building structures with high thermal inertia.

These limitations were first identified and systematically analysed in the
Ph.D. dissertation of Lo Brano at the Università degli Studi di Palermo, using
a case study of a historic building in Marsala (Sicily, Southern Italy). The
building, typical of the Mediterranean architectural heritage, featured
sandstone walls up to 900 mm thick with inside-outside plaster --- a
construction typology for which the ASHRAE standard coefficient tables, developed
primarily for North American lightweight construction, are simply not available.

Using control system analysis tools (step response, Bode diagrams, pole-zero
plots), the research demonstrated that:

1. For a sampling period $\Delta = 1$ h and 15 poles, the step response
   diverges to $10^8$ instead of converging to the correct steady-state value
   of 50°C. **Increasing the number of poles makes the results worse, not
   better.** This is the most counterintuitive finding.

2. The discrete-time transfer function exhibits **non-minimum phase zeros**
   (zeros outside the unit circle in the complex plane), which cause the
   initial transient response to move in the opposite direction of the
   excitation --- a physically impossible behaviour for a passive thermal
   system.

3. The Bode diagram magnitude at low frequencies exceeds 50 dB for
   $G_{15}(z)$ with $\Delta = 1$ h, implying a steady-state gain greater than
   300 --- whereas a thermal system can never amplify its input signal.

These findings, published in two peer-reviewed journal articles
[@Beccali2005zone; @Beccali2005reliable], led to three key guidelines:

- *The best model is not always the one with the largest number of poles.*
  A model with only 5 poles can represent the system as well as one with 15,
  without any numerical problem.
- *Increasing the sampling period is always a good way to eliminate numerical
  problems.* With $\Delta = 3$ h, all models behave correctly.
- *An automatic procedure for the optimal selection of significant poles and
  residues (Procedure I) gives a dramatic improvement* to the quality of
  simulations, achieving PME < 1% even for $\Delta = 1$ h with massive walls.

Despite the importance of these improvements, no open-source implementation
has been available to the scientific community until now. Existing tools that
implement the TFM are either **proprietary and closed-source** (TRNSYS, DOE-2),
**limited to pre-computed coefficient tables** (ASHRAE Handbook), or **not
available as standalone libraries**, being embedded in monolithic simulation
environments. `wall-ctf` fills this gap by providing a transparent, validated,
and extensible Python implementation.

# Method

The algorithm implemented in `wall-ctf` follows the mathematical framework
described in @Beccali2005zone and @Beccali2005reliable:

1. **Thermal transmission matrix**: For each homogeneous layer, the 2×2
   Laplace-domain transmission matrix is built from the Carslaw-Jaeger
   solution of the heat equation [@Carslaw1959]. The matrix elements involve
   hyperbolic functions of $L\sqrt{s/\alpha}$, where $L$ is the layer
   thickness and $\alpha$ is the thermal diffusivity. For a multilayer wall,
   the overall matrix is the ordered product of all layer matrices, including
   surface resistances.

2. **Root finding**: The poles of the transfer function (zeros of the $B(s)$
   element of the transmission matrix) are located on the negative real axis
   using an adaptive-step search algorithm with the substitution
   $\sqrt{s} = j\delta$, as described in the Ph.D. dissertation of Lo Brano.
   The step size is automatically adjusted based on the thermophysical
   properties of the layers to ensure no root is missed.

3. **Heaviside expansion**: The ramp response is decomposed into partial
   fractions, yielding the coefficients $C_0$, $C_1$, and the residues at
   each pole. The Mitalas instruction ($C_1 = -\sum \text{res}_n$) is applied
   to ensure numerical stability.

4. **Procedure I** (optimal pole/residue selection): Residues are sorted by
   magnitude and only those above a significance threshold
   ($\sigma = 10^{-10}$) are retained, following the approach described in
   @Beccali2005reliable. This prevents the numerical problems that arise from
   insignificant poles and dramatically improves the reliability of the
   simulation for all wall types.

5. **Z-domain coefficients**: The denominator polynomial is computed by
   iterative multiplication of $(1 - e^{s_n \Delta} z^{-1})$ factors, and
   the numerator is obtained from the ramp-response convolution. The number
   of effective coefficients is automatically capped based on the significance
   of the denominator terms.

6. **Fourier validation**: An independent harmonic solution based on the
   complex thermal quadrupole provides a reference for assessing the quality
   of the selected pole/residue set. This is not an alternative simulation
   method (it must be recomputed for every different input signal) but rather
   a quality check to verify that the truncated transfer function faithfully
   represents the wall's thermal dynamics.

The library supports parallel computation of multiple walls via Python's
`concurrent.futures` and provides a command-line interface for batch
processing with JSON input/output.

# Validation

The CTF output has been validated against the independent Fourier
(harmonic analysis) solution for various wall typologies. \autoref{fig:validation}
shows the comparison for a heavy concrete wall (plaster 2 cm + concrete 25 cm
+ plaster 2 cm), yielding a Percentage Mean Error (PME) of 0.15%.

![Comparison between the Z-transform (CTF) output and the Fourier reference
solution for a heavy concrete wall. The two curves are virtually
indistinguishable (PME = 0.15%), confirming that the selected set of poles
and residues accurately represents the wall's thermal dynamics. Once computed,
these coefficients can be reused with any arbitrary input
signal.\label{fig:validation}](../docs/ctf_vs_fourier.png)

Additional tests with lightweight walls (wood, 10 cm), two-layer walls
(brick + concrete), and various sampling periods confirm PME values
consistently below 5% for standard constructions.

# Acknowledgements

The algorithm was originally developed as part of the Ph.D. dissertation
of Valerio Lo Brano at the Università degli Studi di Palermo (Department of
Energy and Environmental Research - DREAM), under the supervision of
Prof. Giorgio Beccali and Prof. Aldo Orioli. The author acknowledges the
contribution of Prof. Maurizio Cellura to the journal publications that
formalised the method.

# References
