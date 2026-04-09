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

`wall-ctf` is a Python library for computing the Conduction Transfer Function
(CTF) coefficients of multilayer building walls. The CTF method (also called
Transfer Function Method, TFM) underpins programs like TRNSYS, DOE-2, BLAST
and TARP, and is the procedure recommended by the ASHRAE Handbook of
Fundamentals for cooling and heating load calculations [@Mitalas1967;
@ASHRAE2021].

Heat conduction through a wall is governed by the Fourier equation, a partial
differential equation that is hard to solve directly in the time domain. The
TFM sidesteps this difficulty by moving to the Laplace domain, where the
problem becomes algebraic. Following Carslaw and Jaeger [@Carslaw1959], each
homogeneous layer is described by a 2×2 transmission matrix; the overall wall
matrix is just the product of the individual ones. A Heaviside partial-fraction
expansion followed by the Z-transform turns the result into a ratio of
polynomials in $z^{-1}$. The coefficients of those polynomials are the CTF set.

The point is that these coefficients, once obtained, belong to the wall -- not
to the input signal. Run a year-long simulation at hourly resolution (8760
steps) and it takes a fraction of a second, because the expensive part was
already done when the coefficients were computed.

# Statement of need

Mitalas and Stephenson published the TFM in the early 1970s [@Mitalas1967].
More than fifty years later, it remains the backbone of the most widely
distributed building energy simulation tools. Yet, to the best of the author's
knowledge, there does not appear to be a standalone, open-source library that
lets a researcher compute CTF coefficients from arbitrary wall thermophysical
data and get the result as a reusable set of numbers.

Commercial programs (TRNSYS, DOE-2) embed the TFM in proprietary code that
cannot be inspected or called in isolation. EnergyPlus computes CTFs internally
with a state-space method, but the logic sits deep inside the simulation engine
and is not exposed as a separate module. The ASHRAE Handbook provides
pre-computed tables for a limited catalogue of North American wall typologies --
useful in the US, much less so for the sandstone-and-plaster walls of a
Sicilian palazzo or, more generally, for any wall that is not in the table.

There is also a numerical problem, and it is not obvious. When the original
algorithm is applied to massive walls -- say, 0.4 m of sandstone or more, the
kind you find all over Southern Europe -- adding more poles to the transfer
function makes the results worse, not better. The Percentage Mean Error can
jump from below 1% to over 1000%. The author ran into this during his Ph.D.
at the Università degli Studi di Palermo and spent considerable time tracking
down the cause: numerical ill-conditioning that plants non-minimum phase zeros
in the discrete transfer function [@Beccali2005zone; @Beccali2005reliable].
The research led to a few practical rules:

- The best model is not necessarily the one with the most poles. Five poles
  can be enough where fifteen fail.
- A larger sampling period (2--3 h instead of 1 h) removes most numerical
  difficulties, at the cost of time resolution.
- Sorting residues by magnitude and discarding the negligible ones (Procedure
  I) restores accuracy even at 1 h sampling, without manual tuning.

`wall-ctf` packages these insights into a tool aimed at building physicists,
energy engineers, and anyone who needs CTF coefficients without writing the
algorithm from scratch.

# State of the field

At the time of writing, there does not seem to be a pip-installable Python
package dedicated to CTF coefficient computation. A search of PyPI, GitHub,
CRAN and the MATLAB File Exchange turned up only one related project: FastCTF
[@Khalighi2021], a small C++ program that computes conduction transfer
functions via Padé approximants -- a different approach from the
Mitalas-Stephenson Z-transform method used here. FastCTF is not packaged as a
library, has no Python bindings, and has seen little activity since 2023.
General-purpose wrappers like `eppy` or `honeybee-energy` talk to EnergyPlus
but do not expose its internal CTF computation.

Compared with FastCTF and with what is buried inside commercial simulators,
`wall-ctf` offers a pip-installable library with no external dependencies
beyond NumPy, the full Mitalas-Stephenson algorithm (Heaviside expansion,
Mitalas instruction, ramp interpolation), Procedure I for automatic pole
selection, a built-in Fourier quality check, and parallel computation of
multiple walls.

# Software design

The code follows the mathematical derivation in @Beccali2005zone from start to
finish: build the transmission matrices, find the roots of $B(s)$, expand into
partial fractions, convert to the Z-domain, extract numerator and denominator
coefficients. Every function carries a reference to the relevant equation in
the source papers, so a reader can check the maths against the implementation
line by line.

Performance-sensitive loops are written with NumPy. The Fourier validation, for
instance, evaluates all harmonics at once through batched 2×2 complex matrix
products instead of looping in Python. The denominator polynomial is built by
iterative factor multiplication -- O($n_p^2$) -- rather than the
exponential-time recursion of the original VB.NET code. When several walls need
processing, `concurrent.futures` distributes the work across CPU cores.

From the user's side, a wall is a list of layers (Python objects or a JSON
file) and `compute_ctf()` returns the coefficients. Procedure I runs by
default, so the user does not need to guess how many poles are appropriate. A
command-line interface handles batch jobs.

One design choice deserves mention. The original algorithm was written in
Imperial units, where the root-finding step size of 0.01 worked fine. In SI
the roots are much closer together, and that same step misses every other one.
`wall-ctf` computes an adaptive step from the layer diffusivities and
thicknesses. This is not documented in the original publications because it
was not an issue in Imperial; it only surfaced during the Python rewrite.

# Research impact statement

The method behind `wall-ctf` has been in use since 2005. The two journal papers
[@Beccali2005zone; @Beccali2005reliable] have been cited in the building
simulation literature and contributed to the discussion on how the TFM performs
with massive European walls -- a question that matters for energy retrofitting
of historical buildings, a growing concern across Southern Europe.

The original software, THELDA (later CATI2005), was written in VB.NET and was
used to produce a database of Z-transform coefficients for Mediterranean
building typologies. That code was never released publicly. `wall-ctf` is a
complete rewrite in Python that makes the same algorithm available as an open
library, so that other groups can reproduce and build on the results.

# AI usage disclosure

The mathematical algorithm was developed and validated entirely in the author's
Ph.D. dissertation and in the referenced journal publications, without
assistance from AI tools.

# References
