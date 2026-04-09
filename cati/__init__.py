"""CATI - Conduction Transfer Function coefficients for multilayer walls.

Computes the Z-transform transfer function coefficients for heat conduction
through multilayer wall assemblies, using the method described in:

    [1] G. Beccali, M. Cellura, V. Lo Brano, A. Orioli,
        "Single thermal zone balance solved by Transfer Function Method",
        Energy and Buildings, 37 (2005) 1268-1277.
        DOI: 10.1016/j.enbuild.2005.02.010

    [2] G. Beccali, M. Cellura, V. Lo Brano, A. Orioli,
        "Is the transfer function method reliable in a European building
        context? A theoretical analysis and a case study in the south
        of Italy", Applied Thermal Engineering, 25 (2005) 341-357.
        DOI: 10.1016/j.applthermaleng.2004.06.010

Originally developed as part of the Ph.D. dissertation of Valerio Lo Brano
at the Università degli Studi di Palermo (DREAM).

License: CC BY-NC 4.0 (non-commercial use).
For commercial licensing contact: valerio.lobrano@unipa.it

Copyright (c) 2005-2026 Valerio Lo Brano, Università degli Studi di Palermo.
"""

from cati.wall import Wall, Layer
from cati.ctf import compute_ctf, compute_ctf_batch, CTFResult

__all__ = ["Wall", "Layer", "compute_ctf", "compute_ctf_batch", "CTFResult"]
__version__ = "1.0.8"
__author__ = "Valerio Lo Brano"
__email__ = "valerio.lobrano@unipa.it"
