"""Wall and layer data model for multilayer wall assemblies.

Copyright (c) 2005-2026 Valerio Lo Brano, Universita degli Studi di Palermo.
License: CC BY-NC 4.0 (non-commercial use).
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class Layer:
    """A single layer in a wall assembly.

    For material layers: specify thickness, density, specific_heat, conductivity.
    For air gaps: set thickness=0 and specify resistance.
    For surface resistances (first/last layer): set thickness=0 and specify resistance.

    All values in SI units:
        thickness: m
        density: kg/m³
        specific_heat: J/(kg·K)
        conductivity: W/(m·K)
        resistance: m²·K/W (used for air gaps and surface resistances)
    """

    name: str = ""
    thickness: float = 0.0
    density: float = 0.0
    specific_heat: float = 0.0
    conductivity: float = 0.0
    resistance: float = 0.0

    @property
    def is_surface_resistance(self) -> bool:
        return self.thickness == 0.0 and self.conductivity == 0.0

    @property
    def thermal_resistance(self) -> float:
        if self.thickness == 0.0:
            return self.resistance
        return self.thickness / self.conductivity

    @property
    def thermal_diffusivity(self) -> float:
        if self.thickness == 0.0 or self.density == 0.0:
            return 0.0
        return self.conductivity / (self.density * self.specific_heat)


@dataclass
class Wall:
    """A multilayer wall assembly.

    The first layer must be the external surface resistance,
    the last layer must be the internal surface resistance,
    and intermediate layers are material layers or air gaps
    (ordered from external to internal).
    """

    layers: list[Layer] = field(default_factory=list)
    name: str = ""

    @property
    def n_layers(self) -> int:
        return len(self.layers)

    @property
    def material_layers(self) -> list[Layer]:
        """Return only the material/air-gap layers (excluding surface resistances)."""
        return self.layers[1:-1]

    @property
    def external_surface_resistance(self) -> float:
        return self.layers[0].resistance

    @property
    def internal_surface_resistance(self) -> float:
        return self.layers[-1].resistance

    @property
    def total_thermal_resistance(self) -> float:
        return sum(layer.thermal_resistance for layer in self.layers)

    @property
    def thermal_transmittance(self) -> float:
        return 1.0 / self.total_thermal_resistance

    def validate(self) -> None:
        if len(self.layers) < 3:
            raise ValueError(
                "A wall must have at least 3 layers: "
                "external surface resistance, one material layer, "
                "internal surface resistance."
            )
        if not self.layers[0].is_surface_resistance:
            raise ValueError("First layer must be external surface resistance.")
        if not self.layers[-1].is_surface_resistance:
            raise ValueError("Last layer must be internal surface resistance.")

    @classmethod
    def from_dict(cls, data: dict) -> Wall:
        layers = [Layer(**layer_data) for layer_data in data["layers"]]
        return cls(layers=layers, name=data.get("name", ""))

    @classmethod
    def from_json(cls, path: str | Path) -> Wall:
        with open(path) as f:
            data = json.load(f)
        return cls.from_dict(data)

    def to_dict(self) -> dict:
        from dataclasses import asdict

        return {"name": self.name, "layers": [asdict(layer) for layer in self.layers]}
