"""
Fluorophore properties and database for fluorescence assays.
"""

from typing import Dict, Optional
from dataclasses import dataclass


@dataclass
class Fluorophore:
    """Represents a fluorophore with its optical properties."""

    name: str
    excitation_max: float  # nm
    emission_max: float    # nm
    extinction_coeff: Optional[float] = None  # M^-1 cm^-1
    quantum_yield: Optional[float] = None
    brightness: Optional[float] = None  # extinction_coeff * quantum_yield

    def __post_init__(self):
        """Calculate brightness if not provided."""
        if (self.brightness is None and
            self.extinction_coeff is not None and
            self.quantum_yield is not None):
            self.brightness = self.extinction_coeff * self.quantum_yield


class FluorophoreDB:
    """Database of common fluorophores used in biological assays."""

    def __init__(self):
        self._fluorophores = self._load_default_fluorophores()

    def _load_default_fluorophores(self) -> Dict[str, Fluorophore]:
        """Load default fluorophore database."""
        fluorophores = {
            # Green fluorophores
            "FITC": Fluorophore(
                name="FITC",
                excitation_max=495,
                emission_max=519,
                extinction_coeff=68000,
                quantum_yield=0.85
            ),
            "GFP": Fluorophore(
                name="GFP",
                excitation_max=488,
                emission_max=507,
                extinction_coeff=55000,
                quantum_yield=0.60
            ),
            "Alexa488": Fluorophore(
                name="Alexa Fluor 488",
                excitation_max=495,
                emission_max=519,
                extinction_coeff=71000,
                quantum_yield=0.92
            ),

            # Red fluorophores
            "RFP": Fluorophore(
                name="RFP",
                excitation_max=558,
                emission_max=583,
                extinction_coeff=100000,
                quantum_yield=0.25
            ),
            "Alexa594": Fluorophore(
                name="Alexa Fluor 594",
                excitation_max=590,
                emission_max=617,
                extinction_coeff=92000,
                quantum_yield=0.66
            ),
            "Cy5": Fluorophore(
                name="Cy5",
                excitation_max=649,
                emission_max=670,
                extinction_coeff=250000,
                quantum_yield=0.27
            ),

            # Blue fluorophores
            "DAPI": Fluorophore(
                name="DAPI",
                excitation_max=358,
                emission_max=461,
                extinction_coeff=27000,
                quantum_yield=0.28
            ),
            "Hoechst": Fluorophore(
                name="Hoechst 33342",
                excitation_max=350,
                emission_max=461,
                extinction_coeff=46000,
                quantum_yield=0.38
            ),

            # Common assay dyes
            "Calcein": Fluorophore(
                name="Calcein",
                excitation_max=495,
                emission_max=515,
                extinction_coeff=75000,
                quantum_yield=0.77
            ),
            "Resazurin": Fluorophore(
                name="Resazurin",
                excitation_max=563,
                emission_max=587,
                extinction_coeff=54000,
                quantum_yield=0.074
            ),
        }

        return fluorophores

    def get_fluorophore(self, name: str) -> Optional[Fluorophore]:
        """Get a fluorophore by name."""
        return self._fluorophores.get(name.upper())

    def add_fluorophore(self, fluorophore: Fluorophore) -> None:
        """Add a custom fluorophore to the database."""
        self._fluorophores[fluorophore.name.upper()] = fluorophore

    def list_fluorophores(self) -> Dict[str, Fluorophore]:
        """Get all available fluorophores."""
        return self._fluorophores.copy()

    def find_by_wavelength(self,
                          excitation: Optional[float] = None,
                          emission: Optional[float] = None,
                          tolerance: float = 20) -> Dict[str, Fluorophore]:
        """
        Find fluorophores by excitation/emission wavelengths.

        Parameters
        ----------
        excitation : float, optional
            Target excitation wavelength (nm)
        emission : float, optional
            Target emission wavelength (nm)
        tolerance : float, default 20
            Wavelength tolerance (nm)

        Returns
        -------
        Dict[str, Fluorophore]
            Matching fluorophores
        """
        matches = {}

        for name, fluorophore in self._fluorophores.items():
            match = True

            if excitation is not None:
                if abs(fluorophore.excitation_max - excitation) > tolerance:
                    match = False

            if emission is not None:
                if abs(fluorophore.emission_max - emission) > tolerance:
                    match = False

            if match:
                matches[name] = fluorophore

        return matches


# Global fluorophore database instance
fluorophore_db = FluorophoreDB()
