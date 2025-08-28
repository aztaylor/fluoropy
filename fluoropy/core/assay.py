"""
Base assay classes for fluorescence-based experiments.
"""

from typing import Dict, List, Optional, Union
from abc import ABC, abstractmethod
import numpy as np

from .plate import Plate, Well
from .fluorophore import Fluorophore, fluorophore_db


class Assay(ABC):
    """Base class for fluorescence assays."""

    def __init__(self,
                 name: str,
                 fluorophore: Union[str, Fluorophore],
                 plate: Optional[Plate] = None):
        """
        Initialize an assay.

        Parameters
        ----------
        name : str
            Name of the assay
        fluorophore : str or Fluorophore
            Fluorophore used in the assay
        plate : Plate, optional
            Plate containing the assay data
        """
        self.name = name
        self.plate = plate

        # Handle fluorophore input
        if isinstance(fluorophore, str):
            self.fluorophore = fluorophore_db.get_fluorophore(fluorophore)
            if self.fluorophore is None:
                raise ValueError(f"Unknown fluorophore: {fluorophore}")
        else:
            self.fluorophore = fluorophore

        self.metadata = {}

    @abstractmethod
    def analyze(self) -> Dict:
        """Analyze the assay data. Must be implemented by subclasses."""
        pass

    def add_plate(self, plate: Plate) -> None:
        """Add a plate to the assay."""
        self.plate = plate

    def get_control_wells(self, control_type: str = "control") -> List[Well]:
        """Get control wells from the plate."""
        if self.plate is None:
            return []
        return self.plate.get_wells_by_type(control_type)

    def calculate_background(self,
                           blank_wells: Optional[List[str]] = None) -> float:
        """
        Calculate background fluorescence from blank wells.

        Parameters
        ----------
        blank_wells : list of str, optional
            List of blank well positions. If None, uses wells marked as 'blank'

        Returns
        -------
        float
            Mean background fluorescence
        """
        if self.plate is None:
            return 0.0

        if blank_wells is None:
            blanks = self.plate.get_wells_by_type("blank")
        else:
            blanks = [self.plate.get_well(pos) for pos in blank_wells
                     if self.plate.get_well(pos) is not None]

        if not blanks:
            return 0.0

        blank_values = []
        for well in blanks:
            if well.fluorescence is not None:
                if isinstance(well.fluorescence, list):
                    blank_values.append(well.fluorescence[-1])  # Latest timepoint
                else:
                    blank_values.append(well.fluorescence)

        return np.mean(blank_values) if blank_values else 0.0

    def subtract_background(self, background: Optional[float] = None) -> None:
        """
        Subtract background from all wells.

        Parameters
        ----------
        background : float, optional
            Background value to subtract. If None, calculates from blank wells
        """
        if self.plate is None:
            return

        if background is None:
            background = self.calculate_background()

        for well in self.plate.wells.values():
            if well.fluorescence is not None:
                if isinstance(well.fluorescence, list):
                    well.fluorescence = [f - background for f in well.fluorescence]
                else:
                    well.fluorescence = well.fluorescence - background

    def __repr__(self) -> str:
        return f"Assay({self.name}, {self.fluorophore.name if self.fluorophore else 'No fluorophore'})"


class EndpointAssay(Assay):
    """Simple endpoint fluorescence assay."""

    def analyze(self) -> Dict:
        """
        Analyze endpoint assay data.

        Returns
        -------
        Dict
            Analysis results including statistics and background info
        """
        if self.plate is None:
            return {"error": "No plate data available"}

        # Get sample wells
        sample_wells = self.plate.get_wells_by_type("sample")
        controls = self.plate.get_wells_by_type("control")

        # Calculate statistics
        sample_values = [w.fluorescence for w in sample_wells
                        if w.fluorescence is not None]
        control_values = [w.fluorescence for w in controls
                         if w.fluorescence is not None]

        results = {
            "assay_type": "endpoint",
            "fluorophore": self.fluorophore.name if self.fluorophore else None,
            "n_samples": len(sample_values),
            "n_controls": len(control_values),
            "background": self.calculate_background(),
        }

        if sample_values:
            results.update({
                "sample_mean": np.mean(sample_values),
                "sample_std": np.std(sample_values),
                "sample_cv": np.std(sample_values) / np.mean(sample_values) * 100
            })

        if control_values:
            results.update({
                "control_mean": np.mean(control_values),
                "control_std": np.std(control_values),
            })

        return results
