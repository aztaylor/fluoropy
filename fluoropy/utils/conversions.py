"""
Unit conversion and data transformation utilities.
"""

from typing import Optional, Union

import numpy as np

# Each unit as a multiple of its family's base unit.
_MOLAR = {"M": 1.0, "mM": 1e-3, "µM": 1e-6, "nM": 1e-9, "pM": 1e-12}  # mol/L
_MASS = {  # g/L
    "g/L": 1.0,
    "mg/mL": 1.0,
    "mg/L": 1e-3,
    "µg/mL": 1e-3,
    "µg/L": 1e-6,
    "ng/mL": 1e-6,
    "ng/L": 1e-9,
    "pg/mL": 1e-9,
}


def _canonical_unit(unit: str) -> str:
    """Spell micro one way: the micro sign, Greek mu and plain 'u' all match."""
    unit = unit.strip().replace("\u03bc", "µ")  # Greek mu -> micro sign
    if unit.startswith("u") and unit[1:2] in ("M", "g"):
        unit = "µ" + unit[1:]
    return unit


def convert_units(
    value: Union[float, np.ndarray],
    from_unit: str,
    to_unit: str,
    molecular_weight: Optional[float] = None,
) -> Union[float, np.ndarray]:
    """
    Convert between concentration units used in fluorescence assays.

    Molar units (M, mM, µM, nM, pM) convert among themselves, as do mass
    concentrations (g/L, mg/mL, mg/L, µg/mL, µg/L, ng/mL, ng/L, pg/mL).
    Converting between the two families needs ``molecular_weight``. Micro may
    be written with the micro sign, Greek mu or a plain ``u``.

    Parameters
    ----------
    value : float or array-like
        Value(s) to convert
    from_unit : str
        Source unit
    to_unit : str
        Target unit
    molecular_weight : float, optional
        In g/mol. Required only between molar and mass units.

    Returns
    -------
    float or array-like
        Converted value(s)

    Raises
    ------
    ValueError
        For an unknown unit, or for a molar/mass conversion without
        ``molecular_weight``.

    Examples
    --------
    >>> convert_units(1000, "nM", "µM")
    1.0
    >>> convert_units(200, "ng/mL", "ug/mL")
    0.2
    >>> convert_units(1, "uM", "ng/mL", molecular_weight=500)
    500.0
    """
    source, target = _canonical_unit(from_unit), _canonical_unit(to_unit)
    for unit, given in ((source, from_unit), (target, to_unit)):
        if unit not in _MOLAR and unit not in _MASS:
            raise ValueError(
                f"Unknown unit {given!r}. Supported: "
                f"{sorted(_MOLAR)} and {sorted(_MASS)}"
            )

    if source == target:
        return value

    if (source in _MOLAR) == (target in _MOLAR):
        table = _MOLAR if source in _MOLAR else _MASS
        return value * (table[source] / table[target])

    if molecular_weight is None:
        raise ValueError(
            f"Converting {from_unit} to {to_unit} needs molecular_weight (g/mol)"
        )
    # Unit ratio first, so matching prefixes cancel exactly (1e-6 / 1e-6).
    if source in _MOLAR:  # mol/L * g/mol = g/L
        return value * (_MOLAR[source] / _MASS[target]) * molecular_weight
    return value * (_MASS[source] / _MOLAR[target]) / molecular_weight
