"""
Plate and Well classes for managing fluorescence assay data.
"""

from typing import Dict, List, Optional, Tuple, Union, Any
import numpy as np

# ---------------------------------------------------------------------------
# Well roles
# ---------------------------------------------------------------------------
# A well's role is stored once, as a string, and the is_* booleans are derived
# from it. Storing the booleans independently would allow states like
# "negative control that is not a control", and would require callers to keep
# several fields in sync by hand.
#
# NOTE ON POLARITY: "negative control" means the *no-effect* reference and
# "positive control" the *maximal-effect* reference. Neither says anything
# about signal direction. For a repressing construct the negative control (a
# non-targeting guide, say) carries the HIGHEST signal. Nothing in this package
# assumes an ordering between them.

ROLE_SAMPLE = "sample"
ROLE_BLANK = "blank"
ROLE_CONTROL = "control"  # a control whose polarity is unspecified
ROLE_NEGATIVE_CONTROL = "negative_control"
ROLE_POSITIVE_CONTROL = "positive_control"

VALID_ROLES = frozenset({
    ROLE_SAMPLE,
    ROLE_BLANK,
    ROLE_CONTROL,
    ROLE_NEGATIVE_CONTROL,
    ROLE_POSITIVE_CONTROL,
})

#: Accepted spellings that canonicalize to a role in VALID_ROLES.
ROLE_ALIASES = {
    "nc": ROLE_NEGATIVE_CONTROL,
    "neg": ROLE_NEGATIVE_CONTROL,
    "neg_control": ROLE_NEGATIVE_CONTROL,
    "no_effect": ROLE_NEGATIVE_CONTROL,
    "pc": ROLE_POSITIVE_CONTROL,
    "pos": ROLE_POSITIVE_CONTROL,
    "pos_control": ROLE_POSITIVE_CONTROL,
    "max_effect": ROLE_POSITIVE_CONTROL,
    "ctrl": ROLE_CONTROL,
    "test": ROLE_SAMPLE,
    "unknown": ROLE_SAMPLE,
}

#: Roles that count as controls for get_control_wells() and is_control.
CONTROL_ROLES = frozenset({
    ROLE_CONTROL,
    ROLE_NEGATIVE_CONTROL,
    ROLE_POSITIVE_CONTROL,
})


def canonical_role(role: Optional[str]) -> str:
    """
    Normalize a role string, accepting the aliases in ROLE_ALIASES.

    Raises
    ------
    ValueError
        If the role is not recognised. Roles are a closed set on purpose: a
        typo silently creating a new role would quietly drop wells out of
        every control lookup.
    """
    if role is None:
        return ROLE_SAMPLE

    key = str(role).strip().lower().replace("-", "_").replace(" ", "_")

    if key in VALID_ROLES:
        return key
    if key in ROLE_ALIASES:
        return ROLE_ALIASES[key]

    raise ValueError(
        f"Unknown well role {role!r}. Valid roles: {sorted(VALID_ROLES)}; "
        f"accepted aliases: {sorted(ROLE_ALIASES)}"
    )


class Well:
    """
    Represents a single well in a microplate - a simple data container.

    This class stores raw data and metadata for a single well without
    performing any statistical calculations. Statistical analysis is
    handled by the Sample and SampleFrame classes.

    Features:
    - Stores well position and identification
    - Stores sample information and metadata
    - Stores raw time series data for multiple measurement types
    - Supports well exclusion from analysis
    - Simple data container - no statistical methods
    """

    def __init__(self, well_id: str, row: int, column: int):
        """
        Initialize a Well object.

        Parameters
        ----------
        well_id : str
            Well identifier (e.g., 'A1', 'B2')
        row : int
            Row index (0-based, e.g., 0=A, 1=B, etc.)
        column : int
            Column index (0-based, e.g., 0=1, 1=2, etc.)
        """
        # Position information
        self.well_id = well_id
        self.row = row  # 0-based row index
        self.column = column  # 0-based column index
        self.row_letter = chr(ord('A') + row)  # A, B, C, etc.
        self.column_number = column + 1  # 1, 2, 3, etc.

        # Alternative access for backward compatibility
        self.position = well_id

        # Sample information. sample_name identifies which sample this well
        # holds; `sample_type` remains available as an alias.
        self.sample_name: Optional[str] = None
        self.concentration: Optional[float] = None
        self.medium: Optional[str] = None

        # Molecule concentrations (without units in names)
        self.antibiotics: Dict[str, float] = {}  # e.g., {'Kan': 50.0, 'Chlor': 34.0}
        self.inducers: Dict[str, float] = {}  # e.g., {'aTc': 200.0, 'IPTG': 0.5}
        self.other_modifications: Dict[str, float] = {}  # e.g., {'supplement': 1.0}

        # Units for each molecule
        self.antibiotics_units: Dict[str, str] = {}  # e.g., {'Kan': 'µg/mL', 'Chlor': 'µg/mL'}
        self.inducers_units: Dict[str, str] = {}  # e.g., {'aTc': 'ng/mL', 'IPTG': 'mM'}
        self.other_modifications_units: Dict[str, str] = {}  # e.g., {'supplement': 'g/L'}

        # Strain modifications (non-chemical)
        self.strain_modifications: Optional[List[str]] = None

        # Molecule of interest (which molecule's concentration is "the" concentration)
        self.moi: Optional[str] = None

        # Well classification. Stored once as a role; the is_* booleans below
        # are derived, so they cannot disagree with each other.
        self._role: str = ROLE_SAMPLE

        # Exclusion system
        self.exclude: bool = False
        self.exclusion_reason: Optional[str] = None

        # Raw data storage
        self.time_series: Dict[str, np.ndarray] = {}  # Raw time series data
        self.time_points: Optional[np.ndarray] = None

        self.plate_id: Optional[str] = None  # To be set when added to a plate

        # Metadata storage
        self.metadata: Dict[str, Any] = {}

    def __repr__(self) -> str:
        """String representation of the well."""
        excluded_str = " [EXCLUDED]" if self.exclude else ""
        role_str = "" if self.role == ROLE_SAMPLE else f", role={self.role}"
        return (f"Well({self.well_id}, sample={self.sample_name}, "
                f"conc={self.concentration}{role_str}){excluded_str}")

    # ======================================================================
    # IDENTITY
    # ======================================================================

    @property
    def sample_type(self) -> Optional[str]:
        """Alias of :attr:`sample_name`, kept for existing code and notebooks."""
        return self.sample_name

    @sample_type.setter
    def sample_type(self, value: Optional[str]) -> None:
        self.sample_name = value

    # ======================================================================
    # ROLE
    # ======================================================================
    # One stored field, several derived views. Assigning any of the booleans
    # writes through to the role, so existing code that sets well.is_blank
    # keeps working and the fields cannot drift apart.

    @property
    def role(self) -> str:
        """
        What this well is for: sample, blank, control, negative_control or
        positive_control.

        "negative" and "positive" describe the *effect* (none vs maximal), not
        the signal direction -- a repressing construct's negative control
        carries the highest signal.

        Accepts the aliases in :data:`ROLE_ALIASES` on assignment (``'nc'``,
        ``'no_effect'``, ``'max_effect'``, ...) and stores the canonical form.
        """
        return getattr(self, "_role", ROLE_SAMPLE)

    @role.setter
    def role(self, value: Optional[str]) -> None:
        self._role = canonical_role(value)

    def _set_role_flag(self, flag: bool, role_when_true: str,
                       roles_cleared: frozenset) -> None:
        """
        Write-through helper for the boolean views.

        Setting a flag True adopts its role. Setting it False only resets to
        `sample` when the current role is one the flag actually covers, so
        `well.is_blank = False` never silently clears an unrelated role.
        """
        if flag:
            self._role = role_when_true
        elif self.role in roles_cleared:
            self._role = ROLE_SAMPLE

    @property
    def is_blank(self) -> bool:
        return self.role == ROLE_BLANK

    @is_blank.setter
    def is_blank(self, value: bool) -> None:
        self._set_role_flag(bool(value), ROLE_BLANK, frozenset({ROLE_BLANK}))

    @property
    def is_control(self) -> bool:
        """True for any control, whatever its polarity."""
        return self.role in CONTROL_ROLES

    @is_control.setter
    def is_control(self, value: bool) -> None:
        # Assigning True keeps an already-known polarity rather than
        # flattening negative_control back to a bare control.
        if value and self.role in CONTROL_ROLES:
            return
        self._set_role_flag(bool(value), ROLE_CONTROL, CONTROL_ROLES)

    @property
    def is_negative_control(self) -> bool:
        """The no-effect reference. Says nothing about signal direction."""
        return self.role == ROLE_NEGATIVE_CONTROL

    @is_negative_control.setter
    def is_negative_control(self, value: bool) -> None:
        self._set_role_flag(bool(value), ROLE_NEGATIVE_CONTROL,
                            frozenset({ROLE_NEGATIVE_CONTROL}))

    @property
    def is_positive_control(self) -> bool:
        """The maximal-effect reference. Says nothing about signal direction."""
        return self.role == ROLE_POSITIVE_CONTROL

    @is_positive_control.setter
    def is_positive_control(self, value: bool) -> None:
        self._set_role_flag(bool(value), ROLE_POSITIVE_CONTROL,
                            frozenset({ROLE_POSITIVE_CONTROL}))

    # Short aliases, matching how plates get talked about at the bench.
    @property
    def is_nc(self) -> bool:
        return self.is_negative_control

    @is_nc.setter
    def is_nc(self, value: bool) -> None:
        self.is_negative_control = value

    @property
    def is_pc(self) -> bool:
        return self.is_positive_control

    @is_pc.setter
    def is_pc(self, value: bool) -> None:
        self.is_positive_control = value

    # ======================================================================
    # BASIC INFORMATION METHODS
    # ======================================================================

    def set_sample_info(self, sample_type: str,
                        concentration: Optional[float] = None,
                        medium: Optional[str] = None,
                        antibiotics: Optional[Dict[str, float]] = None,
                        inducers: Optional[Dict[str, float]] = None,
                        moi: Optional[str] = None,
                        other_modifications: Optional[Dict[str, float]] = None,
                        strain_modifications: Optional[List[str]] = None,
                        antibiotics_units: Optional[Dict[str, str]] = None,
                        inducers_units: Optional[Dict[str, str]] = None,
                        other_modifications_units: Optional[Dict[str, str]] = None,
                        is_blank: bool = False,
                        is_control: bool = False,
                        role: Optional[str] = None):

        """
        Set sample information for the well.

        Parameters
        ----------
        sample_type : str
            Name of the sample this well holds (alias of ``sample_name``)
        concentration : float, optional
            Concentration of the sample. Mutually exclusive with ``moi``.
        role : str, optional
            Well role: ``'sample'``, ``'blank'``, ``'control'``,
            ``'negative_control'`` or ``'positive_control'``, or any alias in
            ``ROLE_ALIASES`` (``'nc'``, ``'no_effect'``, ``'max_effect'``, ...).
            Takes precedence over ``is_blank`` / ``is_control``, which remain
            for existing callers.
        medium : str, optional
            Growth medium used
        antibiotics : Dict[str, float], optional
            Antibiotic name to concentration mapping (e.g., {'Kan': 50.0, 'Chlor': 34.0})
        inducers : Dict[str, float], optional
            Inducer name to concentration mapping (e.g., {'aTc': 200.0, 'IPTG': 0.5})
        moi : str, optional
            Molecule of interest for this well (e.g., 'aTc') - determines primary concentration
        other_modifications : Dict[str, float], optional
            Other modifications with concentrations (e.g., {'supplement': 1.0})
        strain_modifications : List[str], optional
            List of strain modifications
        antibiotics_units : Dict[str, str], optional
            Units for each antibiotic (e.g., {'Kan': 'µg/mL', 'Chlor': 'µg/mL'})
        inducers_units : Dict[str, str], optional
            Units for each inducer (e.g., {'aTc': 'ng/mL', 'IPTG': 'mM'})
        other_modifications_units : Dict[str, str], optional
            Units for each modification (e.g., {'supplement': 'g/L'})
        is_blank : bool, default False
            Whether this well is a blank control
        is_control : bool, default False
            Whether this well is a control
        """
        self.sample_name = sample_type
        self.medium = medium
        self.moi = moi
        self.strain_modifications = strain_modifications

        if role is not None:
            self.role = role
        else:
            # Legacy flags. is_blank first so an explicit is_control still
            # lands, and so passing both does not depend on argument order.
            self.is_blank = is_blank
            self.is_control = is_control

        # Set molecule concentrations
        if antibiotics is not None:
            if not isinstance(antibiotics, dict):
                raise ValueError(f"Antibiotics should be a dict. Got: {type(antibiotics)}")
            self.antibiotics.update(antibiotics)

        if inducers is not None:
            if not isinstance(inducers, dict):
                raise ValueError(f"Inducers should be a dict. Got: {type(inducers)}")
            self.inducers.update(inducers)

        if other_modifications is not None:
            if not isinstance(other_modifications, dict):
                raise ValueError(f"Other modifications should be a dict. Got: {type(other_modifications)}")
            self.other_modifications.update(other_modifications)

        # Set molecule units
        if antibiotics_units is not None:
            if not isinstance(antibiotics_units, dict):
                raise ValueError(f"Antibiotics units should be a dict. Got: {type(antibiotics_units)}")
            self.antibiotics_units.update(antibiotics_units)

        if inducers_units is not None:
            if not isinstance(inducers_units, dict):
                raise ValueError(f"Inducers units should be a dict. Got: {type(inducers_units)}")
            self.inducers_units.update(inducers_units)

        if other_modifications_units is not None:
            if not isinstance(other_modifications_units, dict):
                raise ValueError(f"Other modifications units should be a dict. Got: {type(other_modifications_units)}")
            self.other_modifications_units.update(other_modifications_units)

        self._set_concentration(concentration)

    def _set_concentration(self, concentration: Optional[float] = None):
        """
        Resolve this well's concentration from exactly one source.

        Priority: an explicit ``concentration`` argument, then ``moi`` (look up
        that molecule's value), then a lone inducer/antibiotic/modification.
        Only one source is used, so the value is never ambiguous.

        Parameters
        ----------
        concentration : float, optional
            Explicit concentration. Mutually exclusive with ``self.moi``.

        Raises
        ------
        ValueError
            If both an explicit concentration and an ``moi`` are given, or if
            ``moi`` names a molecule this well does not have.
        """
        if concentration is not None and self.moi is not None:
            raise ValueError(
                "Cannot provide both 'concentration' and 'moi' parameters. "
                "Please provide only one to avoid ambiguity."
            )

        if concentration is not None:
            # Previously this branch existed but assigned self.concentration to
            # itself, so a concentration passed to set_sample_info() was
            # silently discarded.
            self.concentration = float(concentration)
        elif self.concentration is not None and self.moi is None:
            # Already set directly on the attribute; leave it alone.
            pass
        elif self.moi is not None:
            if self.moi in self.inducers.keys():
                self.concentration = self.inducers[self.moi]
            elif self.moi in self.antibiotics.keys():
                self.concentration = self.antibiotics[self.moi]
            elif self.moi in self.other_modifications.keys():
                self.concentration = self.other_modifications[self.moi]
            else:
                raise ValueError(f"MOI '{self.moi}' not found in inducers, antibiotics, or other modifications for this well.\n Current inducers: {self.inducers}\n Current antibiotics: {self.antibiotics}\n Current modifications: {self.other_modifications}")
        elif len(self.inducers) == 1:
            self.concentration = self.inducers[list(self.inducers.keys())[0]]
        elif len(self.antibiotics) == 1:
            self.concentration = self.antibiotics[list(self.antibiotics.keys())[0]]
        elif len(self.other_modifications) == 1:
            self.concentration = self.other_modifications[list(self.other_modifications.keys())[0]]


    def set_concentration_molecule(self, molecule: str):
        """Set the concentration for this well based on the molecule of
        interest. Can be members of antibiotics, inducers, or other
        modifications which have associated concentration values.

        Parameters
        ----------
        molecule : str
            Name of the molecule to set concentration for (e.g., 'aTc', 'Kan')
        """
        if molecule in self.inducers:
            self.concentration = self.inducers[molecule]
        elif molecule in self.antibiotics:
            self.concentration = self.antibiotics[molecule]
        elif molecule in self.other_modifications:
            self.concentration = self.other_modifications[molecule]
        else:
            raise ValueError(
                f"Molecule '{molecule}' not found in any molecule dictionary.\n"
                f"  Inducers: {self.inducers}\n"
                f"  Antibiotics: {self.antibiotics}\n"
                f"  Other modifications: {self.other_modifications}"
            )
        self.moi = molecule

    def get_concentration(self) -> Optional[float]:
        """Get the concentration for this well."""
        return self.concentration

    @property
    def condition_key(self) -> tuple:
        """Return a hashable key representing this well's experimental condition.

        Returns (medium, frozenset(antibiotics), frozenset(inducers), frozenset(other_modifications)).
        Used for matching blanks to samples in SampleFrame.
        """
        return (
            self.medium,
            frozenset(self.antibiotics.items()) if self.antibiotics else frozenset(),
            frozenset(self.inducers.items()) if self.inducers else frozenset(),
            frozenset(self.other_modifications.items()) if self.other_modifications else frozenset()
        )

    # ======================================================================
    # EXCLUSION METHODS
    # ======================================================================

    def exclude_well(self, reason: str = "Manual exclusion"):
        """
        Exclude this well from analysis.

        Parameters
        ----------
        reason : str, default "Manual exclusion"
            Reason for excluding the well
        """
        self.exclude = True
        self.exclusion_reason = reason

    def include_well(self):
        """Include this well back in analysis."""
        self.exclude = False
        self.exclusion_reason = None

    def is_excluded(self) -> bool:
        """Check if this well is excluded from analysis."""
        return self.exclude

    # ======================================================================
    # DATA STORAGE METHODS
    # ======================================================================

    def add_time_series(self, measurement_type: str, data: Union[List, np.ndarray],
                       time_points: Optional[Union[List, np.ndarray]] = None):
        """
        Add time series data for a specific measurement type.

        Parameters
        ----------
        measurement_type : str
            Type of measurement (e.g., 'OD600', 'GFP', 'fluorescence')
        data : array-like
            Time series data values
        time_points : array-like, optional
            Time points corresponding to the data
        """
        self.time_series[measurement_type] = np.array(data)
        if time_points is not None:
            self.time_points = np.array(time_points)

    def get_measurement(self, measurement_type: str) -> Optional[np.ndarray]:
        """
        Get time series data for a specific measurement type.

        Parameters
        ----------
        measurement_type : str
            Type of measurement to retrieve

        Returns
        -------
        np.ndarray or None
            Time series data for the measurement type, or None if not found
        """
        return self.time_series.get(measurement_type)

    def get_available_measurements(self) -> List[str]:
        """Get list of available measurement types."""
        return list(self.time_series.keys())

    def add_metadata(self, key: str, value: Any):
        """
        Add metadata to the well.

        Parameters
        ----------
        key : str
            Metadata key
        value : Any
            Metadata value
        """
        self.metadata[key] = value

    def get_metadata(self, key: str, default: Any = None) -> Any:
        """
        Get metadata value.

        Parameters
        ----------
        key : str
            Metadata key
        default : Any, optional
            Default value if key not found

        Returns
        -------
        Any
            Metadata value or default
        """
        return self.metadata.get(key, default)
