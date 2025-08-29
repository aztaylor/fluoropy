"""
Using @property for controlled access to private data.
"""

class SecurePlate:
    """Plate class with controlled access to internal data."""

    def __init__(self, plate_format="96"):
        # Private attributes using name mangling
        self.__wells = {}
        self.__format = plate_format
        self.__validation_enabled = True
        self.__modification_count = 0

    # READ-ONLY PROPERTIES
    @property
    def format(self):
        """Read-only access to plate format."""
        return self.__format

    @property
    def well_count(self):
        """Read-only access to number of wells."""
        return len(self.__wells)

    @property
    def modification_count(self):
        """Read-only access to how many times plate was modified."""
        return self.__modification_count

    # CONTROLLED READ/WRITE PROPERTY
    @property
    def validation_enabled(self):
        """Get validation status."""
        return self.__validation_enabled

    @validation_enabled.setter
    def validation_enabled(self, value):
        """Set validation status with type checking."""
        if not isinstance(value, bool):
            raise TypeError("validation_enabled must be a boolean")
        self.__validation_enabled = value
        print(f"Validation {'enabled' if value else 'disabled'}")

    # COMPUTED PROPERTY
    @property
    def well_positions(self):
        """Get list of all well positions (computed each time)."""
        return sorted(self.__wells.keys())

    # METHODS USING PRIVATE DATA
    def add_well(self, position, fluorescence=None):
        """Add a well with validation."""
        if self.__validation_enabled:
            self.__validate_position(position)

        self.__wells[position] = fluorescence
        self.__modification_count += 1

    def get_well(self, position):
        """Get well data."""
        return self.__wells.get(position)

    def __validate_position(self, position):
        """Private validation method."""
        if not isinstance(position, str):
            raise TypeError("Position must be a string")
        if len(position) < 2:
            raise ValueError("Position must be at least 2 characters")

    # PREVENT DIRECT ACCESS TO WELLS DICT
    def get_wells_copy(self):
        """Get a copy of wells data (not the original dict)."""
        return self.__wells.copy()

    def __repr__(self):
        return f"SecurePlate({self.format}, {self.well_count} wells)"


def demonstrate_properties():
    """Show how properties provide controlled access."""

    plate = SecurePlate("96")

    print("=== READ-ONLY PROPERTIES ===")
    print(f"Format: {plate.format}")  # ✅ Can read
    print(f"Well count: {plate.well_count}")  # ✅ Can read

    # Try to modify read-only property - FAILS
    try:
        plate.format = "384"
    except AttributeError as e:
        print(f"❌ Can't set read-only property: {e}")

    print("\n=== CONTROLLED WRITE PROPERTY ===")
    print(f"Validation enabled: {plate.validation_enabled}")  # ✅ Can read

    plate.validation_enabled = False  # ✅ Can write with validation

    # Try invalid value - FAILS
    try:
        plate.validation_enabled = "yes"
    except TypeError as e:
        print(f"❌ Invalid type rejected: {e}")

    print("\n=== COMPUTED PROPERTIES ===")
    plate.add_well("A1", 1000)
    plate.add_well("B2", 1500)
    print(f"Well positions: {plate.well_positions}")  # ✅ Computed fresh each time
    print(f"Modifications: {plate.modification_count}")  # ✅ Tracks changes

    print("\n=== PRIVATE DATA PROTECTION ===")

    # Can't access private wells dict directly
    try:
        wells = plate.__wells
    except AttributeError as e:
        print(f"❌ Can't access private __wells: {e}")

    # But can get a safe copy
    wells_copy = plate.get_wells_copy()
    print(f"✅ Safe copy: {wells_copy}")

    # Modifying the copy doesn't affect the original
    wells_copy["C3"] = 999
    print(f"Original unchanged: {plate.get_wells_copy()}")


if __name__ == "__main__":
    demonstrate_properties()
