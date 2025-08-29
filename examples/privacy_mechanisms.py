"""
Demonstration of Python's privacy mechanisms.
"""

class PrivacyExample:
    """Shows different levels of privacy in Python."""

    def __init__(self):
        # Public attribute
        self.public_data = "Anyone can access this"

        # Protected (convention only)
        self._protected_data = "Please don't access directly"

        # Private (name mangling)
        self.__private_data = "This gets name mangled"

        # Very private method
        self.__setup_internals()

    # PUBLIC METHODS
    def get_data(self):
        """Public method to access data safely."""
        return {
            'public': self.public_data,
            'protected': self._protected_data,
            'private': self.__private_data  # Can access from inside class
        }

    # PROTECTED METHOD (convention)
    def _internal_calculation(self, value):
        """Protected method - please don't call directly."""
        return value * 2 + self.__private_multiplier()

    # PRIVATE METHOD (name mangled)
    def __private_multiplier(self):
        """Private method - name gets mangled."""
        return 42

    def __setup_internals(self):
        """Private setup method."""
        self.__internal_state = "initialized"


def demonstrate_privacy():
    """Show what can and can't be accessed."""

    obj = PrivacyExample()

    print("=== PUBLIC ACCESS ===")
    print("Public data:", obj.public_data)  # ✅ Works fine

    print("\n=== PROTECTED ACCESS (Convention) ===")
    print("Protected data:", obj._protected_data)  # ⚠️ Works but shouldn't
    result = obj._internal_calculation(10)  # ⚠️ Works but discouraged
    print("Protected method result:", result)

    print("\n=== PRIVATE ACCESS ATTEMPTS ===")

    # Try to access private data - FAILS
    try:
        print("Private data:", obj.__private_data)
    except AttributeError as e:
        print(f"❌ Can't access __private_data: {e}")

    # Try to call private method - FAILS
    try:
        obj.__private_multiplier()
    except AttributeError as e:
        print(f"❌ Can't call __private_method: {e}")

    print("\n=== NAME MANGLING REVEALED ===")

    # Show what Python actually did to the names
    print("Object attributes:")
    for attr in dir(obj):
        if 'private' in attr.lower() or 'PrivacyExample' in attr:
            print(f"  {attr}")

    # You CAN still access if you know the mangled name
    mangled_data = obj._PrivacyExample__private_data
    print(f"Mangled access: {mangled_data}")  # ⚠️ Possible but very bad practice

    mangled_result = obj._PrivacyExample__private_multiplier()
    print(f"Mangled method result: {mangled_result}")


if __name__ == "__main__":
    demonstrate_privacy()
