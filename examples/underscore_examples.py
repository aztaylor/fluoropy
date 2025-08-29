"""
Example demonstrating how __all__ works.
"""

# In a module called 'example_module.py'
__all__ = ["public_function", "PublicClass"]

def public_function():
    """This is intended for users."""
    return "I'm public!"

def _private_function():
    """This is internal implementation."""
    return "I'm private!"

def __very_private_function():
    """This is really internal."""
    return "I'm very private!"

class PublicClass:
    """This is for users."""
    pass

class _PrivateClass:
    """This is internal."""
    pass

# What happens when users import:

# Method 1: from example_module import *
# Only gets: public_function, PublicClass
# Does NOT get: _private_function, __very_private_function, _PrivateClass

# Method 2: import example_module
# Can access everything:
# example_module.public_function()      ✅ Works
# example_module._private_function()    ⚠️ Works but shouldn't use
# example_module.__very_private_function() ⚠️ Works but really shouldn't use

def demonstrate_all_behavior():
    """Show how __all__ affects imports."""

    # This is what __all__ controls:
    print("Items in __all__:", __all__)

    # All items in the module (including private ones):
    all_items = [name for name in globals() if not name.startswith('__')]
    print("All items in module:", all_items)

    # What gets imported with "from module import *":
    star_import_items = [name for name in globals() if name in __all__]
    print("Available with 'import *':", star_import_items)

if __name__ == "__main__":
    demonstrate_all_behavior()
