"""
Overview of development dependencies and their purposes.
"""

def explain_dev_dependencies():
    """Explain what each development tool does."""
    
    dev_tools = {
        "pytest>=6.0": {
            "purpose": "Testing framework",
            "what_it_does": "Runs your tests, discovers test files, provides fixtures",
            "commands": [
                "pytest",                    # Run all tests
                "pytest tests/",             # Run tests in specific directory  
                "pytest -v",                 # Verbose output
                "pytest --tb=short",         # Short traceback format
            ],
            "example_usage": """
            # In your tests/test_plate.py
            import pytest
            import fluoropy
            
            def test_plate_creation():
                plate = fluoropy.Plate("96")
                assert plate.format == "96"
            
            def test_invalid_format():
                with pytest.raises(ValueError):
                    fluoropy.Plate("invalid")
            """
        },
        
        "pytest-cov": {
            "purpose": "Test coverage measurement",
            "what_it_does": "Shows how much of your code is covered by tests",
            "commands": [
                "pytest --cov=fluoropy",              # Basic coverage
                "pytest --cov=fluoropy --cov-report=html",  # HTML report
                "pytest --cov=fluoropy --cov-report=term-missing",  # Missing lines
            ],
            "example_usage": """
            # Terminal output shows:
            # Name               Stmts   Miss  Cover   Missing
            # fluoropy/__init__.py    10      2    80%   15-16
            # fluoropy/core/plate.py  50      5    90%   45, 67-70
            """
        },
        
        "black": {
            "purpose": "Code formatting",
            "what_it_does": "Automatically formats your Python code consistently",
            "commands": [
                "black .",                   # Format all Python files
                "black fluoropy/",           # Format specific directory
                "black --check .",           # Check if formatting needed (CI)
                "black --diff .",            # Show what would change
            ],
            "example_usage": """
            # Before black:
            def function(x,y,z):
                return x+y*z
            
            # After black:
            def function(x, y, z):
                return x + y * z
            """
        },
        
        "isort": {
            "purpose": "Import sorting",
            "what_it_does": "Organizes your import statements consistently",
            "commands": [
                "isort .",                   # Sort all imports
                "isort fluoropy/",           # Sort specific directory
                "isort --check-only .",      # Check if sorting needed
                "isort --diff .",            # Show what would change
            ],
            "example_usage": """
            # Before isort:
            import os
            import fluoropy
            import numpy as np
            from typing import Dict
            import pandas as pd
            
            # After isort:
            import os
            from typing import Dict
            
            import numpy as np
            import pandas as pd
            
            import fluoropy
            """
        },
        
        "flake8": {
            "purpose": "Code linting",
            "what_it_does": "Finds style issues, unused imports, syntax errors",
            "commands": [
                "flake8 .",                  # Lint all files
                "flake8 fluoropy/",          # Lint specific directory
                "flake8 --statistics .",     # Show error statistics
            ],
            "example_usage": """
            # Flake8 output:
            # fluoropy/core/plate.py:15:80: E501 line too long (85 > 79 characters)
            # fluoropy/core/plate.py:23:1: F401 'numpy' imported but unused
            # fluoropy/core/plate.py:45:5: E303 too many blank lines (3)
            """
        },
        
        "mypy": {
            "purpose": "Static type checking",
            "what_it_does": "Checks type hints and finds type-related errors",
            "commands": [
                "mypy fluoropy/",            # Type check package
                "mypy --strict fluoropy/",   # Strict type checking
                "mypy --install-types",      # Install missing type stubs
            ],
            "example_usage": """
            # Type hints in your code:
            def add_well(self, well: Well) -> None:
                self.wells[well.position] = well
            
            # Mypy catches errors:
            # error: Argument 1 to "add_well" has incompatible type "str"; expected "Well"
            """
        },
        
        "pre-commit": {
            "purpose": "Git hooks automation",
            "what_it_does": "Runs checks automatically before commits",
            "commands": [
                "pre-commit install",        # Install git hooks
                "pre-commit run --all-files", # Run on all files
                "pre-commit autoupdate",     # Update hook versions
            ],
            "example_usage": """
            # .pre-commit-config.yaml
            repos:
            - repo: https://github.com/psf/black
              rev: 23.3.0
              hooks:
              - id: black
            - repo: https://github.com/pycqa/isort
              rev: 5.12.0
              hooks:
              - id: isort
            """
        }
    }
    
    return dev_tools


def explain_docs_dependencies():
    """Explain documentation tools."""
    
    docs_tools = {
        "sphinx>=4.0": {
            "purpose": "Documentation generator",
            "what_it_does": "Creates beautiful documentation websites from your docstrings",
            "commands": [
                "sphinx-quickstart docs",   # Initialize docs
                "sphinx-build docs docs/_build", # Build documentation
                "sphinx-autogen docs/api.rst",   # Generate API docs
            ],
            "example_usage": """
            # Your docstrings become web pages:
            def add_well(self, well: Well) -> None:
                '''
                Add a well to the plate.
                
                Parameters
                ----------
                well : Well
                    The well object to add
                    
                Examples
                --------
                >>> plate = Plate("96")
                >>> well = Well("A1", 1000)
                >>> plate.add_well(well)
                '''
            """
        },
        
        "sphinx-rtd-theme": {
            "purpose": "Documentation theme",
            "what_it_does": "Makes your docs look like ReadTheDocs (professional)",
            "commands": [
                # Added to conf.py: html_theme = 'sphinx_rtd_theme'
            ],
            "example_usage": "Provides the clean, professional look you see on readthedocs.org"
        },
        
        "myst-parser": {
            "purpose": "Markdown support in Sphinx",
            "what_it_does": "Let you write documentation in Markdown instead of just reStructuredText",
            "commands": [
                # Added to conf.py: extensions = ['myst_parser']
            ],
            "example_usage": """
            # You can write docs in Markdown:
            # Getting Started
            
            ## Installation
            ```bash
            pip install fluoropy
            ```
            
            ## Quick Example
            ```python
            import fluoropy
            plate = fluoropy.Plate("96")
            ```
            """
        }
    }
    
    return docs_tools


def explain_viz_dependencies():
    """Explain visualization tools."""
    
    viz_tools = {
        "matplotlib>=3.5.0": {
            "purpose": "Basic plotting",
            "what_it_does": "Create static plots, heatmaps, dose-response curves",
            "example_usage": """
            import matplotlib.pyplot as plt
            import numpy as np
            
            # Plot dose-response curve
            concentrations = [0, 1, 5, 10, 50, 100]
            fluorescence = [50, 200, 500, 700, 900, 950]
            
            plt.semilogx(concentrations[1:], fluorescence[1:])
            plt.xlabel('Concentration (µM)')
            plt.ylabel('Fluorescence (RFU)')
            plt.title('Dose-Response Curve')
            """
        },
        
        "plotly>=5.0.0": {
            "purpose": "Interactive plotting",
            "what_it_does": "Create interactive plots, 3D visualizations, web-ready charts",
            "example_usage": """
            import plotly.express as px
            import pandas as pd
            
            # Interactive plate heatmap
            df = plate.get_fluorescence_data()
            fig = px.imshow(plate.get_plate_matrix(), 
                           title="Plate Heatmap",
                           color_continuous_scale="viridis")
            fig.show()  # Opens in browser
            """
        },
        
        "seaborn>=0.11.0": {
            "purpose": "Statistical plotting",
            "what_it_does": "Beautiful statistical plots, distributions, correlations",
            "example_usage": """
            import seaborn as sns
            import pandas as pd
            
            # Plot dose-response with confidence intervals
            df = experiment.get_combined_dataframe()
            sns.lineplot(data=df, x='concentration', y='fluorescence', 
                        estimator='mean', ci=95)
            plt.xscale('log')
            """
        }
    }
    
    return viz_tools


if __name__ == "__main__":
    print("📋 See the source code for detailed explanations of each dependency!")
    print("\n🧪 Dev tools: Testing, formatting, linting, type checking")
    print("📚 Docs tools: Generate beautiful documentation websites") 
    print("📊 Viz tools: Create plots and visualizations")
    
    # Show how to install them
    print("\n💡 To install:")
    print("pip install -e '.[dev]'     # Development tools")
    print("pip install -e '.[docs]'    # Documentation tools") 
    print("pip install -e '.[viz]'     # Visualization tools")
    print("pip install -e '.[dev,docs,viz]'  # Everything!")
