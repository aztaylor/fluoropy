# fluoropy

A Python package supplying tools for experimental design and data analysis for
plate reader based fluorescence assays. While the package is focused on biological
experiments that rely on time series fluorescence readout, It should be
adaptable to other similar experiments.

Currently, the code base include tools for plate reader experiment analysis,
including endpoint analysis, timeseries analysis, data normalization,
statistical analysis, and visualization. I plan to add tools for more experimental pipelines in the future.

## Installation

```bash
pip install fluoropy
```

## Usage

```python
import fluoropy

# Your usage examples here
```

## Development

1. Clone the repository:
```bash
git clone https://github.com/aztaylor/fluoropy.git
cd fluoropy
```

2. Create a virtual environment:
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

3. Install in development mode:
```bash
pip install -e ".[dev]"
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.
