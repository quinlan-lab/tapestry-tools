# tapestry-tools Development Guide

## Project Overview

A bioinformatics tool suite for analyzing DNA methylation patterns, e.g., identifying imprinted loci using founder-phased haplotype data.

## Build/Lint/Test Commands

### Installation (one-time)
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
source $HOME/.local/bin/env
uv tool install git+https://github.com/quinlan-lab/tapestry-tools.git
```

### Development Setup
```bash
uv sync
```

### Running Tools

#### Command-line tools (installed via pyproject.toml):
- `compute-delta-methylation` - Compute methylation differences between haplotypes
- `call-imprinted-loci` - Identify imprinted loci from methylation data
- `compute-methylation-loci` - Compute methylation at specific loci
- `clean-coords` - Convert coordinate formats
- `lift-over` - Lift coordinates between genome builds

#### Example: Run with testing mode (2 samples, faster)
```bash
compute-delta-methylation --sample-meth-beds <file> --delta-meth-bed <output> --testing
```

#### Example: Production run
```bash
compute-delta-methylation --sample-meth-beds <file> --delta-meth-bed <output>
```

### Type Checking & Linting

No explicit ruff configuration exists yet. For type checking and linting:
```bash
# Type checking (mypy not configured, but recommended)
mypy src/

# Linting (if ruff is installed)
ruff check src/
ruff format src/ --check

# Alternatively, use Python's built-in compiler check
python -m py_compile src/tapestry_tools/*.py
```

### Testing

**No formal test suite exists yet.** To verify functionality:
1. Run tools with `--testing` flag (uses 2 samples)
2. Check outputs match expected formats
3. Manual verification of intermediate dataframes

## Code Style Guidelines

### Polars Chaining Style

Always use multi-line chaining for Polars operations:

✅ CORRECT:
```python
import polars as pl

df = (
    pl
    .read_csv("data.csv")
    .select([
        pl.col("col1"),
        pl.col("col2")
    ])
    .filter(pl.col("col1") > 0)
    .groupby("col2")
    .agg(pl.col("col1").mean())
)
```

❌ INCORRECT:
```python
df = pl.read_csv("data.csv").select([pl.col("col1"), pl.col("col2")]).filter(...)
```

**Rules:**
- Start with `pl` on its own line
- Each method call on a new line
- Opening parenthesis on same line as `pl`
- Makes code readable and easier to debug

### Imports

- Standard library imports first (`os`, `sys`, `pathlib`, `argparse`)
- External imports (`polars`, `bioframe`, `rich_argparse`, `tqdm`)
- Relative imports for local project modules only

Example:
```python
import argparse
from rich_argparse import RichHelpFormatter
from pathlib import Path
import polars as pl
from tqdm import tqdm

from .read_data import read_dataframe_from_bed
```

### Type Hints

- Use type hints for all function signatures
- Use `pl.DataFrame` for Polars DataFrames
- Use `Path` from pathlib for file paths
- Example:
```python
def read_tapestry(bed) -> pl.DataFrame:
def write_dataframe_to_bed(df: pl.DataFrame, file_path: str, source: str):
```

### Naming Conventions

- **Functions/Methods:** snake_case (`compute_delta_methylation`, `read_dataframe_from_bed`)
- **Classes:** PascalCase (not currently used extensively)
- **Variables:** snake_case
- **Constants:** UPPER_SNAKE_CASE
- **Private:** prefix with underscore (`_private_function`)

### Error Handling

- Use assertions for internal consistency checks: `assert df.columns == [...]`
- Use meaningful error messages with `print()` or `sys.stderr.write()`
- Exit with non-zero status on errors: `sys.exit(1)`
- Check file existence before reading: `if Path(file).exists()`

### File Format Conventions

**BED Files:**
- Tab-separated values
- Header line starts with `#`
- Comment lines start with `##`
- Use `write_dataframe_to_bed()` for writing

Example header:
```python
"##source='_filename.py'\n#" + "\t".join(columns) + "\n"
```

### Main Function Pattern

All CLI tools follow this structure:
```python
def main():
    class RichDefaultsFormatter(RichHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
        pass
    
    parser = argparse.ArgumentParser(
        description="Tool description",
        formatter_class=RichDefaultsFormatter
    )
    # Add arguments...
    
    if len(sys.argv) == 1:
        parser.print_help(sys.stdout)
        sys.exit(1)
    
    args = parser.parse_args()
    # Execute logic...
```

### Testing Mode

Many functions accept a `testing` parameter that limits processing to 2 samples for faster iteration. Always use `--testing` flag during development.

## Module Structure

**Core modules in `src/tapestry_tools/`:**
- `read_data.py` - File reading utilities (BED, tapestry format)
- `write_data.py` - File writing utilities
- `methylation.py` - Methylation computation logic
- `imprinting.py` - Imprinting identification logic
- `version_sort.py` - Chromosome-aware sorting
- `tile.py` - Genome tile generation
- `liftover.py` - Coordinate conversion between genome builds

**Experiments directory:** Jupyter notebooks and exploration scripts

## Important Notes

- Requires Python >=3.13
- Uses `uv` for dependency management
- Chief dependencies: polars, bioframe, rich-argparse, tqdm
- Reference genome: hg38
- Tile size default: 1000 bp
