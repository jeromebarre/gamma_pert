# Perturbator - Emission Perturbation Generator

A Python tool for generating spatially-correlated random perturbations for atmospheric emission fields, designed for ensemble data assimilation applications.

## Overview

The **Perturbator** tool (`gamma_pert.py`) creates smooth, spatially-correlated random perturbation fields using gamma distributions and Gaussian kernel smoothing. This is particularly useful for:

- Generating ensemble members for emission uncertainty quantification
- Creating scaling factors for atmospheric chemistry models
- Perturbing emission inventories (e.g., CEDS - Community Emissions Data System)

## Features

- **Gamma-distributed perturbations**: Uses gamma distributions that transition from Gaussian-like (low errors) to lognormal-like (high errors) shapes
- **Spatial correlation**: Applies Gaussian kernel smoothing to create horizontally correlated perturbation fields
- **Multi-sector support**: Apply different perturbation magnitudes and correlation lengths to different emission sectors
- **Multi-species support**: Process multiple chemical species simultaneously
- **Ensemble generation**: Generate multiple ensemble members in a single run
- **NetCDF I/O**: Read and write standard NetCDF emission files with compression support

## Installation

### Dependencies

```bash
pip install numpy scipy xarray pyyaml
```

### Required Python packages:
- `numpy` - Numerical operations
- `scipy` - Signal processing and interpolation
- `xarray` - NetCDF file handling
- `pyyaml` - YAML configuration parsing

## Usage

### Basic Command

```bash
python gamma_pert.py -i input.yaml
```

### Command Line Arguments

| Argument | Description |
|----------|-------------|
| `-i`, `--yaml_file` | **Required.** Path to the YAML configuration file |

## Configuration

The tool is configured via a YAML input file. Below are the configuration parameters:

### YAML Configuration Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `template file` | string | NetCDF file used as a template for grid dimensions |
| `emission files` | list | List of emission NetCDF files to perturb |
| `species list` | list | List of chemical species prefixes (e.g., `['NO', 'CO', 'CH2O']`) |
| `sector list` | list | List of emission sector suffixes (e.g., `['_agr', '_ene', '_ind']`) |
| `sector pert` | list | Relative perturbation magnitude for each sector (standard deviation) |
| `sector hcor` | list | Horizontal correlation length for each sector (in km) |
| `geo dims` | list | Names of longitude and latitude dimensions `[lon_dim, lat_dim]` |
| `time dim` | string | Name of the time dimension |
| `members` | int | Number of ensemble members to generate |
| `scaling factors out` | bool | If `True`, output scaling factor fields as separate variables |
| `only scaling` | bool | If `True`, only output scaling factors (drop original perturbed fields) |
| `inpath` | string | Input directory path |
| `outpath` | string | Output directory path |
| `domain` | string | Domain type (currently only `"global"` is supported) |

### Example YAML Configuration

```yaml
template file: NO_anthro_CEDS_x3600_y1800_t12.2019.nc
emission files: [NO_anthro_CEDS_x3600_y1800_t12.2019.nc, CO_anthro_CEDS_x3600_y1800_t12.2019.nc]
species list: ['NO', 'CO']
sector list: [_agr, _ene, _ind, _rco, _shp, _slv, _tra, _wst]
sector pert: [0.75, 0.3, 0.3, 0.45, 0.40, 0.5, 0.35, 0.55]  # relative error (not %)
sector hcor: [500, 500, 500, 500, 1000, 500, 500, 500]       # correlation length in km
geo dims: [lon, lat]
time dim: time
members: 32
scaling factors out: True
only scaling: True
inpath: /path/to/input/
outpath: /path/to/output/
domain: global
```

## Algorithm Details

### 1. Gamma Distribution Sampling (`gamma_scaled`)

Samples from a gamma distribution with:
- **Mean = 1** (multiplicative scaling factor)
- **Standard deviation = err** (specified perturbation magnitude)

The gamma distribution naturally produces positive values, making it suitable for emission scaling factors. The shape transitions from:
- **Gaussian-like** for low errors (< 0.3)
- **Lognormal-like** for high errors (> 1.0)

### 2. Random Sample Generation (`BuildRandomSample`)

1. Creates a coarse random field at reduced resolution (proportional to kernel width)
2. Samples from the scaled gamma distribution
3. Upsamples to full resolution using bilinear interpolation

### 3. Gaussian Kernel Smoothing (`Kernel`, `Smooth`)

1. Generates an n-dimensional Gaussian kernel with specified standard deviation
2. Normalizes the kernel to preserve the mean
3. Applies FFT-based convolution for efficient smoothing

### 4. Member Perturbation (`BuildMemberPert`)

Combines the random sampling and smoothing steps:
1. Builds a random sample at appropriate resolution
2. Applies Gaussian smoothing with the specified correlation length
3. Returns the final perturbation field

## Output Files

Output files follow the naming convention:
```
{base_name}_pert{XXX}.nc
```
Where `XXX` is the zero-padded ensemble member number (001, 002, etc.).

### Output Variables

For each emission variable `{species}{sector}`:
- **Perturbed emissions**: `{species}{sector}` (if `only scaling: False`)
- **Scaling factors**: `{species}{sector}_pert` (if `scaling factors out: True`)

### Compression

Output NetCDF files are compressed with:
- **zlib compression** (complevel=5)
- **Shuffle filter** enabled
- **Chunking** optimized for typical access patterns

## Mathematical Background

### Perturbation Statistics

For a perturbation with relative error $\sigma$:

- **Shape parameter**: $\alpha = 1/\sigma^2$
- **Scale parameter**: $\theta = \sigma^2$
- **Mean**: $\mu = \alpha \cdot \theta = 1$
- **Variance**: $\text{Var} = \alpha \cdot \theta^2 = \sigma^2$

### Spatial Correlation

The Gaussian smoothing kernel creates perturbations with approximate horizontal correlation length $L_h$:

$$L_h \approx 2\sigma_{\text{kernel}}$$

Where $\sigma_{\text{kernel}}$ is the standard deviation of the Gaussian kernel in grid points.

## Examples

### Example 1: Basic NOx Perturbations

```yaml
template file: NOx_emissions.nc
emission files: [NOx_emissions.nc]
species list: ['NOx']
sector list: [_agr, _ene, _tra]
sector pert: [0.5, 0.3, 0.4]
sector hcor: [200, 300, 150]
geo dims: [lon, lat]
time dim: time
members: 10
scaling factors out: True
only scaling: False
inpath: ./data/
outpath: ./output/
domain: global
```

### Example 2: Multi-Species CEDS Emissions

See `input_CEDS_NO.yaml` and `input_CEDS_NO_hourly.yaml` for complete examples with multiple species and all CEDS sectors.

## Emission Sectors

The tool supports standard CEDS emission sectors:

| Sector Code | Description |
|-------------|-------------|
| `_agr` | Agriculture |
| `_ene` | Energy |
| `_ind` | Industry |
| `_rco` | Residential and Commercial |
| `_shp` | Shipping |
| `_slv` | Solvents |
| `_tra` | Transportation |
| `_wst` | Waste |

## Limitations

- Currently only supports **global** domain (limited area domains not yet implemented)
- Assumes **regular lat-lon grids**
- Time dimension must be consistent across all input files

## License

This project is provided as-is for research purposes.

## Author

Jerome Barre

## References

- CEDS (Community Emissions Data System): https://github.com/JGCRI/CEDS
- Gamma distribution for emission uncertainties: Various atmospheric chemistry data assimilation literature
