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
- **Gaussian-like** for low errors (< 30%)
- **Lognormal-like** for high errors (> 100%)

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

### Gamma Distribution Formulation

The gamma distribution is a two-parameter family of continuous probability distributions. The probability density function (PDF) is:

$$f(x; \alpha, \theta) = \frac{x^{\alpha-1} e^{-x/\theta}}{\theta^\alpha \Gamma(\alpha)}, \quad x > 0$$

Where:
- $\alpha > 0$ is the **shape parameter**
- $\theta > 0$ is the **scale parameter**
- $\Gamma(\alpha)$ is the gamma function

### Parameterization for Unit Mean

For emission perturbations, we require scaling factors with **mean = 1** and a specified **relative error** (standard deviation) $\sigma$. The gamma distribution moments are:

$$\mathbb{E}[X] = \alpha \theta$$
$$\text{Var}(X) = \alpha \theta^2$$

To achieve mean = 1 with standard deviation = $\sigma$, we solve:

$$\alpha \theta = 1 \quad \text{(mean constraint)}$$
$$\alpha \theta^2 = \sigma^2 \quad \text{(variance constraint)}$$

which gives:

$$\alpha = \frac{1}{\sigma^2}, \quad \theta = \sigma^2$$

With the derived parameters, the mean is:

$$\mathbb{E}[X] = \alpha \theta = \frac{1}{\sigma^2} \cdot \sigma^2 = 1$$

And the variance is:

$$\text{Var}(X) = \alpha \theta^2 = \frac{1}{\sigma^2} \cdot \sigma^4 = \sigma^2$$

This ensures that:
1. The **average perturbation factor is 1** (no systematic bias in emissions)
2. The **spread is controlled by the specified error** $\sigma$

### Gaussian Approximation (Low Error Regime)

For **small errors** ($\sigma < 0.3$), the shape parameter $\alpha = 1/\sigma^2$ becomes large:

| $\sigma$ | $\alpha$ |
|----------|----------|
| 0.1 | 100 |
| 0.2 | 25 |
| 0.3 | 11.1 |

By the **Central Limit Theorem**, as $\alpha \to \infty$, the gamma distribution converges to a Gaussian:

$$\text{Gamma}(\alpha, \theta) \xrightarrow{\alpha \to \infty} \mathcal{N}(\alpha\theta, \alpha\theta^2) = \mathcal{N}(1, \sigma^2)$$

The **skewness** of the gamma distribution is:

$$\gamma_1 = \frac{2}{\sqrt{\alpha}} = 2\sigma$$

For $\sigma = 0.1$: skewness = 0.2 (nearly symmetric, Gaussian-like)
For $\sigma = 0.3$: skewness = 0.6 (slightly asymmetric)

### Lognormal Approximation (High Error Regime)

For **large errors** ($\sigma > 1.0$), the shape parameter $\alpha = 1/\sigma^2$ becomes small:

| $\sigma$ | $\alpha$ |
|----------|----------|
| 1.0 | 1.0 |
| 1.5 | 0.44 |
| 2.0 | 0.25 |

When $\alpha < 1$, the gamma PDF has an asymptote at $x = 0$ and exhibits heavy right-skewness, similar to a lognormal distribution.

The lognormal distribution $\text{LogNormal}(\mu_{\ln}, \sigma_{\ln})$ has:

$$\mathbb{E}[X] = e^{\mu_{\ln} + \sigma_{\ln}^2/2}$$
$$\text{Var}(X) = (e^{\sigma_{\ln}^2} - 1) e^{2\mu_{\ln} + \sigma_{\ln}^2}$$

For a lognormal with mean = 1 and variance = $\sigma^2$:

$$\mu_{\ln} = -\frac{1}{2}\ln(1 + \sigma^2)$$
$$\sigma_{\ln} = \sqrt{\ln(1 + \sigma^2)}$$

Both gamma and lognormal share key properties in this regime:
- **Strictly positive** values (essential for emission scaling)
- **Heavy right tail** (allows large positive perturbations)
- **Mode < Mean < Median** ordering

### Distribution Shape Summary

| Error Regime | $\sigma$ | $\alpha$ | Skewness | Distribution Shape |
|--------------|----------|----------|----------|-------------------|
| Low | < 0.3 | > 11 | < 0.6 | ≈ Gaussian (symmetric) |
| Medium | 0.3 - 1.0 | 1 - 11 | 0.6 - 2.0 | Intermediate |
| High | > 1.0 | < 1 | > 2.0 | ≈ Lognormal (right-skewed) |

### Why Gamma Over Gaussian or Lognormal?

1. **Always positive**: Unlike Gaussian, gamma never produces negative scaling factors
2. **Flexible shape**: Single parameterization smoothly transitions between Gaussian-like and lognormal-like
3. **Simple mean control**: Easy to ensure mean = 1 for unbiased perturbations
4. **Computational efficiency**: Direct sampling without rejection or transformation

### Spatial Correlation

#### Correlation Length Conversion

The horizontal correlation length is specified in **kilometers** in the YAML configuration (`sector hcor`). The code converts this to grid points using the equatorial Earth circumference:

$$\sigma_{\text{grid}} = \frac{L_h}{R_{\text{eq}}}$$

Where:
- $L_h$ is the correlation length in km
- $R_{\text{eq}} = 40075 / n_{\text{lon}}$ is the approximate grid resolution in km at the equator
- $n_{\text{lon}}$ is the number of longitude grid points

For example, with a 0.5° grid ($n_{\text{lon}} = 720$):
- Grid resolution at equator: $R_{\text{eq}} \approx 55.7$ km
- A 500 km correlation length → $\sigma_{\text{grid}} \approx 9$ grid points

#### Gaussian Kernel Formulation

The Gaussian smoothing kernel in 2D is:

$$K(x, y) = \frac{1}{2\pi\sigma^2} \exp\left(-\frac{x^2 + y^2}{2\sigma^2}\right)$$

The kernel is normalized ($\sum K = 1$) to preserve the mean of the perturbation field during convolution. The effective correlation length is approximately:

$$L_h \approx 2\sigma_{\text{kernel}}$$

#### ⚠️ Polar Singularity Limitation

A kernel with fixed grid-point width represents a much smaller physical distance near the poles. A 500 km correlation at the equator becomes effectively ~250 km at 60°N/S and ~85 km at 80°N/S.

For applications requiring strictly **uniform, isotropic spatial correlations** on the sphere, the **JEDI (Joint Effort for Data assimilation Integration)** framework, for example, provides more sophisticated tools.

These tools could properly account for the spherical geometry and maintain isotropic correlation lengths regardless of latitude. However, they require:
- Full JEDI software stack installation
- Model interface implementation (e.g., FV3, MPAS, etc.)
- More setup complexity

For many practical applications**, especially in the tropics and mid-latitudes where most emissions occur, the simplified lat-lon approach in this tool provides adequate results with minimal setup overhead.

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

## Current Limitations

- Currently only supports **global** domain (limited area domains not yet implemented)
- Assumes **regular lat-lon grids**
- Time dimension must be consistent across all input files

## Author

Jerome Barre

## References

- CEDS (Community Emissions Data System): https://github.com/JGCRI/CEDS
