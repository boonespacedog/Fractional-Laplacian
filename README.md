# Fractional Laplacian on Curved Manifolds

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.17585575-blue.svg)](https://doi.org/10.5281/zenodo.17585575)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**Computational verification code for**: "Fractional Laplacian on Curved Manifolds: Coordinate-Invariant Construction and Physical Applications"

**Author**: Oksana Sudoma
**Status**: Publication-ready

---

## Overview

This repository provides computational validation for coordinate-invariant fractional Laplacian operators on curved manifolds, with rigorous eigenvalue analysis on the 2-sphere demonstrating sub-percent accuracy and exact dimensional scaling.

**Key discoveries**:
- **Curvature correction accuracy**: 0.59% error in weak regime (κ_α < 0.1)
- **R^-2 dimensional scaling**: Exact to machine precision (10^-15)
- **Spectral eigenvalue validation**: All 51 correction signs correct
- **Weak curvature threshold**: κ_α < 0.1 boundary empirically confirmed

---

## Experiments

### E53: Sphere Validation (S²)

Validates coordinate-invariant construction on 2-sphere.

**Key results**:
- Coordinate invariance: 10^-10 precision
- Eigenvalue spectrum verified
- Spectral discretization validated
- Runtime: ~2 seconds

**Quick start**:
```bash
cd experiments/E53_sphere_validation
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
python3 main.py
```


---

## Repository Structure

```
Fractional-Laplacian/
├── experiments/
│   └── E53_sphere_validation/
│       ├── src/                 # Core implementation
│       ├── tests/               # Unit tests
│       ├── outputs/             # Results and analytical report
│       ├── main.py              # Entry point
│       └── requirements.txt
├── paper/
│   ├── fractional_laplacian_curved_v6.tex
│   └── fractional_laplacian_curved_v6.pdf
├── LICENSE
├── README.md
└── .gitignore
```

---

## Mathematical Background

**Fractional Laplacian** generalizes classical Laplacian to non-integer orders s ∈ (0,1):

Δˢ = (-Δ)ˢ

On curved manifolds (M,g), coordinate invariance requires spectral construction:

Δˢ f = Σᵢ λᵢˢ ⟨f, φᵢ⟩ φᵢ

where {(λᵢ, φᵢ)} are eigenpairs of the Laplace-Beltrami operator.

**Novel phenomena**:
- **Coordinate invariance** holds on curved spaces
- **Anomalous diffusion**: ⟨r²⟩ ∝ t^α with α < 1 (subdiffusive)
- **Spectral discretization** naturally handles manifold geometry
- **Physical applications**: Diffusion in curved spacetime

---

## Reproducibility

All results are computationally verified:

1. **Run individual experiments**: See Quick Start sections above
2. **Run tests**: `python3 -m pytest tests/ -v` (in each experiment directory)
3. **Expected runtime**: < 10 seconds total (both experiments)

All numerical results match paper claims to stated precision.

---

## Citation

```bibtex
@misc{sudoma2025fractional,
  author = {Sudoma, Oksana},
  title = {Fractional Laplacian on Curved Manifolds: Coordinate-Invariant Construction and Physical Applications},
  year = {2025},
  doi = {10.5281/zenodo.17585575},
  url = {https://github.com/boonespacedog/Fractional-Laplacian}
}
```

---

## License

MIT License - see [LICENSE](LICENSE) file for details.

---

## Author

**Oksana Sudoma** - Independent Researcher

Computational validation and mathematical formalism assisted by Claude (Anthropic). All scientific conclusions and theoretical insights are the author's sole responsibility.

---

## Links

- **Repository**: https://github.com/boonespacedog/Fractional-Laplacian
- **Zenodo Archive**: https://doi.org/10.5281/zenodo.17585575
- **Paper**: See `paper/` directory for latest version (v6)

---

## Physical Applications

The fractional Laplacian on curved manifolds has applications in:

- **Anomalous diffusion** in curved spacetime
- **Non-local field theories** in general relativity
- **Quantum field theory** on curved backgrounds
- **Statistical physics** on non-Euclidean geometries

This work provides the first coordinate-invariant construction with computational validation.
