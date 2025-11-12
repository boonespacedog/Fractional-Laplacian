"""
📄 File: test_structure.py
Purpose: Structure tests for E53 - verify code correctness
Created: November 11, 2025
Used by: pytest

Discovery-based testing: We DISCOVER properties, not assert hardcoded values
"""

import pytest
import numpy as np
import sys
from pathlib import Path

# Add src to path
sys.path.append(str(Path(__file__).parent.parent))

from src.core import SphericalLaplacian
from src.compute import (
    compute_exact_spectrum_sphere,
    compute_fractional_spectrum_exact,
    compute_curvature_correction,
    compute_weak_curvature_parameter
)


class TestSphericalLaplacianStructure:
    """Test structural properties of SphericalLaplacian class."""

    def test_initialization(self):
        """Test object initialization with various parameters."""
        print("\n" + "="*60)
        print("TESTING: SphericalLaplacian Initialization")
        print("="*60)

        # Valid initialization
        laplacian = SphericalLaplacian(R=2.0, alpha=0.5, l_max=50)

        print(f"✓ Created laplacian with R={laplacian.R}, α={laplacian.alpha}, l_max={laplacian.l_max}")

        # Test parameter storage
        assert laplacian.R == 2.0, "R not stored correctly"
        assert laplacian.alpha == 0.5, "alpha not stored correctly"
        assert laplacian.l_max == 50, "l_max not stored correctly"

        # Test invalid parameters
        with pytest.raises(ValueError):
            SphericalLaplacian(R=-1.0)  # Negative radius
            print("✓ Correctly rejected negative radius")

        with pytest.raises(ValueError):
            SphericalLaplacian(alpha=1.5)  # alpha out of range
            print("✓ Correctly rejected α outside (0,1)")

        print("\n✅ All initialization tests passed")

    def test_eigenvalue_computation(self):
        """Discover eigenvalue formula structure."""
        print("\n" + "="*60)
        print("DISCOVERY: Eigenvalue Formula Structure")
        print("="*60)

        laplacian = SphericalLaplacian(R=1.0, alpha=0.5, l_max=10)

        # Test exact eigenvalues
        for l in [0, 1, 2, 5, 10]:
            lambda_l = laplacian.eigenvalue_exact(l)
            expected_formula = l * (l + 1) / laplacian.R**2

            print(f"ℓ={l:2d}: λ_ℓ = {lambda_l:.6f}")
            print(f"      Formula: ℓ(ℓ+1)/R² = {l}×{l+1}/{laplacian.R}² = {expected_formula:.6f}")

            # Structure test: formula matches computation
            assert abs(lambda_l - expected_formula) < 1e-10, f"Formula mismatch for l={l}"

        print("\n✅ Eigenvalue formula structure confirmed: λ_ℓ = ℓ(ℓ+1)/R²")

    def test_fractional_power(self):
        """Discover fractional eigenvalue behavior."""
        print("\n" + "="*60)
        print("DISCOVERY: Fractional Eigenvalue Behavior")
        print("="*60)

        for alpha in [0.1, 0.5, 0.9]:
            laplacian = SphericalLaplacian(R=1.0, alpha=alpha, l_max=5)
            print(f"\nα = {alpha}:")
            print("-" * 40)

            for l in [1, 2, 3]:
                exact = laplacian.eigenvalue_exact(l)
                fractional = laplacian.eigenvalue_fractional_exact(l)
                expected = exact ** alpha

                print(f"  ℓ={l}: λ_ℓ = {exact:.3f}, λ_ℓ^α = {fractional:.6f}")
                print(f"       Check: ({exact:.3f})^{alpha} = {expected:.6f}")

                # Structure test: power law applied correctly
                assert abs(fractional - expected) < 1e-10, f"Power law error for α={alpha}, l={l}"

        print("\n✅ Fractional power structure confirmed: λ_ℓ^α = [ℓ(ℓ+1)/R²]^α")

    def test_correction_term_structure(self):
        """Discover correction term formula structure."""
        print("\n" + "="*60)
        print("DISCOVERY: Correction Term Formula")
        print("="*60)

        laplacian = SphericalLaplacian(R=2.0, alpha=0.7, l_max=10)

        print(f"Parameters: R={laplacian.R}, α={laplacian.alpha}")
        print("\nCorrection formula: C(ℓ,R,α) = α(α-1)/12R² · ℓ(ℓ+1)")

        # Calculate prefactor
        prefactor = laplacian.alpha * (laplacian.alpha - 1) / 12.0
        print(f"Prefactor: α(α-1)/12 = {laplacian.alpha}×({laplacian.alpha}-1)/12 = {prefactor:.6f}")

        # Note: negative for 0 < α < 1
        if 0 < laplacian.alpha < 1:
            print(f"  → Negative (as expected for 0<α<1)")

        print("\nCorrection values:")
        for l in [0, 1, 2, 5, 10]:
            correction = laplacian.compute_correction_term(l)
            expected = prefactor * l * (l + 1) / laplacian.R**2

            print(f"  ℓ={l:2d}: C = {correction:.8f}")

            # Structure test: formula implementation
            assert abs(correction - expected) < 1e-10, f"Formula mismatch for l={l}"

        print("\n✅ Correction formula structure confirmed")

    def test_array_properties(self):
        """Test array computation and caching."""
        print("\n" + "="*60)
        print("TESTING: Array Properties and Caching")
        print("="*60)

        laplacian = SphericalLaplacian(R=1.5, alpha=0.6, l_max=20)

        # Test array shapes
        exact_arr = laplacian.eigenvalues_exact
        fractional_arr = laplacian.eigenvalues_fractional_exact
        corrected_arr = laplacian.eigenvalues_fractional_corrected
        corrections_arr = laplacian.correction_terms

        print(f"Array shapes (l_max={laplacian.l_max}):")
        print(f"  Exact eigenvalues:      {exact_arr.shape} = ({laplacian.l_max + 1},)")
        print(f"  Fractional exact:       {fractional_arr.shape}")
        print(f"  Fractional corrected:   {corrected_arr.shape}")
        print(f"  Correction terms:       {corrections_arr.shape}")

        # All should have same shape
        assert exact_arr.shape == (laplacian.l_max + 1,), "Wrong shape for exact eigenvalues"
        assert fractional_arr.shape == exact_arr.shape, "Shape mismatch"
        assert corrected_arr.shape == exact_arr.shape, "Shape mismatch"
        assert corrections_arr.shape == exact_arr.shape, "Shape mismatch"

        # Test caching (should return same object)
        exact_arr2 = laplacian.eigenvalues_exact
        assert exact_arr is exact_arr2, "Caching not working"
        print("\n✓ Caching working correctly")

        print("\n✅ All array property tests passed")


class TestComputeFunctions:
    """Test compute module functions."""

    def test_exact_spectrum_computation(self):
        """Test exact eigenvalue spectrum computation."""
        print("\n" + "="*60)
        print("TESTING: Exact Spectrum Computation")
        print("="*60)

        R = 2.5
        l_max = 10
        spectrum = compute_exact_spectrum_sphere(R, l_max)

        print(f"Computed spectrum for R={R}, l_max={l_max}:")
        print(f"  Shape: {spectrum.shape}")
        print(f"  First 5 values: {spectrum[:5]}")

        # Structure tests
        assert spectrum.shape == (l_max + 1,), "Wrong spectrum shape"
        assert spectrum[0] == 0, "λ_0 should be 0"

        # Check monotonicity
        differences = np.diff(spectrum)
        assert np.all(differences >= 0), "Eigenvalues should be non-decreasing"
        print("✓ Eigenvalues are non-decreasing")

        # Check formula
        for l in range(min(5, l_max + 1)):
            expected = l * (l + 1) / R**2
            assert abs(spectrum[l] - expected) < 1e-10, f"Formula error at l={l}"
        print("✓ Formula verified for first 5 values")

        print("\n✅ Exact spectrum computation tests passed")

    def test_curvature_parameter(self):
        """Discover curvature parameter regimes."""
        print("\n" + "="*60)
        print("DISCOVERY: Curvature Parameter Regimes")
        print("="*60)

        R_values = [0.5, 1.0, 2.0, 5.0, 10.0]
        l_test = 10

        print(f"Testing ℓ={l_test} with various radii:")
        print("-" * 50)
        print("   R     κ_α    Regime")
        print("-" * 50)

        for R in R_values:
            kappa = compute_weak_curvature_parameter(R, l_test)

            # Determine regime
            if kappa < 0.1:
                regime = "WEAK     (κ < 0.1)"
            elif kappa <= 1.0:
                regime = "INTERMEDIATE"
            else:
                regime = "STRONG   (κ > 1.0)"

            print(f"{R:5.1f}  {kappa:8.4f}  {regime}")

        print("-" * 50)
        print("Formula: κ_α = ℓ(ℓ+1)/R²")
        print(f"For ℓ={l_test}: κ_α = {l_test*(l_test+1)}/R²")

        # Discover transition radii
        R_weak_boundary = np.sqrt(l_test * (l_test + 1) / 0.1)
        R_strong_boundary = np.sqrt(l_test * (l_test + 1) / 1.0)

        print(f"\nTransition radii for ℓ={l_test}:")
        print(f"  Weak → Intermediate:   R = {R_weak_boundary:.2f}")
        print(f"  Intermediate → Strong: R = {R_strong_boundary:.2f}")

        print("\n✅ Curvature parameter structure discovered")

    def test_fractional_spectrum_handling(self):
        """Test fractional spectrum with edge cases."""
        print("\n" + "="*60)
        print("TESTING: Fractional Spectrum Edge Cases")
        print("="*60)

        eigenvalues = np.array([0.0, 1.0, 2.0, 5.0, 10.0])
        alpha = 0.5

        fractional = compute_fractional_spectrum_exact(eigenvalues, alpha)

        print(f"Input eigenvalues: {eigenvalues}")
        print(f"Fractional order α = {alpha}")
        print(f"Output λ^α: {fractional}")

        # Structure tests
        assert fractional.shape == eigenvalues.shape, "Shape preserved"
        assert fractional[0] == 0.0, "Zero eigenvalue remains zero"
        assert np.all(fractional[1:] > 0), "Positive eigenvalues remain positive"

        # Check power law
        for i in range(1, len(eigenvalues)):
            expected = eigenvalues[i] ** alpha
            assert abs(fractional[i] - expected) < 1e-10, f"Power law error at i={i}"

        print("\n✓ Zero eigenvalue handled correctly")
        print("✓ Power law applied correctly")
        print("\n✅ Fractional spectrum tests passed")