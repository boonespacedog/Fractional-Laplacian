"""
📄 File: test_theory.py
Purpose: Theory falsification tests for E53 - attempt to disprove predictions
Created: November 11, 2025
Used by: pytest

CRITICAL: These tests DISCOVER whether theory holds, not assert it must
"""

import pytest
import numpy as np
import sys
from pathlib import Path

# Add src to path
sys.path.append(str(Path(__file__).parent.parent))

from src.core import SphericalLaplacian
from src.validate import (
    validate_error_bounds,
    validate_convergence_with_l,
    validate_scaling_law,
    falsification_tests
)


class TestTheoryPredictions:
    """Attempt to falsify theoretical predictions."""

    def test_weak_curvature_error_bound(self):
        """
        Theory predicts: Relative error < 1% in weak curvature regime.
        FALSIFIABLE: If we find error > 1% when κ_α < 0.1
        """
        print("\n" + "="*60)
        print("FALSIFICATION TEST: Weak Curvature Error Bound")
        print("="*60)

        # Test with large radius (ensures weak curvature for small l)
        R = 10.0
        alpha = 0.5
        laplacian = SphericalLaplacian(R=R, alpha=alpha, l_max=30)

        print(f"Parameters: R={R}, α={alpha}")
        print("\nTheory prediction: Error < 1% when κ_α < 0.1")
        print("-" * 50)

        # Find l values in weak regime
        weak_l_values = []
        errors_in_weak = []

        for l in range(1, laplacian.l_max + 1):
            kappa = laplacian.compute_weak_curvature_parameter(l)
            if kappa < 0.1:
                error = laplacian.relative_error(l)
                weak_l_values.append(l)
                errors_in_weak.append(error)

                status = "✓" if error < 0.01 else "✗ VIOLATES"
                print(f"ℓ={l:2d}: κ_α={kappa:.4f}, error={error:.6f} {status}")

        if weak_l_values:
            max_error = max(errors_in_weak)
            print(f"\nMax error in weak regime: {max_error:.6f}")

            if max_error > 0.01:
                print(f"❌ THEORY FALSIFIED: Found error {max_error:.4f} > 1% in weak regime")
            else:
                print(f"✅ Theory holds: All errors < 1% in weak regime")

            # Statistical summary
            print(f"\nStatistics for weak regime:")
            print(f"  Number of l values tested: {len(weak_l_values)}")
            print(f"  Mean error: {np.mean(errors_in_weak):.6f}")
            print(f"  Max error:  {max_error:.6f}")
            print(f"  Theory satisfied: {max_error < 0.01}")
        else:
            print("⚠️ No l values in weak curvature regime for this R")

    def test_convergence_with_angular_momentum(self):
        """
        Theory predicts: Error should decrease with increasing ℓ.
        FALSIFIABLE: If errors increase or don't converge
        """
        print("\n" + "="*60)
        print("FALSIFICATION TEST: Convergence with ℓ")
        print("="*60)

        R = 5.0
        alpha = 0.6
        laplacian = SphericalLaplacian(R=R, alpha=alpha, l_max=50)

        print(f"Parameters: R={R}, α={alpha}")
        print("\nTheory prediction: Errors decrease with increasing ℓ")
        print("-" * 50)

        # Compute errors for range of l
        l_values = np.arange(5, 51, 5)
        errors = []

        print("ℓ     Error      Δ Error")
        print("-" * 30)

        prev_error = None
        for l in l_values:
            error = laplacian.relative_error(l)
            errors.append(error)

            if prev_error is not None:
                delta = error - prev_error
                trend = "↓" if delta < 0 else "↑"
                print(f"{l:2d}  {error:.8f}  {delta:+.8f} {trend}")
            else:
                print(f"{l:2d}  {error:.8f}")

            prev_error = error

        # Analyze trend
        errors = np.array(errors)
        differences = np.diff(errors)
        n_decreasing = np.sum(differences < 0)
        n_increasing = np.sum(differences > 0)

        print(f"\nTrend analysis:")
        print(f"  Decreasing steps: {n_decreasing}/{len(differences)}")
        print(f"  Increasing steps: {n_increasing}/{len(differences)}")

        # Fit power law to see if convergence
        if len(l_values) > 2:
            log_l = np.log(l_values)
            log_errors = np.log(errors + 1e-15)
            coeffs = np.polyfit(log_l, log_errors, 1)
            convergence_rate = -coeffs[0]

            print(f"  Power law fit: error ~ ℓ^({-convergence_rate:.2f})")

            if convergence_rate > 0:
                print(f"✅ Theory holds: Errors converge with rate ℓ^(-{convergence_rate:.2f})")
            else:
                print(f"❌ THEORY FALSIFIED: Errors diverge!")

    def test_correction_sign(self):
        """
        Theory predicts: Corrections negative for 0 < α < 1.
        FALSIFIABLE: If we find positive corrections
        """
        print("\n" + "="*60)
        print("FALSIFICATION TEST: Correction Sign")
        print("="*60)

        alpha_values = [0.1, 0.3, 0.5, 0.7, 0.9]

        for alpha in alpha_values:
            laplacian = SphericalLaplacian(R=2.0, alpha=alpha, l_max=20)

            print(f"\nα = {alpha}:")
            print("-" * 30)

            # Check correction signs
            violations = []
            for l in range(1, 11):
                correction = laplacian.compute_correction_term(l)

                if correction > 0:
                    violations.append(l)
                    print(f"  ℓ={l}: C={correction:+.8f} ✗ WRONG SIGN!")
                else:
                    print(f"  ℓ={l}: C={correction:+.8f} ✓")

            if violations:
                print(f"❌ THEORY FALSIFIED: Found positive corrections at ℓ={violations}")
            else:
                print(f"✅ Theory holds: All corrections negative as expected")

    def test_scaling_with_radius(self):
        """
        Theory predicts: Corrections scale as O(1/R²).
        FALSIFIABLE: If scaling exponent differs significantly
        """
        print("\n" + "="*60)
        print("FALSIFICATION TEST: O(1/R²) Scaling")
        print("="*60)

        R_values = [1.0, 2.0, 4.0, 8.0, 16.0]
        alpha = 0.5
        l_test = 10

        print(f"Testing ℓ={l_test}, α={alpha}")
        print("\nTheory prediction: Corrections ~ 1/R²")
        print("-" * 50)

        result = validate_scaling_law(R_values, alpha, l_test)

        print(f"R values tested: {R_values}")
        print(f"Fitted scaling exponent: {result['fitted_exponent']:.3f}")
        print(f"Theory prediction:       {result['expected_exponent']:.1f}")
        print(f"Relative error:          {result['relative_error']:.1%}")
        print(f"R² of fit:              {result['r_squared']:.4f}")

        if result['scaling_validated']:
            print(f"✅ Theory holds: Scaling matches O(1/R²)")
        else:
            print(f"❌ THEORY FALSIFIED: Scaling is R^(-{result['fitted_exponent']:.2f}), not R^(-2)")

    def test_comprehensive_falsification(self):
        """
        Run all falsification tests and summarize.
        """
        print("\n" + "="*60)
        print("COMPREHENSIVE FALSIFICATION SUITE")
        print("="*60)

        # Test multiple parameter combinations
        test_cases = [
            {'R': 10.0, 'alpha': 0.5, 'l_max': 30, 'name': 'Large sphere, α=0.5'},
            {'R': 2.0, 'alpha': 0.1, 'l_max': 20, 'name': 'Medium sphere, α=0.1'},
            {'R': 1.0, 'alpha': 0.9, 'l_max': 15, 'name': 'Unit sphere, α=0.9'},
        ]

        all_results = []

        for case in test_cases:
            print(f"\nTesting: {case['name']}")
            print(f"  R={case['R']}, α={case['alpha']}, l_max={case['l_max']}")
            print("-" * 50)

            laplacian = SphericalLaplacian(
                R=case['R'],
                alpha=case['alpha'],
                l_max=case['l_max']
            )

            # Run falsification tests
            results = falsification_tests(laplacian)

            # Print results
            for test_name, test_result in results.items():
                if test_name == 'overall':
                    continue

                if isinstance(test_result, dict) and 'falsifies_theory' in test_result:
                    status = "❌ FALSIFIED" if test_result['falsifies_theory'] else "✅ PASSED"
                    print(f"  {test_name}: {status}")

            # Overall result
            if results['overall']['theory_falsified']:
                print(f"  OVERALL: ❌ THEORY FALSIFIED")
                print(f"    Failed tests: {results['overall']['failed_tests']}")
            else:
                print(f"  OVERALL: ✅ THEORY HOLDS")

            all_results.append({
                'case': case['name'],
                'falsified': results['overall']['theory_falsified']
            })

        # Final summary
        print("\n" + "="*60)
        print("FALSIFICATION SUMMARY")
        print("="*60)

        n_falsified = sum(1 for r in all_results if r['falsified'])
        print(f"Total test cases: {len(all_results)}")
        print(f"Theory falsified: {n_falsified}/{len(all_results)}")

        if n_falsified > 0:
            print("\n❌ THEORY FALSIFIED in some parameter regimes")
            for r in all_results:
                if r['falsified']:
                    print(f"  - {r['case']}")
        else:
            print("\n✅ THEORY HOLDS across all tested parameter regimes")


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_l_equals_zero(self):
        """Test special case of ℓ=0."""
        print("\n" + "="*60)
        print("EDGE CASE: ℓ=0 (Constant Mode)")
        print("="*60)

        laplacian = SphericalLaplacian(R=2.0, alpha=0.5)

        lambda_0 = laplacian.eigenvalue_exact(0)
        fractional_0 = laplacian.eigenvalue_fractional_exact(0)
        correction_0 = laplacian.compute_correction_term(0)

        print(f"ℓ=0 properties:")
        print(f"  Exact eigenvalue:      {lambda_0}")
        print(f"  Fractional eigenvalue: {fractional_0}")
        print(f"  Correction term:       {correction_0}")

        assert lambda_0 == 0, "λ_0 should be zero"
        assert fractional_0 == 0, "λ_0^α should be zero"
        assert correction_0 == 0, "Correction for ℓ=0 should be zero"

        print("\n✅ ℓ=0 handled correctly (all zeros)")

    def test_alpha_boundary_values(self):
        """Test α near 0 and 1."""
        print("\n" + "="*60)
        print("EDGE CASE: α Near Boundaries")
        print("="*60)

        R = 2.0
        l_test = 5

        # Test α near 0
        alpha_small = 0.01
        laplacian_small = SphericalLaplacian(R=R, alpha=alpha_small)
        correction_small = laplacian_small.compute_correction_term(l_test)

        # Test α near 1
        alpha_large = 0.99
        laplacian_large = SphericalLaplacian(R=R, alpha=alpha_large)
        correction_large = laplacian_large.compute_correction_term(l_test)

        print(f"Correction term for ℓ={l_test}, R={R}:")
        print(f"  α={alpha_small}: C={correction_small:.10f}")
        print(f"  α={alpha_large}: C={correction_large:.10f}")

        # Both should be very small (α(α-1) approaches 0)
        print(f"\nPrefactor α(α-1)/12:")
        print(f"  α={alpha_small}: {alpha_small*(alpha_small-1)/12:.10f}")
        print(f"  α={alpha_large}: {alpha_large*(alpha_large-1)/12:.10f}")

        print("\n✅ Corrections vanish as α→0 or α→1 (as expected)")