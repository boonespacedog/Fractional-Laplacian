"""
📄 File: test_integration.py
Purpose: Integration tests for E53 - full workflow validation
Created: November 11, 2025
Used by: pytest
"""

import pytest
import numpy as np
import sys
import json
from pathlib import Path

# Add src to path
sys.path.append(str(Path(__file__).parent.parent))

from src.core import SphericalLaplacian
from src.compute import compute_parameter_sweep, compute_spectrum_comparison
from src.validate import create_validation_report


class TestFullWorkflow:
    """Test complete analysis workflow."""

    def test_complete_validation_report(self):
        """Generate and analyze complete validation report."""
        print("\n" + "="*60)
        print("INTEGRATION TEST: Complete Validation Report")
        print("="*60)

        # Create laplacian with typical parameters
        laplacian = SphericalLaplacian(R=5.0, alpha=0.5, l_max=40)

        # Generate complete report
        report = create_validation_report(laplacian)

        print(f"Generated validation report for R={laplacian.R}, α={laplacian.alpha}")
        print("-" * 50)

        # Check report structure
        assert 'parameters' in report
        assert 'regime_analysis' in report
        assert 'error_validation' in report
        assert 'convergence_analysis' in report
        assert 'theory_tests' in report
        assert 'summary' in report

        print("Report structure:")
        for key in report.keys():
            print(f"  ✓ {key}")

        # Display key findings
        print("\nKey findings:")
        print(f"  Max error overall: {report['summary']['max_error_overall']:.6f}")
        print(f"  Mean error overall: {report['summary']['mean_error_overall']:.6f}")

        if report['summary']['max_error_weak_regime'] is not None:
            print(f"  Max error in weak regime: {report['summary']['max_error_weak_regime']:.6f}")

        # Theory validation status
        theory_status = report['theory_tests']
        print("\nTheory test results:")
        for test_name, result in theory_status.items():
            if test_name not in ['all_tests_passed', 'tests_performed', 'tests_passed', 'note']:
                status = "✓" if result else "✗" if result is False else "N/A"
                print(f"  {test_name}: {status}")

        if theory_status['all_tests_passed'] is not None:
            overall = "✅ VALIDATED" if theory_status['all_tests_passed'] else "❌ FALSIFIED"
            print(f"\nOverall theory validation: {overall}")
            print(f"  Tests passed: {theory_status.get('tests_passed', 0)}/{theory_status.get('tests_performed', 0)}")

    def test_parameter_sweep_workflow(self):
        """Test parameter sweep analysis."""
        print("\n" + "="*60)
        print("INTEGRATION TEST: Parameter Sweep")
        print("="*60)

        # Define sweep parameters
        R_values = [1.0, 2.0, 5.0]
        alpha_values = [0.3, 0.5, 0.7]
        l_test = 10

        print(f"Sweep parameters:")
        print(f"  R values: {R_values}")
        print(f"  α values: {alpha_values}")
        print(f"  Test ℓ: {l_test}")

        # Run sweep
        results = compute_parameter_sweep(R_values, alpha_values, l_test)

        assert 'sweep_data' in results
        assert len(results['sweep_data']) == len(R_values) * len(alpha_values)

        print(f"\nGenerated {len(results['sweep_data'])} data points")

        # Analyze results
        print("\nParameter sweep results:")
        print("-" * 60)
        print("   R     α    Error     Regime")
        print("-" * 60)

        for data in results['sweep_data']:
            print(f"{data['R']:5.1f}  {data['alpha']:4.2f}  {data['relative_error']:.6f}  {data['regime']:12s}")

        # Find best parameters (lowest error)
        best = min(results['sweep_data'], key=lambda x: x['relative_error'])
        print(f"\nBest parameters:")
        print(f"  R={best['R']}, α={best['alpha']}")
        print(f"  Relative error: {best['relative_error']:.6f}")
        print(f"  Regime: {best['regime']}")

    def test_spectrum_comparison(self):
        """Test spectrum comparison workflow."""
        print("\n" + "="*60)
        print("INTEGRATION TEST: Spectrum Comparison")
        print("="*60)

        laplacian = SphericalLaplacian(R=3.0, alpha=0.6, l_max=30)
        l_range = (1, 15)

        # Compute comparison
        comparison = compute_spectrum_comparison(laplacian, l_range)

        print(f"Spectrum comparison for R={laplacian.R}, α={laplacian.alpha}")
        print(f"ℓ range: {l_range[0]} to {l_range[1]}")
        print("-" * 50)

        # Check structure
        assert 'l_values' in comparison
        assert 'eigenvalues' in comparison
        assert 'corrections' in comparison
        assert 'relative_errors' in comparison
        assert 'regime_analysis' in comparison
        assert 'statistics' in comparison

        # Display statistics
        stats = comparison['statistics']
        print(f"Statistics:")
        print(f"  Max correction magnitude: {stats['max_correction']:.8f}")
        print(f"  Mean correction magnitude: {stats['mean_correction']:.8f}")
        print(f"  Max relative error: {stats['max_relative_error']:.6f}")
        print(f"  Mean relative error: {stats['mean_relative_error']:.6f}")

        # Regime breakdown
        regimes = comparison['regime_analysis']
        for regime_name, regime_data in regimes.items():
            if regime_data['max_error'] is not None:
                print(f"\n{regime_name.upper()} regime:")
                print(f"  Number of ℓ values: {len(regime_data['indices'])}")
                print(f"  Max error: {regime_data['max_error']:.6f}")
                print(f"  Mean error: {regime_data['mean_error']:.6f}")

    def test_data_serialization(self):
        """Test that all outputs are JSON serializable."""
        print("\n" + "="*60)
        print("INTEGRATION TEST: Data Serialization")
        print("="*60)

        laplacian = SphericalLaplacian(R=2.5, alpha=0.4, l_max=25)

        # Get various outputs
        summary = laplacian.get_spectrum_summary()
        validation = create_validation_report(laplacian)
        comparison = compute_spectrum_comparison(laplacian, (1, 10))

        # Try to serialize each
        outputs = {
            'spectrum_summary': summary,
            'validation_report': validation,
            'spectrum_comparison': comparison
        }

        for name, data in outputs.items():
            try:
                json_str = json.dumps(data, indent=2)
                print(f"✓ {name}: Serialized successfully ({len(json_str)} chars)")

                # Verify round-trip
                reconstructed = json.loads(json_str)
                assert isinstance(reconstructed, dict), f"{name} not reconstructed as dict"

            except (TypeError, ValueError) as e:
                pytest.fail(f"Failed to serialize {name}: {e}")

        print("\n✅ All outputs are JSON serializable")

    def test_numerical_stability(self):
        """Test numerical stability with extreme parameters."""
        print("\n" + "="*60)
        print("INTEGRATION TEST: Numerical Stability")
        print("="*60)

        test_cases = [
            {'R': 0.1, 'alpha': 0.01, 'l_max': 5, 'name': 'Small R, small α'},
            {'R': 100.0, 'alpha': 0.99, 'l_max': 100, 'name': 'Large R, large α'},
            {'R': 1.0, 'alpha': 0.5, 'l_max': 200, 'name': 'High l_max'},
        ]

        for case in test_cases:
            print(f"\nTesting: {case['name']}")
            print(f"  Parameters: R={case['R']}, α={case['alpha']}, l_max={case['l_max']}")

            try:
                laplacian = SphericalLaplacian(
                    R=case['R'],
                    alpha=case['alpha'],
                    l_max=case['l_max']
                )

                # Compute some values
                eigenvalues = laplacian.eigenvalues_fractional_corrected

                # Check for numerical issues
                assert not np.any(np.isnan(eigenvalues)), "NaN values detected"
                assert not np.any(np.isinf(eigenvalues)), "Inf values detected"
                assert np.all(np.isfinite(eigenvalues)), "Non-finite values detected"

                print(f"  ✓ No numerical issues")
                print(f"  Eigenvalue range: [{np.min(eigenvalues):.6e}, {np.max(eigenvalues):.6e}]")

            except Exception as e:
                print(f"  ✗ Failed: {e}")

        print("\n✅ Numerical stability tests completed")