"""
📄 File: main.py
Purpose: Main script for E53 Sphere Validation experiment
Created: November 11, 2025
Used by: Command line execution

Usage:
    python main.py              # Run default analysis
    python main.py --R 5.0      # Custom radius
    python main.py --alpha 0.6  # Custom fractional order
    python main.py --plots      # Generate all plots
"""

import numpy as np
import json
import argparse
from pathlib import Path
from datetime import datetime


class NumpyEncoder(json.JSONEncoder):
    """Custom JSON encoder for numpy types."""
    def default(self, obj):
        if isinstance(obj, (np.bool_, np.bool)):
            return bool(obj)
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)

from src.core import SphericalLaplacian
from src.compute import (
    compute_convergence_analysis,
    compute_parameter_sweep,
    compute_spectrum_comparison
)
from src.validate import (
    create_validation_report,
    falsification_tests,
    validate_scaling_law
)
from src.visualize import (
    plot_error_vs_l,
    plot_correction_magnitude,
    plot_convergence_analysis,
    plot_parameter_sweep,
    create_summary_figure
)


def ensure_output_dirs():
    """Create output directories if they don't exist."""
    dirs = ['outputs', 'outputs/plots', 'outputs/data']
    for dir_path in dirs:
        Path(dir_path).mkdir(parents=True, exist_ok=True)


def run_main_analysis(R: float = 5.0, alpha: float = 0.5, l_max: int = 50):
    """
    Run main analysis for given parameters.

    Args:
        R: Sphere radius
        alpha: Fractional order
        l_max: Maximum angular momentum
    """
    print("=" * 70)
    print("E53: SPHERE VALIDATION (S²)")
    print("Pure Analytical Validation of Fractional Laplacian")
    print("=" * 70)

    # Create Laplacian
    print(f"\nInitializing SphericalLaplacian:")
    print(f"  Radius R = {R}")
    print(f"  Fractional order α = {alpha}")
    print(f"  Maximum ℓ = {l_max}")

    laplacian = SphericalLaplacian(R=R, alpha=alpha, l_max=l_max)

    # === DISCOVERY PHASE ===
    print("\n" + "="*70)
    print("DISCOVERY PHASE: Learning System Properties")
    print("="*70)

    # Discover eigenvalue structure
    print("\n1. Eigenvalue Structure:")
    print("-" * 40)
    for l in [0, 1, 2, 5, 10, 20]:
        if l <= l_max:
            exact = laplacian.eigenvalue_exact(l)
            fractional = laplacian.eigenvalue_fractional_exact(l)
            corrected = laplacian.eigenvalue_fractional_corrected(l)
            correction = laplacian.compute_correction_term(l)

            print(f"ℓ={l:2d}:")
            print(f"  λ_ℓ = {exact:.6f} (exact)")
            print(f"  λ_ℓ^α = {fractional:.6f} (fractional)")
            print(f"  Correction = {correction:.8f}")
            print(f"  Corrected = {corrected:.6f}")

    # Discover regime boundaries
    print("\n2. Curvature Regime Discovery:")
    print("-" * 40)
    l_weak_boundary = int(np.sqrt(0.1 * R**2))
    l_strong_boundary = int(np.sqrt(1.0 * R**2))

    print(f"For R={R}:")
    print(f"  Weak curvature (κ<0.1): ℓ < {l_weak_boundary}")
    print(f"  Intermediate: {l_weak_boundary} ≤ ℓ ≤ {l_strong_boundary}")
    print(f"  Strong curvature (κ>1): ℓ > {l_strong_boundary}")

    # === VALIDATION PHASE ===
    print("\n" + "="*70)
    print("VALIDATION PHASE: Testing Theory Predictions")
    print("="*70)

    # Generate validation report
    report = create_validation_report(laplacian)

    print("\n1. Error Analysis:")
    print("-" * 40)
    print(f"  Max relative error: {report['summary']['max_error_overall']:.6f}")
    print(f"  Mean relative error: {report['summary']['mean_error_overall']:.6f}")

    if report['summary']['max_error_weak_regime'] is not None:
        print(f"  Max error in weak regime: {report['summary']['max_error_weak_regime']:.6f}")

    print("\n2. Theory Tests:")
    print("-" * 40)
    theory_tests = report['theory_tests']
    for test_name, result in theory_tests.items():
        if test_name not in ['all_tests_passed', 'tests_performed', 'tests_passed', 'note']:
            if result is True:
                status = "✅ PASSED"
            elif result is False:
                status = "❌ FAILED"
            else:
                status = "⚪ N/A"
            print(f"  {test_name}: {status}")

    # === FALSIFICATION PHASE ===
    print("\n" + "="*70)
    print("FALSIFICATION PHASE: Attempting to Disprove Theory")
    print("="*70)

    falsification = falsification_tests(laplacian)

    for test_name, test_result in falsification.items():
        if test_name == 'overall':
            continue
        if isinstance(test_result, dict) and 'falsifies_theory' in test_result:
            if test_result['falsifies_theory']:
                print(f"  {test_name}: ❌ FALSIFIES THEORY")
                print(f"    {test_result.get('description', '')}")
            else:
                print(f"  {test_name}: ✅ Theory holds")

    # Overall result
    if falsification['overall']['theory_falsified']:
        print(f"\n❌ THEORY FALSIFIED")
        print(f"  Failed tests: {falsification['overall']['failed_tests']}")
    else:
        print(f"\n✅ THEORY VALIDATED - All falsification tests passed")

    # === SCALING ANALYSIS ===
    print("\n" + "="*70)
    print("SCALING ANALYSIS: O(1/R²) Convergence")
    print("="*70)

    R_values = [R/2, R, R*2, R*4]
    l_test = min(10, l_max//2)

    scaling_result = validate_scaling_law(R_values, alpha, l_test)

    print(f"Testing with R values: {R_values}")
    print(f"Fitted scaling: R^(-{scaling_result['fitted_exponent']:.3f})")
    print(f"Theory prediction: R^(-2)")
    print(f"Relative error: {scaling_result['relative_error']:.1%}")

    if scaling_result['scaling_validated']:
        print("✅ Scaling matches theoretical prediction")
    else:
        print("❌ Scaling deviates from theory")

    # Save results
    save_results(laplacian, report, falsification, scaling_result)

    return laplacian, report


def save_results(laplacian, report, falsification, scaling_result):
    """Save all results to files."""
    ensure_output_dirs()

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    # Save validation report
    report_file = f"outputs/data/validation_report_{timestamp}.json"
    with open(report_file, 'w') as f:
        json.dump(report, f, indent=2, cls=NumpyEncoder)
    print(f"\n✓ Saved validation report to {report_file}")

    # Save falsification results
    falsification_file = f"outputs/data/falsification_{timestamp}.json"
    with open(falsification_file, 'w') as f:
        json.dump(falsification, f, indent=2, cls=NumpyEncoder)
    print(f"✓ Saved falsification results to {falsification_file}")

    # Save scaling analysis
    scaling_file = f"outputs/data/scaling_analysis_{timestamp}.json"
    with open(scaling_file, 'w') as f:
        json.dump(scaling_result, f, indent=2, cls=NumpyEncoder)
    print(f"✓ Saved scaling analysis to {scaling_file}")

    # Save spectrum summary
    summary = laplacian.get_spectrum_summary()
    summary_file = f"outputs/data/spectrum_summary_{timestamp}.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2, cls=NumpyEncoder)
    print(f"✓ Saved spectrum summary to {summary_file}")


def generate_plots(laplacian):
    """Generate all visualization plots."""
    ensure_output_dirs()

    print("\nGenerating plots...")

    # 1. Error vs l plot
    fig1 = plot_error_vs_l(laplacian, save_path="outputs/plots/error_vs_l.png")
    print("✓ Generated error vs ℓ plot")

    # 2. Correction magnitude plot
    fig2 = plot_correction_magnitude(laplacian, save_path="outputs/plots/correction_magnitude.png")
    print("✓ Generated correction magnitude plot")

    # 3. Convergence analysis
    R_values = [1.0, 2.0, 4.0, 8.0]
    fig3 = plot_convergence_analysis(
        R_values, laplacian.alpha, l_fixed=10,
        save_path="outputs/plots/convergence_analysis.png"
    )
    print("✓ Generated convergence analysis plot")

    # 4. Parameter sweep
    R_sweep = [1.0, 2.0, 3.0, 5.0, 10.0]
    alpha_sweep = [0.1, 0.3, 0.5, 0.7, 0.9]
    fig4 = plot_parameter_sweep(
        R_sweep, alpha_sweep, l_test=10,
        save_path="outputs/plots/parameter_sweep.png"
    )
    print("✓ Generated parameter sweep heatmap")

    # 5. Summary figure
    fig5 = create_summary_figure(laplacian, save_path="outputs/plots/summary.png")
    print("✓ Generated summary figure")

    print("\nAll plots saved to outputs/plots/")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="E53 Sphere Validation - Analytical validation of fractional Laplacian"
    )
    parser.add_argument('--R', type=float, default=5.0,
                       help='Sphere radius (default: 5.0)')
    parser.add_argument('--alpha', type=float, default=0.5,
                       help='Fractional order (default: 0.5)')
    parser.add_argument('--l_max', type=int, default=50,
                       help='Maximum angular momentum (default: 50)')
    parser.add_argument('--plots', action='store_true',
                       help='Generate visualization plots')
    parser.add_argument('--sweep', action='store_true',
                       help='Run parameter sweep analysis')

    args = parser.parse_args()

    # Validate inputs
    if args.R <= 0:
        print("Error: R must be positive")
        return 1

    if not 0 < args.alpha < 1:
        print("Error: alpha must be in (0,1)")
        return 1

    if args.l_max < 0:
        print("Error: l_max must be non-negative")
        return 1

    # Run main analysis
    laplacian, report = run_main_analysis(args.R, args.alpha, args.l_max)

    # Generate plots if requested
    if args.plots:
        generate_plots(laplacian)

    # Run parameter sweep if requested
    if args.sweep:
        print("\n" + "="*70)
        print("PARAMETER SWEEP ANALYSIS")
        print("="*70)

        R_values = [1.0, 2.0, 5.0, 10.0]
        alpha_values = [0.1, 0.3, 0.5, 0.7, 0.9]

        sweep_results = compute_parameter_sweep(R_values, alpha_values, l_test=10)

        # Find optimal parameters
        best = min(sweep_results['sweep_data'], key=lambda x: x['relative_error'])
        print(f"\nOptimal parameters (lowest error):")
        print(f"  R = {best['R']}")
        print(f"  α = {best['alpha']}")
        print(f"  Relative error = {best['relative_error']:.6f}")
        print(f"  Regime = {best['regime']}")

        # Save sweep results
        ensure_output_dirs()
        sweep_file = "outputs/data/parameter_sweep.json"
        with open(sweep_file, 'w') as f:
            json.dump(sweep_results, f, indent=2, cls=NumpyEncoder)
        print(f"\n✓ Saved parameter sweep to {sweep_file}")

    print("\n" + "="*70)
    print("E53 EXPERIMENT COMPLETE")
    print("="*70)

    # Final summary
    if report['summary']['theory_validated']:
        print("✅ Theory validated: Corrections improve fractional eigenvalue accuracy")
    else:
        print("❌ Theory not fully validated in tested regime")

    return 0


if __name__ == "__main__":
    exit(main())