"""
E53 Sphere Validation S² - Pure analytical validation of fractional Laplacian.
"""

from .core import SphericalLaplacian
from .compute import (
    compute_exact_spectrum_sphere,
    compute_fractional_spectrum_exact,
    compute_fractional_spectrum_corrected,
    compute_curvature_correction,
    compute_weak_curvature_parameter,
    compute_convergence_analysis,
    compute_spectrum_comparison,
    compute_parameter_sweep
)
from .validate import (
    validate_weak_curvature_regime,
    validate_error_bounds,
    validate_convergence_with_l,
    create_validation_report,
    validate_scaling_law,
    falsification_tests
)
from .visualize import (
    plot_error_vs_l,
    plot_correction_magnitude,
    plot_convergence_analysis,
    plot_parameter_sweep,
    plot_eigenvalue_comparison,
    create_summary_figure
)

__all__ = [
    # Core
    'SphericalLaplacian',

    # Compute
    'compute_exact_spectrum_sphere',
    'compute_fractional_spectrum_exact',
    'compute_fractional_spectrum_corrected',
    'compute_curvature_correction',
    'compute_weak_curvature_parameter',
    'compute_convergence_analysis',
    'compute_spectrum_comparison',
    'compute_parameter_sweep',

    # Validate
    'validate_weak_curvature_regime',
    'validate_error_bounds',
    'validate_convergence_with_l',
    'create_validation_report',
    'validate_scaling_law',
    'falsification_tests',

    # Visualize
    'plot_error_vs_l',
    'plot_correction_magnitude',
    'plot_convergence_analysis',
    'plot_parameter_sweep',
    'plot_eigenvalue_comparison',
    'create_summary_figure'
]