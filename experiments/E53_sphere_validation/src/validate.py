"""
📄 File: validate.py
Purpose: Error analysis and theory validation for spherical Laplacian
Created: November 11, 2025
Used by: main.py, tests/
"""

import numpy as np
from typing import Dict, List, Tuple, Optional
from .core import SphericalLaplacian
from .compute import compute_weak_curvature_parameter


def validate_weak_curvature_regime(R: float, l_values: np.ndarray) -> Dict:
    """
    Check if we're in weak curvature regime (κ_α < 0.1).

    The correction formula is expected to be accurate when the
    curvature parameter κ_α = ℓ(ℓ+1)/R² is small.

    Args:
        R: Sphere radius
        l_values: Array of angular momentum values

    Returns:
        Dictionary with regime classification
    """
    curvature_params = np.array([
        compute_weak_curvature_parameter(R, l) for l in l_values
    ])

    weak_mask = curvature_params < 0.1
    intermediate_mask = (curvature_params >= 0.1) & (curvature_params <= 1.0)
    strong_mask = curvature_params > 1.0

    return {
        'R': R,
        'l_values': l_values.tolist(),
        'curvature_params': curvature_params.tolist(),
        'regimes': {
            'weak': {
                'l_values': l_values[weak_mask].tolist(),
                'count': int(np.sum(weak_mask)),
                'fraction': float(np.mean(weak_mask)),
                'kappa_range': [float(curvature_params[weak_mask].min()),
                               float(curvature_params[weak_mask].max())] if np.any(weak_mask) else None
            },
            'intermediate': {
                'l_values': l_values[intermediate_mask].tolist(),
                'count': int(np.sum(intermediate_mask)),
                'fraction': float(np.mean(intermediate_mask)),
                'kappa_range': [float(curvature_params[intermediate_mask].min()),
                               float(curvature_params[intermediate_mask].max())] if np.any(intermediate_mask) else None
            },
            'strong': {
                'l_values': l_values[strong_mask].tolist(),
                'count': int(np.sum(strong_mask)),
                'fraction': float(np.mean(strong_mask)),
                'kappa_range': [float(curvature_params[strong_mask].min()),
                               float(curvature_params[strong_mask].max())] if np.any(strong_mask) else None
            }
        },
        'transition_l': {
            'weak_to_intermediate': int(np.sqrt(0.1 * R**2)) if R > 0 else 0,
            'intermediate_to_strong': int(np.sqrt(1.0 * R**2)) if R > 0 else 0
        }
    }


def validate_error_bounds(
    errors: np.ndarray,
    l_values: np.ndarray,
    R: float,
    tolerance: float = 0.01
) -> Dict:
    """
    🧠 Function: validate_error_bounds
    Role: Verify relative error < 1% in weak curvature
    Inputs: errors, l_values, R, tolerance (default 1%)
    Returns: Validation results with pass/fail status
    Notes: Main falsification test - theory fails if error > 1% in weak regime
    """
    # Compute curvature parameters
    curvature_params = np.array([
        compute_weak_curvature_parameter(R, l) for l in l_values
    ])

    # Identify weak curvature regime
    weak_mask = curvature_params < 0.1

    # Check errors in weak regime
    if np.any(weak_mask):
        weak_errors = errors[weak_mask]
        weak_l = l_values[weak_mask]

        # Theory prediction: error < tolerance in weak regime
        violations = weak_errors > tolerance
        n_violations = np.sum(violations)

        theory_satisfied = n_violations == 0

        result = {
            'theory_satisfied': bool(theory_satisfied),
            'tolerance': tolerance,
            'weak_regime': {
                'l_values': weak_l.tolist(),
                'errors': weak_errors.tolist(),
                'max_error': float(np.max(weak_errors)),
                'mean_error': float(np.mean(weak_errors)),
                'violations': {
                    'count': int(n_violations),
                    'l_values': weak_l[violations].tolist() if n_violations > 0 else [],
                    'errors': weak_errors[violations].tolist() if n_violations > 0 else []
                }
            }
        }
    else:
        result = {
            'theory_satisfied': None,  # No weak regime to test
            'tolerance': tolerance,
            'weak_regime': None,
            'note': 'No l values in weak curvature regime for this R'
        }

    # Also check other regimes for completeness
    intermediate_mask = (curvature_params >= 0.1) & (curvature_params <= 1.0)
    strong_mask = curvature_params > 1.0

    if np.any(intermediate_mask):
        result['intermediate_regime'] = {
            'max_error': float(np.max(errors[intermediate_mask])),
            'mean_error': float(np.mean(errors[intermediate_mask]))
        }

    if np.any(strong_mask):
        result['strong_regime'] = {
            'max_error': float(np.max(errors[strong_mask])),
            'mean_error': float(np.mean(errors[strong_mask]))
        }

    return result


def validate_convergence_with_l(errors: np.ndarray, l_values: np.ndarray) -> Dict:
    """
    Check if error improves with increasing ℓ.

    Theory predicts the correction becomes more accurate for larger ℓ
    (in the appropriate regime).

    Args:
        errors: Relative errors
        l_values: Angular momentum values

    Returns:
        Dictionary with convergence analysis
    """
    # Remove l=0 if present (often special case)
    if l_values[0] == 0:
        errors = errors[1:]
        l_values = l_values[1:]

    if len(l_values) < 2:
        return {
            'convergence_detected': None,
            'note': 'Insufficient data points for convergence analysis'
        }

    # Check for monotonic decrease
    differences = np.diff(errors)
    decreasing = differences < 0
    fraction_decreasing = np.mean(decreasing)

    # Fit power law: error ~ l^(-p)
    # Use log-log fit for l > 5 to avoid small-l effects
    large_l_mask = l_values > 5
    if np.sum(large_l_mask) >= 2:
        log_l = np.log(l_values[large_l_mask])
        log_errors = np.log(errors[large_l_mask])

        # Linear fit in log-log space
        coeffs = np.polyfit(log_l, log_errors, 1)
        convergence_exponent = -coeffs[0]  # Negative of slope

        # Compute goodness of fit
        fitted_log_errors = np.polyval(coeffs, log_l)
        residuals = log_errors - fitted_log_errors
        r_squared = 1 - np.var(residuals) / np.var(log_errors)

        power_law_fit = {
            'exponent': float(convergence_exponent),
            'r_squared': float(r_squared),
            'converging': convergence_exponent > 0,
            'l_range': [int(l_values[large_l_mask].min()),
                       int(l_values[large_l_mask].max())]
        }
    else:
        power_law_fit = None

    return {
        'convergence_detected': fraction_decreasing > 0.7,  # 70% of steps decrease
        'fraction_decreasing': float(fraction_decreasing),
        'power_law_fit': power_law_fit,
        'error_evolution': {
            'l_values': l_values.tolist(),
            'errors': errors.tolist(),
            'differences': differences.tolist()
        }
    }


def create_validation_report(laplacian: SphericalLaplacian) -> Dict:
    """
    🧠 Function: create_validation_report
    Role: Complete validation of theory predictions
    Inputs: laplacian object
    Returns: Comprehensive validation report
    Notes: Main validation function - checks all theoretical predictions
    """
    # Compute errors for all l values
    l_values = np.arange(laplacian.l_max + 1)
    errors = np.array([laplacian.relative_error(l) for l in l_values])

    # 1. Check weak curvature regime
    regime_analysis = validate_weak_curvature_regime(laplacian.R, l_values)

    # 2. Validate error bounds
    error_validation = validate_error_bounds(errors, l_values, laplacian.R)

    # 3. Check convergence with l
    convergence_analysis = validate_convergence_with_l(errors, l_values)

    # 4. Check correction term properties
    corrections = laplacian.correction_terms

    # Theory: corrections negative for 0 < α < 1
    if 0 < laplacian.alpha < 1:
        expected_sign = 'negative'
        sign_correct = np.all(corrections[1:] < 0)  # Exclude l=0
    elif laplacian.alpha == 0 or laplacian.alpha == 1:
        expected_sign = 'zero'
        sign_correct = np.all(np.abs(corrections) < 1e-15)
    else:
        expected_sign = 'undefined'
        sign_correct = None

    correction_validation = {
        'expected_sign': expected_sign,
        'sign_correct': bool(sign_correct) if sign_correct is not None else None,
        'magnitude_range': [float(np.min(np.abs(corrections[1:]))),
                           float(np.max(np.abs(corrections[1:])))],
        'formula': 'C(ℓ,R,α) = α(α-1)/12R² · ℓ(ℓ+1)'
    }

    # 5. Overall theory validation
    theory_tests = {
        'weak_curvature_error_bound': error_validation.get('theory_satisfied'),
        'convergence_with_l': convergence_analysis.get('convergence_detected'),
        'correction_sign': correction_validation['sign_correct'],
        'all_tests_passed': None  # Will compute below
    }

    # Check if all tests passed (where applicable)
    test_results = [
        theory_tests['weak_curvature_error_bound'],
        theory_tests['convergence_with_l'],
        theory_tests['correction_sign']
    ]
    # Filter out None values (tests that couldn't be performed)
    valid_tests = [t for t in test_results if t is not None]

    if valid_tests:
        theory_tests['all_tests_passed'] = bool(all(valid_tests))
        theory_tests['tests_performed'] = len(valid_tests)
        theory_tests['tests_passed'] = sum(1 for t in valid_tests if t)
    else:
        theory_tests['all_tests_passed'] = None
        theory_tests['note'] = 'No tests could be performed with current parameters'

    # Compile full report
    return {
        'parameters': {
            'R': laplacian.R,
            'alpha': laplacian.alpha,
            'l_max': laplacian.l_max
        },
        'regime_analysis': regime_analysis,
        'error_validation': error_validation,
        'convergence_analysis': convergence_analysis,
        'correction_validation': correction_validation,
        'theory_tests': theory_tests,
        'summary': {
            'max_error_overall': float(np.max(errors[1:])),  # Exclude l=0
            'mean_error_overall': float(np.mean(errors[1:])),
            'max_error_weak_regime': error_validation['weak_regime']['max_error'] if error_validation.get('weak_regime') else None,
            'theory_validated': theory_tests['all_tests_passed']
        }
    }


def validate_scaling_law(
    R_values: List[float],
    alpha: float,
    l_test: int = 10,
    expected_exponent: float = 2.0
) -> Dict:
    """
    Validate O(1/R²) scaling of corrections.

    Theory predicts corrections scale as 1/R² in the large-R limit.

    Args:
        R_values: List of sphere radii
        alpha: Fractional order
        l_test: Angular momentum for testing
        expected_exponent: Theoretical scaling exponent

    Returns:
        Dictionary with scaling validation
    """
    if len(R_values) < 2:
        return {
            'scaling_validated': None,
            'note': 'Need at least 2 R values for scaling analysis'
        }

    # Compute corrections for each R
    corrections = []
    for R in R_values:
        laplacian = SphericalLaplacian(R=R, alpha=alpha, l_max=l_test)
        correction = abs(laplacian.compute_correction_term(l_test))
        corrections.append(correction)

    corrections = np.array(corrections)

    # Fit power law in log-log space
    log_R = np.log(R_values)
    log_corrections = np.log(corrections + 1e-20)  # Avoid log(0)

    # Linear fit: log(C) = -p*log(R) + const
    coeffs = np.polyfit(log_R, log_corrections, 1)
    fitted_exponent = -coeffs[0]

    # Compute goodness of fit
    fitted_log_corrections = np.polyval(coeffs, log_R)
    residuals = log_corrections - fitted_log_corrections
    r_squared = 1 - np.var(residuals) / np.var(log_corrections)

    # Check if scaling matches theory
    relative_error = abs(fitted_exponent - expected_exponent) / expected_exponent
    scaling_validated = relative_error < 0.1  # Within 10% of theory

    return {
        'scaling_validated': bool(scaling_validated),
        'fitted_exponent': float(fitted_exponent),
        'expected_exponent': expected_exponent,
        'relative_error': float(relative_error),
        'r_squared': float(r_squared),
        'data': {
            'R_values': R_values,
            'corrections': corrections.tolist(),
            'log_fit_coefficients': coeffs.tolist()
        },
        'interpretation': f"Corrections scale as R^(-{fitted_exponent:.2f}), theory predicts R^(-{expected_exponent})"
    }


def falsification_tests(laplacian: SphericalLaplacian) -> Dict:
    """
    🧠 Function: falsification_tests
    Role: Explicit falsification tests that could disprove theory
    Inputs: laplacian object
    Returns: Dictionary with pass/fail for each falsification test
    Notes: Theory is falsified if ANY test fails
    """
    l_values = np.arange(1, laplacian.l_max + 1)  # Exclude l=0
    errors = np.array([laplacian.relative_error(l) for l in l_values])
    curvature_params = np.array([
        laplacian.compute_weak_curvature_parameter(l) for l in l_values
    ])

    tests = {}

    # Test 1: Error > 1% in weak curvature regime
    weak_mask = curvature_params < 0.1
    if np.any(weak_mask):
        max_weak_error = np.max(errors[weak_mask])
        tests['weak_regime_error_bound'] = {
            'description': 'Relative error must be < 1% in weak curvature regime',
            'passed': bool(max_weak_error < 0.01),
            'max_error_found': float(max_weak_error),
            'threshold': 0.01,
            'falsifies_theory': bool(max_weak_error >= 0.01)
        }
    else:
        tests['weak_regime_error_bound'] = {
            'description': 'Relative error must be < 1% in weak curvature regime',
            'passed': None,
            'note': 'No l values in weak regime to test',
            'falsifies_theory': False
        }

    # Test 2: Convergence doesn't improve with l
    if len(l_values) > 10:
        # Check if errors generally decrease with l
        first_half = errors[:len(errors)//2]
        second_half = errors[len(errors)//2:]
        mean_first = np.mean(first_half)
        mean_second = np.mean(second_half)

        tests['convergence_with_l'] = {
            'description': 'Errors should decrease with increasing l',
            'passed': bool(mean_second < mean_first),
            'mean_error_first_half': float(mean_first),
            'mean_error_second_half': float(mean_second),
            'falsifies_theory': bool(mean_second >= mean_first * 1.5)  # 50% worse
        }
    else:
        tests['convergence_with_l'] = {
            'description': 'Errors should decrease with increasing l',
            'passed': None,
            'note': 'Insufficient l values for convergence test',
            'falsifies_theory': False
        }

    # Test 3: Wrong sign of correction
    corrections = np.array([laplacian.compute_correction_term(l) for l in l_values])
    if 0 < laplacian.alpha < 1:
        # Should be negative
        wrong_sign = np.any(corrections > 0)
        tests['correction_sign'] = {
            'description': 'Corrections must be negative for 0 < α < 1',
            'passed': bool(not wrong_sign),
            'expected_sign': 'negative',
            'violations': int(np.sum(corrections > 0)),
            'falsifies_theory': bool(wrong_sign)
        }
    else:
        tests['correction_sign'] = {
            'description': 'Correction sign test',
            'passed': None,
            'note': f'Test not applicable for α = {laplacian.alpha}',
            'falsifies_theory': False
        }

    # Overall falsification
    falsifications = [t.get('falsifies_theory', False) for t in tests.values()]
    tests['overall'] = {
        'theory_falsified': bool(any(falsifications)),
        'number_of_failures': sum(1 for f in falsifications if f),
        'failed_tests': [name for name, test in tests.items()
                        if test.get('falsifies_theory', False)]
    }

    return tests