"""
📄 File: compute.py
Purpose: Eigenvalue and correction computations for spherical Laplacian
Created: November 11, 2025
Used by: validate.py, main.py
"""

import numpy as np
from typing import Tuple, Dict, List, Optional
from .core import SphericalLaplacian


def compute_exact_spectrum_sphere(R: float, l_max: int = 50) -> np.ndarray:
    """
    Compute exact eigenvalues λ_ℓ = ℓ(ℓ+1)/R² for ℓ=0,1,...,l_max.

    These are the textbook eigenvalues of the Laplace-Beltrami operator
    on a sphere of radius R. Each eigenvalue corresponds to spherical
    harmonics Y_ℓm with degeneracy 2ℓ+1.

    Args:
        R: Sphere radius
        l_max: Maximum angular momentum

    Returns:
        Array of eigenvalues λ_ℓ

    Theory: Standard spherical harmonic eigenvalues
    """
    if R <= 0:
        raise ValueError(f"Radius must be positive, got {R}")
    if l_max < 0:
        raise ValueError(f"l_max must be non-negative, got {l_max}")

    l_values = np.arange(l_max + 1)
    eigenvalues = l_values * (l_values + 1) / (R ** 2)

    return eigenvalues


def compute_fractional_spectrum_exact(eigenvalues: np.ndarray, alpha: float) -> np.ndarray:
    """
    🧠 Function: compute_fractional_spectrum_exact
    Role: Raise eigenvalues to fractional power
    Inputs: eigenvalues array, fractional order α
    Returns: λ^α for each eigenvalue
    Notes: Simple power law, no geometry correction
    """
    if not 0 < alpha < 1:
        raise ValueError(f"Fractional order α must be in (0,1), got {alpha}")

    # Handle zero eigenvalue (l=0 case)
    result = np.zeros_like(eigenvalues)
    nonzero_mask = eigenvalues > 0
    result[nonzero_mask] = eigenvalues[nonzero_mask] ** alpha

    return result


def compute_curvature_correction(l: int, R: float, alpha: float) -> float:
    """
    Curvature correction term from paper (lines 306-322).

    C(ℓ,R,α) = α(α-1)/12R² · ℓ(ℓ+1) + O(1/R⁴)

    This is the leading-order correction due to sphere curvature.
    Note: Correction is negative for 0<α<1 (reduces eigenvalue).

    Args:
        l: Angular momentum quantum number
        R: Sphere radius
        alpha: Fractional order

    Returns:
        C(ℓ,R,α): Correction term
    """
    if l < 0:
        raise ValueError(f"Angular momentum must be non-negative, got {l}")
    if R <= 0:
        raise ValueError(f"Radius must be positive, got {R}")
    if not 0 < alpha < 1:
        raise ValueError(f"Fractional order α must be in (0,1), got {alpha}")

    # C = α(α-1)/12 · (1/R²) · ℓ(ℓ+1)
    prefactor = alpha * (alpha - 1) / 12.0
    curvature_term = 1.0 / (R ** 2)
    angular_term = l * (l + 1)

    return prefactor * curvature_term * angular_term


def compute_fractional_spectrum_corrected(R: float, alpha: float, l_max: int) -> np.ndarray:
    """
    🧠 Function: compute_fractional_spectrum_corrected
    Role: Apply correction formula to get improved approximation
    Inputs: R (radius), α (order), l_max
    Returns: Array of corrected eigenvalues
    Notes: Main theoretical prediction being tested
    """
    # First compute exact fractional eigenvalues
    exact_eigenvalues = compute_exact_spectrum_sphere(R, l_max)
    fractional_exact = compute_fractional_spectrum_exact(exact_eigenvalues, alpha)

    # Then add corrections
    corrections = np.array([
        compute_curvature_correction(l, R, alpha) for l in range(l_max + 1)
    ])

    return fractional_exact + corrections


def compute_weak_curvature_parameter(R: float, l: int) -> float:
    """
    Compute curvature parameter κ_α = 1/R² · ℓ(ℓ+1).

    Regimes:
    - Weak curvature: κ_α < 0.1
    - Intermediate: 0.1 ≤ κ_α ≤ 1.0
    - Strong curvature: κ_α > 1.0

    Theory predicts correction formula valid in weak regime.

    Args:
        R: Sphere radius
        l: Angular momentum

    Returns:
        κ_α: Curvature parameter
    """
    if R <= 0:
        raise ValueError(f"Radius must be positive, got {R}")
    if l < 0:
        raise ValueError(f"Angular momentum must be non-negative, got {l}")

    return l * (l + 1) / (R ** 2)


def compute_convergence_analysis(
    R_values: List[float],
    alpha: float,
    l_fixed: int = 10
) -> Dict:
    """
    Analyze convergence with varying sphere radius.

    Theory predicts O(1/R²) convergence to flat space limit.

    Args:
        R_values: List of sphere radii to test
        alpha: Fractional order
        l_fixed: Fixed angular momentum for comparison

    Returns:
        Dictionary with convergence analysis
    """
    results = {
        'R_values': R_values,
        'alpha': alpha,
        'l_fixed': l_fixed,
        'corrections': [],
        'relative_corrections': [],
        'curvature_params': []
    }

    for R in R_values:
        # Compute correction magnitude
        correction = compute_curvature_correction(l_fixed, R, alpha)
        results['corrections'].append(correction)

        # Compute relative correction (compared to exact fractional)
        exact_eigenvalue = l_fixed * (l_fixed + 1) / (R ** 2)
        fractional_exact = exact_eigenvalue ** alpha if exact_eigenvalue > 0 else 0
        if abs(fractional_exact) > 1e-15:
            relative = abs(correction) / abs(fractional_exact)
        else:
            relative = float('inf')
        results['relative_corrections'].append(relative)

        # Compute curvature parameter
        kappa = compute_weak_curvature_parameter(R, l_fixed)
        results['curvature_params'].append(kappa)

    # Analyze scaling
    log_R = np.log(R_values)
    log_corrections = np.log(np.abs(results['corrections']))

    # Fit power law: correction ~ R^(-p)
    # log(correction) = -p * log(R) + const
    if len(R_values) > 1:
        # Linear fit in log-log space
        coeffs = np.polyfit(log_R, log_corrections, 1)
        power_law_exponent = -coeffs[0]  # Negative of slope
        results['scaling'] = {
            'fitted_exponent': power_law_exponent,
            'theoretical_exponent': 2.0,  # O(1/R²)
            'relative_error': abs(power_law_exponent - 2.0) / 2.0
        }
    else:
        results['scaling'] = None

    return results


def compute_spectrum_comparison(
    laplacian: SphericalLaplacian,
    l_range: Optional[Tuple[int, int]] = None
) -> Dict:
    """
    🧠 Function: compute_spectrum_comparison
    Role: Comprehensive comparison of exact vs corrected spectra
    Inputs: laplacian object, optional l_range
    Returns: Dictionary with detailed comparison
    Notes: Main analysis function for validation
    """
    if l_range is None:
        l_start, l_end = 0, laplacian.l_max
    else:
        l_start, l_end = l_range

    l_values = np.arange(l_start, l_end + 1)

    # Compute all spectra
    exact = []
    fractional_exact = []
    fractional_corrected = []
    corrections = []
    relative_errors = []
    curvature_params = []

    for l in l_values:
        # Exact eigenvalue
        lambda_l = laplacian.eigenvalue_exact(l)
        exact.append(lambda_l)

        # Fractional (no correction)
        lambda_alpha = laplacian.eigenvalue_fractional_exact(l)
        fractional_exact.append(lambda_alpha)

        # Fractional (with correction)
        lambda_corrected = laplacian.eigenvalue_fractional_corrected(l)
        fractional_corrected.append(lambda_corrected)

        # Correction term
        correction = laplacian.compute_correction_term(l)
        corrections.append(correction)

        # Relative error
        if abs(lambda_alpha) > 1e-15:
            rel_error = abs(lambda_corrected - lambda_alpha) / abs(lambda_alpha)
        else:
            rel_error = 0.0 if abs(lambda_corrected) < 1e-15 else float('inf')
        relative_errors.append(rel_error)

        # Curvature parameter
        kappa = laplacian.compute_weak_curvature_parameter(l)
        curvature_params.append(kappa)

    # Convert to arrays
    exact = np.array(exact)
    fractional_exact = np.array(fractional_exact)
    fractional_corrected = np.array(fractional_corrected)
    corrections = np.array(corrections)
    relative_errors = np.array(relative_errors)
    curvature_params = np.array(curvature_params)

    # Identify regimes
    weak_mask = curvature_params < 0.1
    strong_mask = curvature_params > 1.0
    intermediate_mask = ~weak_mask & ~strong_mask

    return {
        'l_values': l_values.tolist(),
        'eigenvalues': {
            'exact': exact.tolist(),
            'fractional_exact': fractional_exact.tolist(),
            'fractional_corrected': fractional_corrected.tolist()
        },
        'corrections': corrections.tolist(),
        'relative_errors': relative_errors.tolist(),
        'curvature_params': curvature_params.tolist(),
        'regime_analysis': {
            'weak': {
                'indices': np.where(weak_mask)[0].tolist(),
                'max_error': float(np.max(relative_errors[weak_mask])) if np.any(weak_mask) else None,
                'mean_error': float(np.mean(relative_errors[weak_mask])) if np.any(weak_mask) else None
            },
            'intermediate': {
                'indices': np.where(intermediate_mask)[0].tolist(),
                'max_error': float(np.max(relative_errors[intermediate_mask])) if np.any(intermediate_mask) else None,
                'mean_error': float(np.mean(relative_errors[intermediate_mask])) if np.any(intermediate_mask) else None
            },
            'strong': {
                'indices': np.where(strong_mask)[0].tolist(),
                'max_error': float(np.max(relative_errors[strong_mask])) if np.any(strong_mask) else None,
                'mean_error': float(np.mean(relative_errors[strong_mask])) if np.any(strong_mask) else None
            }
        },
        'statistics': {
            'max_correction': float(np.max(np.abs(corrections))),
            'mean_correction': float(np.mean(np.abs(corrections))),
            'max_relative_error': float(np.max(relative_errors[1:])) if len(relative_errors) > 1 else 0.0,  # Exclude l=0
            'mean_relative_error': float(np.mean(relative_errors[1:])) if len(relative_errors) > 1 else 0.0
        }
    }


def compute_parameter_sweep(
    R_values: List[float],
    alpha_values: List[float],
    l_test: int = 10
) -> Dict:
    """
    Sweep over multiple radii and fractional orders.

    Args:
        R_values: List of sphere radii
        alpha_values: List of fractional orders
        l_test: Angular momentum for testing

    Returns:
        Dictionary with sweep results
    """
    results = {
        'R_values': R_values,
        'alpha_values': alpha_values,
        'l_test': l_test,
        'sweep_data': []
    }

    for R in R_values:
        for alpha in alpha_values:
            # Create Laplacian
            laplacian = SphericalLaplacian(R=R, alpha=alpha, l_max=l_test + 10)

            # Compute key quantities for test l
            exact = laplacian.eigenvalue_fractional_exact(l_test)
            corrected = laplacian.eigenvalue_fractional_corrected(l_test)
            correction = laplacian.compute_correction_term(l_test)
            rel_error = laplacian.relative_error(l_test)
            kappa = laplacian.compute_weak_curvature_parameter(l_test)

            results['sweep_data'].append({
                'R': R,
                'alpha': alpha,
                'eigenvalue_exact': exact,
                'eigenvalue_corrected': corrected,
                'correction': correction,
                'relative_error': rel_error,
                'curvature_param': kappa,
                'regime': 'weak' if kappa < 0.1 else ('strong' if kappa > 1.0 else 'intermediate')
            })

    return results