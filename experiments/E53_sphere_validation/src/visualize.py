"""
📄 File: visualize.py
Purpose: Visualization of eigenvalue corrections and error analysis
Created: November 11, 2025
Used by: main.py
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from typing import Dict, List, Optional, Tuple
from .core import SphericalLaplacian
from .compute import compute_spectrum_comparison, compute_convergence_analysis


def plot_error_vs_l(
    laplacian: SphericalLaplacian,
    save_path: Optional[str] = None,
    show_regimes: bool = True
) -> plt.Figure:
    """
    Plot relative error vs angular momentum ℓ.

    Shows how the correction accuracy varies with ℓ and highlights
    different curvature regimes (weak/intermediate/strong).

    Args:
        laplacian: SphericalLaplacian object
        save_path: Path to save figure (optional)
        show_regimes: Whether to show regime boundaries

    Returns:
        matplotlib Figure
    """
    # Compute errors
    l_values = np.arange(1, laplacian.l_max + 1)  # Skip l=0
    errors = np.array([laplacian.relative_error(l) for l in l_values])
    curvature_params = np.array([
        laplacian.compute_weak_curvature_parameter(l) for l in l_values
    ])

    # Create figure
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    # Top panel: Relative error
    ax1.semilogy(l_values, errors, 'b.-', label='Relative Error', markersize=4)
    ax1.axhline(y=0.01, color='r', linestyle='--', label='1% Threshold', alpha=0.7)
    ax1.set_ylabel('Relative Error', fontsize=12)
    ax1.legend(loc='best')
    ax1.grid(True, alpha=0.3)
    ax1.set_title(f'Fractional Laplacian on S² (R={laplacian.R}, α={laplacian.alpha})', fontsize=14)

    # Add regime shading if requested
    if show_regimes:
        # Find transition points
        l_weak_boundary = int(np.sqrt(0.1 * laplacian.R**2))
        l_strong_boundary = int(np.sqrt(1.0 * laplacian.R**2))

        # Shade regions
        if l_weak_boundary > 0:
            ax1.axvspan(0, l_weak_boundary, alpha=0.2, color='green', label='Weak κ<0.1')
        if l_weak_boundary < l_strong_boundary and l_strong_boundary <= laplacian.l_max:
            ax1.axvspan(l_weak_boundary, l_strong_boundary, alpha=0.2, color='orange', label='Intermediate')
        if l_strong_boundary <= laplacian.l_max:
            ax1.axvspan(l_strong_boundary, laplacian.l_max, alpha=0.2, color='red', label='Strong κ>1')

    # Bottom panel: Curvature parameter
    ax2.semilogy(l_values, curvature_params, 'g.-', markersize=4)
    ax2.axhline(y=0.1, color='b', linestyle='--', label='κ=0.1 (weak boundary)', alpha=0.7)
    ax2.axhline(y=1.0, color='r', linestyle='--', label='κ=1.0 (strong boundary)', alpha=0.7)
    ax2.set_xlabel('Angular Momentum ℓ', fontsize=12)
    ax2.set_ylabel('Curvature Parameter κ_α', fontsize=12)
    ax2.legend(loc='best')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Saved error plot to {save_path}")

    return fig


def plot_correction_magnitude(
    laplacian: SphericalLaplacian,
    save_path: Optional[str] = None
) -> plt.Figure:
    """
    🧠 Function: plot_correction_magnitude
    Role: Visualize correction term magnitude and its scaling
    Inputs: laplacian, optional save_path
    Returns: matplotlib Figure
    Notes: Shows both absolute and relative correction magnitudes
    """
    # Compute corrections
    l_values = np.arange(1, laplacian.l_max + 1)
    corrections = np.array([laplacian.compute_correction_term(l) for l in l_values])
    exact_fractional = np.array([laplacian.eigenvalue_fractional_exact(l) for l in l_values])

    # Relative corrections
    relative_corrections = np.abs(corrections / exact_fractional)

    # Create figure
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

    # Top: Absolute correction
    ax1.plot(l_values, -corrections, 'b.-', markersize=4)  # Negative because corrections are negative
    ax1.set_xlabel('Angular Momentum ℓ', fontsize=12)
    ax1.set_ylabel('|Correction Term|', fontsize=12)
    ax1.set_title(f'Curvature Correction Magnitude (R={laplacian.R}, α={laplacian.alpha})', fontsize=14)
    ax1.grid(True, alpha=0.3)
    ax1.text(0.05, 0.95, f'C(ℓ,R,α) = α(α-1)/12R² · ℓ(ℓ+1)',
            transform=ax1.transAxes, fontsize=10,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="yellow", alpha=0.3))

    # Bottom: Relative correction
    ax2.semilogy(l_values, relative_corrections, 'r.-', markersize=4)
    ax2.set_xlabel('Angular Momentum ℓ', fontsize=12)
    ax2.set_ylabel('|Correction| / |Exact|', fontsize=12)
    ax2.set_title('Relative Correction Magnitude', fontsize=14)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Saved correction plot to {save_path}")

    return fig


def plot_convergence_analysis(
    R_values: List[float],
    alpha: float,
    l_fixed: int = 10,
    save_path: Optional[str] = None
) -> plt.Figure:
    """
    Plot convergence analysis with varying radius.

    Shows O(1/R²) scaling of corrections.

    Args:
        R_values: List of sphere radii
        alpha: Fractional order
        l_fixed: Fixed angular momentum
        save_path: Path to save figure

    Returns:
        matplotlib Figure
    """
    # Compute convergence data
    results = compute_convergence_analysis(R_values, alpha, l_fixed)

    # Create figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Left: Correction vs R
    ax1.loglog(R_values, np.abs(results['corrections']), 'bo-', label='Actual', markersize=8)

    # Add theoretical scaling
    R_array = np.array(R_values)
    theoretical = results['corrections'][0] * (R_values[0] / R_array)**2
    ax1.loglog(R_array, np.abs(theoretical), 'r--', label='O(1/R²) Theory', alpha=0.7)

    ax1.set_xlabel('Radius R', fontsize=12)
    ax1.set_ylabel('|Correction|', fontsize=12)
    ax1.set_title(f'Scaling of Correction Term (ℓ={l_fixed}, α={alpha})', fontsize=14)
    ax1.legend()
    ax1.grid(True, alpha=0.3, which='both')

    # Right: Relative correction vs curvature parameter
    ax2.loglog(results['curvature_params'], results['relative_corrections'], 'go-', markersize=8)
    ax2.axvline(x=0.1, color='b', linestyle='--', alpha=0.5, label='Weak regime boundary')
    ax2.axvline(x=1.0, color='r', linestyle='--', alpha=0.5, label='Strong regime boundary')
    ax2.set_xlabel('Curvature Parameter κ_α', fontsize=12)
    ax2.set_ylabel('Relative Correction', fontsize=12)
    ax2.set_title('Relative Correction vs Curvature', fontsize=14)
    ax2.legend()
    ax2.grid(True, alpha=0.3, which='both')

    # Add scaling exponent info if available
    if results.get('scaling'):
        scaling_info = results['scaling']
        text = f"Fitted exponent: {scaling_info['fitted_exponent']:.2f}\n"
        text += f"Theory exponent: {scaling_info['theoretical_exponent']:.1f}\n"
        text += f"Relative error: {scaling_info['relative_error']:.1%}"
        ax1.text(0.05, 0.05, text, transform=ax1.transAxes,
                fontsize=10, bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Saved convergence plot to {save_path}")

    return fig


def plot_parameter_sweep(
    R_values: List[float],
    alpha_values: List[float],
    l_test: int = 10,
    save_path: Optional[str] = None
) -> plt.Figure:
    """
    🧠 Function: plot_parameter_sweep
    Role: 2D heatmap of errors across parameter space
    Inputs: R_values, alpha_values, l_test, save_path
    Returns: matplotlib Figure with heatmap
    Notes: Useful for identifying optimal parameter regimes
    """
    # Create meshgrid
    R_mesh, alpha_mesh = np.meshgrid(R_values, alpha_values)
    error_mesh = np.zeros_like(R_mesh)

    # Compute errors
    for i, alpha in enumerate(alpha_values):
        for j, R in enumerate(R_values):
            laplacian = SphericalLaplacian(R=R, alpha=alpha, l_max=l_test + 5)
            error_mesh[i, j] = laplacian.relative_error(l_test)

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 8))

    # Heatmap with log scale
    im = ax.contourf(R_mesh, alpha_mesh, np.log10(error_mesh + 1e-10),
                     levels=20, cmap='viridis')
    contours = ax.contour(R_mesh, alpha_mesh, np.log10(error_mesh + 1e-10),
                          levels=[-3, -2, -1], colors='white', linewidths=2)
    ax.clabel(contours, inline=True, fontsize=10, fmt='10^%d')

    # Add 1% error boundary
    error_1pct = ax.contour(R_mesh, alpha_mesh, error_mesh,
                           levels=[0.01], colors='red', linewidths=3,
                           linestyles='--', label='1% Error')

    # Labels and formatting
    ax.set_xlabel('Radius R', fontsize=12)
    ax.set_ylabel('Fractional Order α', fontsize=12)
    ax.set_title(f'Relative Error Heatmap (ℓ={l_test})', fontsize=14)

    # Colorbar
    cbar = plt.colorbar(im, ax=ax, label='log₁₀(Relative Error)')

    # Add text annotation for regimes
    kappa = l_test * (l_test + 1) / R_mesh**2
    # Mark weak regime boundary
    R_weak = np.sqrt(l_test * (l_test + 1) / 0.1)
    ax.axvline(x=R_weak, color='white', linestyle=':', alpha=0.5)
    ax.text(R_weak * 1.1, 0.9, 'κ=0.1', fontsize=10, color='white')

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Saved parameter sweep to {save_path}")

    return fig


def plot_eigenvalue_comparison(
    laplacian: SphericalLaplacian,
    l_range: Tuple[int, int] = (1, 20),
    save_path: Optional[str] = None
) -> plt.Figure:
    """
    Compare exact vs corrected fractional eigenvalues.

    Args:
        laplacian: SphericalLaplacian object
        l_range: Range of l values to plot
        save_path: Path to save figure

    Returns:
        matplotlib Figure
    """
    # Get comparison data
    comparison = compute_spectrum_comparison(laplacian, l_range)

    l_values = np.array(comparison['l_values'])
    exact = np.array(comparison['eigenvalues']['fractional_exact'])
    corrected = np.array(comparison['eigenvalues']['fractional_corrected'])
    corrections = np.array(comparison['corrections'])

    # Create figure
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    # Top: Eigenvalues
    ax1.plot(l_values, exact, 'b.-', label='Exact λ^α', markersize=6)
    ax1.plot(l_values, corrected, 'r.--', label='Corrected λ^α + C', markersize=6, alpha=0.7)
    ax1.set_ylabel('Eigenvalue', fontsize=12)
    ax1.set_title(f'Fractional Eigenvalue Comparison (R={laplacian.R}, α={laplacian.alpha})',
                 fontsize=14)
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Bottom: Difference
    ax2.plot(l_values, corrections, 'g.-', markersize=6)
    ax2.set_xlabel('Angular Momentum ℓ', fontsize=12)
    ax2.set_ylabel('Correction C(ℓ,R,α)', fontsize=12)
    ax2.set_title('Curvature Correction Term', fontsize=14)
    ax2.grid(True, alpha=0.3)
    ax2.axhline(y=0, color='black', linestyle='-', alpha=0.3)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Saved eigenvalue comparison to {save_path}")

    return fig


def create_summary_figure(
    laplacian: SphericalLaplacian,
    save_path: Optional[str] = None
) -> plt.Figure:
    """
    🧠 Function: create_summary_figure
    Role: Create comprehensive 4-panel summary figure
    Inputs: laplacian, save_path
    Returns: matplotlib Figure with 4 subplots
    Notes: Main visualization for paper/presentation
    """
    # Compute all needed data
    l_values = np.arange(1, laplacian.l_max + 1)
    errors = np.array([laplacian.relative_error(l) for l in l_values])
    corrections = np.array([laplacian.compute_correction_term(l) for l in l_values])
    exact_fractional = np.array([laplacian.eigenvalue_fractional_exact(l) for l in l_values])
    corrected_fractional = np.array([laplacian.eigenvalue_fractional_corrected(l) for l in l_values])
    curvature_params = np.array([laplacian.compute_weak_curvature_parameter(l) for l in l_values])

    # Create 2x2 subplot figure
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 10))

    # Panel 1: Error vs l with regime shading
    ax1.semilogy(l_values, errors, 'b.-', markersize=4)
    ax1.axhline(y=0.01, color='r', linestyle='--', label='1% Threshold', alpha=0.7)

    # Shade regimes
    weak_mask = curvature_params < 0.1
    if np.any(weak_mask):
        l_weak = l_values[weak_mask]
        ax1.axvspan(l_weak.min(), l_weak.max(), alpha=0.2, color='green', label='Weak κ<0.1')

    ax1.set_xlabel('Angular Momentum ℓ')
    ax1.set_ylabel('Relative Error')
    ax1.set_title('Error Analysis')
    ax1.legend(loc='best', fontsize=9)
    ax1.grid(True, alpha=0.3)

    # Panel 2: Eigenvalue comparison
    ax2.plot(l_values[:20], exact_fractional[:20], 'b.-', label='Exact λ^α', markersize=6)
    ax2.plot(l_values[:20], corrected_fractional[:20], 'r.--', label='Corrected', markersize=6, alpha=0.7)
    ax2.set_xlabel('Angular Momentum ℓ')
    ax2.set_ylabel('Fractional Eigenvalue')
    ax2.set_title('Eigenvalue Comparison')
    ax2.legend(loc='best', fontsize=9)
    ax2.grid(True, alpha=0.3)

    # Panel 3: Correction scaling
    ax3.loglog(l_values, -corrections, 'g.-', markersize=4)
    # Add power law fit
    coeffs = np.polyfit(np.log(l_values[5:]), np.log(-corrections[5:]), 1)
    fitted = np.exp(np.polyval(coeffs, np.log(l_values)))
    ax3.loglog(l_values, fitted, 'r--', alpha=0.5,
              label=f'Fit: ℓ^{coeffs[0]:.2f}')
    ax3.set_xlabel('Angular Momentum ℓ')
    ax3.set_ylabel('|Correction|')
    ax3.set_title('Correction Scaling')
    ax3.legend(loc='best', fontsize=9)
    ax3.grid(True, alpha=0.3, which='both')

    # Panel 4: Curvature parameter
    ax4.semilogy(l_values, curvature_params, 'purple', linewidth=2)
    ax4.axhline(y=0.1, color='b', linestyle='--', alpha=0.5, label='Weak boundary')
    ax4.axhline(y=1.0, color='r', linestyle='--', alpha=0.5, label='Strong boundary')
    ax4.fill_between(l_values, 0.001, 0.1, alpha=0.2, color='green', label='Weak regime')
    ax4.fill_between(l_values, 0.1, 1.0, alpha=0.2, color='orange', label='Intermediate')
    ax4.fill_between(l_values, 1.0, 100, alpha=0.2, color='red', label='Strong regime')
    ax4.set_xlabel('Angular Momentum ℓ')
    ax4.set_ylabel('κ_α = ℓ(ℓ+1)/R²')
    ax4.set_title('Curvature Regimes')
    ax4.legend(loc='best', fontsize=9)
    ax4.grid(True, alpha=0.3)
    ax4.set_ylim([0.01, 100])

    # Overall title
    fig.suptitle(f'Spherical Fractional Laplacian Analysis (R={laplacian.R}, α={laplacian.alpha})',
                fontsize=16, y=1.02)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Saved summary figure to {save_path}")

    return fig