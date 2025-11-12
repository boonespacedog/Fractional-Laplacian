"""
📄 File: core.py
Purpose: Core objects for fractional Laplacian on sphere S²
Created: November 11, 2025
Used by: compute.py, validate.py, visualize.py, main.py
"""

import numpy as np
from typing import Optional, Tuple, Dict
import warnings


class SphericalLaplacian:
    """
    Fractional Laplacian on 2-sphere with curvature corrections.

    This class represents the fractional Laplace-Beltrami operator (-Δ_S²)^α
    on the unit 2-sphere, including analytical eigenvalues and first-order
    curvature corrections from the paper.

    Theory reference: Paper Section 3.2, lines 306-322
    Mathematical basis: Spherical harmonic eigenvalues ℓ(ℓ+1)/R²
    """

    def __init__(self, R: float = 1.0, alpha: float = 0.5, l_max: int = 50):
        """
        Initialize spherical Laplacian with given parameters.

        Args:
            R: Sphere radius (default: 1.0)
            alpha: Fractional order ∈ (0,1) (default: 0.5)
            l_max: Maximum angular momentum quantum number (default: 50)

        Raises:
            ValueError: If R <= 0, alpha not in (0,1), or l_max < 0
        """
        # === Input Validation ===
        if R <= 0:
            raise ValueError(f"Radius R must be positive, got {R}")
        if not 0 < alpha < 1:
            raise ValueError(f"Fractional order α must be in (0,1), got {alpha}")
        if l_max < 0:
            raise ValueError(f"l_max must be non-negative, got {l_max}")

        # Store parameters
        self.R = R
        self.alpha = alpha
        self.l_max = l_max

        # Lazy-computed properties (cached after first access)
        self._eigenvalues_exact = None
        self._eigenvalues_fractional_exact = None
        self._eigenvalues_fractional_corrected = None
        self._correction_terms = None

    def eigenvalue_exact(self, l: int) -> float:
        """
        Compute exact eigenvalue of Laplace-Beltrami operator on S².

        Formula: λ_ℓ = ℓ(ℓ+1)/R²

        This is the textbook formula for spherical harmonics on a sphere
        of radius R. Each eigenvalue has degeneracy (2ℓ+1).

        Args:
            l: Angular momentum quantum number (ℓ ≥ 0)

        Returns:
            λ_ℓ: Eigenvalue for angular momentum ℓ

        Notes:
            - This is the EXACT analytical formula (no approximation)
            - Eigenvalues are independent of magnetic quantum number m
            - Each λ_ℓ has multiplicity 2ℓ+1 (from m = -ℓ,...,+ℓ)
        """
        if l < 0:
            raise ValueError(f"Angular momentum ℓ must be non-negative, got {l}")

        return float(l * (l + 1)) / (self.R ** 2)

    def eigenvalue_fractional_exact(self, l: int) -> float:
        """
        🧠 Function: eigenvalue_fractional_exact
        Role: Compute exact fractional power of eigenvalue
        Inputs: l (angular momentum)
        Returns: λ_ℓ^α = [ℓ(ℓ+1)/R²]^α
        Notes: This is the "naive" fractional power without corrections
        """
        lambda_l = self.eigenvalue_exact(l)

        # Handle l=0 case (zero eigenvalue)
        if lambda_l == 0:
            return 0.0

        return lambda_l ** self.alpha

    def compute_correction_term(self, l: int) -> float:
        """
        Compute first-order curvature correction.

        Formula from paper (lines 306-322):
        C(ℓ,R,α) = α(α-1)/12 · (1/R²) · ℓ(ℓ+1)

        This correction arises from the non-flat geometry of the sphere
        and becomes significant when κ_α = ℓ(ℓ+1)/R² is not small.

        Args:
            l: Angular momentum quantum number

        Returns:
            C(ℓ,R,α): Curvature correction term

        Notes:
            - Correction is O(1/R²) in the large-R limit
            - For α=0 or α=1, correction vanishes (as expected)
            - Negative for 0<α<1 (reduces eigenvalue magnitude)
        """
        if l < 0:
            raise ValueError(f"Angular momentum ℓ must be non-negative, got {l}")

        # Compute correction: α(α-1)/12 · (1/R²) · ℓ(ℓ+1)
        prefactor = self.alpha * (self.alpha - 1) / 12.0
        curvature_factor = 1.0 / (self.R ** 2)
        angular_factor = float(l * (l + 1))

        return prefactor * curvature_factor * angular_factor

    def eigenvalue_fractional_corrected(self, l: int) -> float:
        """
        🧠 Function: eigenvalue_fractional_corrected
        Role: Apply curvature correction to fractional eigenvalue
        Inputs: l (angular momentum)
        Returns: λ_ℓ^α + C(ℓ,R,α)
        Notes: Main theoretical prediction we're testing
        """
        exact = self.eigenvalue_fractional_exact(l)
        correction = self.compute_correction_term(l)
        return exact + correction

    def relative_error(self, l: int) -> float:
        """
        Compute relative error between exact and corrected eigenvalues.

        Formula: |λ_corrected - λ_exact| / |λ_exact|

        Theory predicts this should be < 1% in weak curvature regime.

        Args:
            l: Angular momentum quantum number

        Returns:
            Relative error (dimensionless)
        """
        exact = self.eigenvalue_fractional_exact(l)
        corrected = self.eigenvalue_fractional_corrected(l)

        # Handle zero eigenvalue case
        if abs(exact) < 1e-15:
            if abs(corrected) < 1e-15:
                return 0.0  # Both zero
            else:
                return float('inf')  # Undefined relative error

        return abs(corrected - exact) / abs(exact)

    def compute_weak_curvature_parameter(self, l: int) -> float:
        """
        🧠 Function: compute_weak_curvature_parameter
        Role: Check if we're in weak curvature regime
        Inputs: l (angular momentum)
        Returns: κ_α = ℓ(ℓ+1)/R²
        Notes: Weak regime when κ_α < 0.1, strong when κ_α > 1
        """
        return float(l * (l + 1)) / (self.R ** 2)

    @property
    def eigenvalues_exact(self) -> np.ndarray:
        """
        Get array of exact eigenvalues (cached).

        Returns:
            Array of λ_ℓ for ℓ = 0,1,...,l_max
        """
        if self._eigenvalues_exact is None:
            self._eigenvalues_exact = np.array([
                self.eigenvalue_exact(l) for l in range(self.l_max + 1)
            ])
        return self._eigenvalues_exact

    @property
    def eigenvalues_fractional_exact(self) -> np.ndarray:
        """
        Get array of exact fractional eigenvalues (cached).

        Returns:
            Array of λ_ℓ^α for ℓ = 0,1,...,l_max
        """
        if self._eigenvalues_fractional_exact is None:
            self._eigenvalues_fractional_exact = np.array([
                self.eigenvalue_fractional_exact(l) for l in range(self.l_max + 1)
            ])
        return self._eigenvalues_fractional_exact

    @property
    def eigenvalues_fractional_corrected(self) -> np.ndarray:
        """
        Get array of corrected fractional eigenvalues (cached).

        Returns:
            Array of λ_ℓ^α + C(ℓ,R,α) for ℓ = 0,1,...,l_max
        """
        if self._eigenvalues_fractional_corrected is None:
            self._eigenvalues_fractional_corrected = np.array([
                self.eigenvalue_fractional_corrected(l) for l in range(self.l_max + 1)
            ])
        return self._eigenvalues_fractional_corrected

    @property
    def correction_terms(self) -> np.ndarray:
        """
        Get array of correction terms (cached).

        Returns:
            Array of C(ℓ,R,α) for ℓ = 0,1,...,l_max
        """
        if self._correction_terms is None:
            self._correction_terms = np.array([
                self.compute_correction_term(l) for l in range(self.l_max + 1)
            ])
        return self._correction_terms

    def get_spectrum_summary(self) -> Dict:
        """
        Get comprehensive summary of eigenvalue spectrum.

        Returns:
            Dictionary with all spectral information
        """
        l_values = np.arange(self.l_max + 1)
        relative_errors = np.array([
            self.relative_error(l) for l in range(self.l_max + 1)
        ])
        curvature_params = np.array([
            self.compute_weak_curvature_parameter(l) for l in range(self.l_max + 1)
        ])

        # Find regime boundaries
        weak_mask = curvature_params < 0.1
        strong_mask = curvature_params > 1.0
        intermediate_mask = ~weak_mask & ~strong_mask

        return {
            'parameters': {
                'R': self.R,
                'alpha': self.alpha,
                'l_max': self.l_max
            },
            'eigenvalues': {
                'exact': self.eigenvalues_exact.tolist(),
                'fractional_exact': self.eigenvalues_fractional_exact.tolist(),
                'fractional_corrected': self.eigenvalues_fractional_corrected.tolist(),
                'corrections': self.correction_terms.tolist()
            },
            'errors': {
                'relative': relative_errors.tolist(),
                'max_relative_error': float(np.max(relative_errors[1:])),  # Exclude l=0
                'mean_relative_error': float(np.mean(relative_errors[1:]))  # Exclude l=0
            },
            'regimes': {
                'weak_curvature': {
                    'l_values': l_values[weak_mask].tolist(),
                    'count': int(np.sum(weak_mask)),
                    'max_error': float(np.max(relative_errors[weak_mask])) if np.any(weak_mask) else None
                },
                'intermediate': {
                    'l_values': l_values[intermediate_mask].tolist(),
                    'count': int(np.sum(intermediate_mask)),
                    'max_error': float(np.max(relative_errors[intermediate_mask])) if np.any(intermediate_mask) else None
                },
                'strong_curvature': {
                    'l_values': l_values[strong_mask].tolist(),
                    'count': int(np.sum(strong_mask)),
                    'max_error': float(np.max(relative_errors[strong_mask])) if np.any(strong_mask) else None
                }
            },
            'theory_validation': {
                'weak_regime_error_bound': 0.01,  # Theory predicts < 1%
                'weak_regime_satisfied': bool(np.all(relative_errors[weak_mask] < 0.01)) if np.any(weak_mask) else None,
                'convergence_order': 'O(κ_α²)'  # Theoretical prediction
            }
        }