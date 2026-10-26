# SPDX-License-Identifier: MIT
r"""
Matrix-valued Mittag-Leffler function evaluation.

Computes E_{α,β}^γ(A) for A ∈ ℂ^{n×n} via three methods:

* ``"eigen"``   — eigendecomposition A = VΛV⁻¹ (fast, fails if non-diag)
* ``"schur"``   — Schur-Parlett recurrence (robust, default)
* ``"contour"`` — Cauchy integral over OPC contour (large n, slow)

The default ``"auto"`` tries eigen, falls back to Schur, then to contour.

Notes
-----
* All three methods agree to machine precision for well-conditioned A.
* For poorly conditioned A (||V||·||V⁻¹|| > 1e12), Schur-Parlett is
  preferred; see Higham (2008), Functions of Matrices, Ch. 4.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import scipy.linalg as sla
from numpy.typing import ArrayLike, NDArray

from mlfp.mlfp import (
    MittagLeffler,
    MittagLefflerError,
    _safe_z_limit,
)


# ============================================================================
# Public class
# ============================================================================
@dataclass(slots=True)
class MittagLefflerMatrix:
    r"""
    Matrix-valued Mittag-Leffler evaluator:  E_{α,β}^γ(A).

    Parameters
    ----------
    alpha : float  (α > 0)
    beta  : float  (default 1.0)
    gamma : float  (default 1.0, Prabhakar parameter)
    method : {"auto", "eigen", "schur", "contour"}  (default "auto")
    tolerance : float  (relative accuracy, default 1e-13)
    block_tol : float  (Schur-Parlett clustering tolerance, default 1e-10)

    Examples
    --------
    >>> import numpy as np
    >>> A = np.array([[-1.0, 0.5], [0.5, -2.0]])
    >>> E = MittagLefflerMatrix(alpha=0.5)(A)
    >>> E.shape
    (2, 2)

    >>> # E_{1,1}(A) = exp(A) — sanity check
    >>> A = np.array([[0.1, 0.2], [0.3, 0.4]])
    >>> E = MittagLefflerMatrix(alpha=1.0, beta=1.0)(A)
    >>> assert np.allclose(E, sla.expm(A), rtol=1e-12)
    """

    alpha: float
    beta: float = 1.0
    gamma: float = 1.0
    method: str = "auto"
    tolerance: float = 1e-13
    block_tol: float = 1e-10

    def __post_init__(self) -> None:
        if not (math.isfinite(self.alpha) and self.alpha > 0.0):
            raise MittagLefflerError("alpha must be positive real.")
        if not math.isfinite(self.beta):
            raise MittagLefflerError("beta must be finite real.")
        if not (math.isfinite(self.gamma) and self.gamma > 0.0):
            raise MittagLefflerError("gamma must be positive real.")
        if self.method not in ("auto", "eigen", "schur", "contour"):
            raise MittagLefflerError(
                f"method must be one of auto/eigen/schur/contour; "
                f"got {self.method!r}"
            )
        if not (0.0 < self.tolerance < 1.0):
            raise MittagLefflerError("tolerance must satisfy 0 < tol < 1.")

    # ------------------------------------------------------------------
    def __call__(self, A: ArrayLike) -> NDArray[np.complex128]:
        """
        Evaluate E^γ_{α,β}(A).

        Parameters
        ----------
        A : array_like, shape (n, n) or (N, n, n)
            Square matrix or batch of square matrices (complex or real).

        Returns
        -------
        ndarray of complex128, same shape as A
        """
        A_arr = np.asarray(A, dtype=np.complex128)
        if A_arr.ndim == 2:
            if A_arr.shape[0] != A_arr.shape[1]:
                raise MittagLefflerError(
                    f"A must be square, got {A_arr.shape}"
                )
            return self._eval_single(A_arr)
        elif A_arr.ndim == 3:
            N, n, m = A_arr.shape
            if n != m:
                raise MittagLefflerError(
                    f"each matrix must be square, got {n}×{m}"
                )
            out = np.empty_like(A_arr)
            for i in range(N):
                out[i] = self._eval_single(A_arr[i])
            return out
        else:
            raise MittagLefflerError(
                f"A must be 2D or 3D, got shape {A_arr.shape}"
            )

    def _mp_series_at(self, z_mp, alpha, beta, gamma):
        """mpmath series for E^γ_{α,β}(z) at current dps (no nested workdps)."""
        import mpmath as _mp
        a_mp = _mp.mpf(alpha)
        b_mp = _mp.mpf(beta)
        g_mp = _mp.mpf(gamma)

        total = _mp.rgamma(b_mp)
        z_power = _mp.mpc(1)
        poch = _mp.mpf(1)
        tol = _mp.mpf(10) ** (-(_mp.mp.dps - 8))

        for k in range(1, 200_000):
            z_power *= z_mp
            poch *= g_mp + k - 1
            term = (
                poch * z_power * _mp.rgamma(a_mp * k + b_mp)
                / _mp.factorial(k)
            )
            total += term
            if term != 0 and abs(term) < abs(total) * tol:
                return total

        raise MittagLefflerError(
            "mpmath series did not converge within 200k terms."
        )

    def _scalar_ml_hp(self, z: complex, dps: int = 25) -> complex:
        """High-precision scalar ML at z (uses mpmath)."""
        import mpmath as _mp
        with _mp.workdps(dps):
            return complex(self._mp_series_at(
                _mp.mpc(z), self.alpha, self.beta, self.gamma,
            ))


    # ------------------------------------------------------------------
    # Dispatcher
    # ------------------------------------------------------------------
    def _eval_single(self, A: np.ndarray) -> np.ndarray:
        method = self.method
        if method == "eigen":
            return self._via_eigen(A)
        if method == "schur":
            return self._via_schur(A)
        if method == "contour":
            return self._via_contour(A)
        # auto
        try:
            return self._via_eigen(A)
        except (np.linalg.LinAlgError, MittagLefflerError):
            pass
        try:
            return self._via_schur(A)
        except (np.linalg.LinAlgError, MittagLefflerError):
            pass
        return self._via_contour(A)

    # ------------------------------------------------------------------
    # Scalar fallback — Mittag-Leffler on scalars
    # ------------------------------------------------------------------
    def _scalar_ml(self, z: complex) -> complex:
        """Evaluate E^γ_{α,β}(z) at a single scalar point."""
        ml = MittagLeffler(
            alpha=self.alpha,
            beta=self.beta,
            gamma=self.gamma,
            tolerance=self.tolerance,
        )
        return complex(ml.evaluate_complex(z))

    # ------------------------------------------------------------------
    # Method 1: eigendecomposition
    # ------------------------------------------------------------------
    def _via_eigen(self, A: np.ndarray) -> np.ndarray:
        """A = VΛV⁻¹ → f(A) = V f(Λ) V⁻¹."""
        n = A.shape[0]
        if n == 1:
            return np.array([[self._scalar_ml(A[0, 0])]], dtype=complex)

        w, V = np.linalg.eig(A)
        cond_V = np.linalg.cond(V)
        if cond_V > 1e12:
            raise MittagLefflerError(
                f"Eigendecomposition ill-conditioned (cond(V)={cond_V:.1e}); "
                "use method='schur'."
            )

        f_w = np.array([self._scalar_ml(wi) for wi in w], dtype=complex)
        return V @ np.diag(f_w) @ np.linalg.inv(V)

    # ------------------------------------------------------------------
    # Method 2: Schur-Parlett
    # ------------------------------------------------------------------
    def _via_schur(self, A: np.ndarray) -> np.ndarray:
        """
        Schur-Parlett algorithm.

        1. T, Q = schur(A)   →   A = QTQᴴ with T upper triangular
        2. Group near-equal diagonal entries into blocks
        3. F_diag = f(block) directly
        4. F_offdiag via Sylvester: T_ii X - X T_jj = R
        5. Return Q F Qᴴ
        """
        n = A.shape[0]
        if n == 1:
            return np.array([[self._scalar_ml(A[0, 0])]], dtype=complex)

        T, Q = sla.schur(A, output="complex")
        diag = np.diag(T).copy()

        # --- Block partition (eigenvalue clustering) ---
        blocks: list[tuple[int, int]] = []
        i = 0
        while i < n:
            j = i
            while (j < n - 1
                   and abs(diag[j + 1] - diag[i]) < self.block_tol):
                j += 1
            blocks.append((i, j))
            i = j + 1

        F = np.zeros_like(T)

        # --- Diagonal blocks: F_ii = f(T_ii) ---
        for (i0, i1) in blocks:
            block = T[i0:i1 + 1, i0:i1 + 1]
            F[i0:i1 + 1, i0:i1 + 1] = self._block_ml(block)

        # --- Off-diagonal blocks: Parlett recurrence ---
        n_b = len(blocks)
        for d in range(1, n_b):
            for i in range(n_b - d):
                j = i + d
                i0, i1 = blocks[i]
                j0, j1 = blocks[j]

                # R = F_ii T_ij − T_ij F_jj + Σ_{i<k<j}(F_ik T_kj − T_ik F_kj)
                R = (F[i0:i1 + 1, i0:i1 + 1] @ T[i0:i1 + 1, j0:j1 + 1]
                     - T[i0:i1 + 1, j0:j1 + 1] @ F[j0:j1 + 1, j0:j1 + 1])

                for k in range(i + 1, j):
                    k0, k1 = blocks[k]
                    R += (F[i0:i1 + 1, k0:k1 + 1] @ T[k0:k1 + 1, j0:j1 + 1]
                          - T[i0:i1 + 1, k0:k1 + 1] @ F[k0:k1 + 1, j0:j1 + 1])

                # Sylvester: T_ii X − X T_jj = R
                # scipy.solve_sylvester(A, B, Q) solves AX + XB = Q
                # So we use B = −T_jj
                F[i0:i1 + 1, j0:j1 + 1] = sla.solve_sylvester(
                    T[i0:i1 + 1, i0:i1 + 1],
                    -T[j0:j1 + 1, j0:j1 + 1],
                    R,
                )

        return Q @ F @ Q.conj().T

    def _block_ml(self, block: np.ndarray) -> np.ndarray:
        """
        E^γ_{α,β}(block) for a small square block.

        Uses Hermite interpolation around the block mean:
            F = Σ_{k=0}^{m-1} f^(k)(c)/k! · N^k,   N = block − cI
        """
        n = block.shape[0]
        if n == 1:
            return np.array([[self._scalar_ml(block[0, 0])]], dtype=complex)

        # Fast path: diagonalizable
        try:
            w, V = np.linalg.eig(block)
            if np.linalg.cond(V) < 1e10:
                f_w = np.array(
                    [self._scalar_ml(wi) for wi in w], dtype=complex,
                )
                return V @ np.diag(f_w) @ np.linalg.inv(V)
        except np.linalg.LinAlgError:
            pass

        # Hermite interpolation for defective blocks
        c = np.mean(np.diag(block))
        N = block - c * np.eye(n, dtype=complex)
        result = np.zeros((n, n), dtype=complex)
        for k in range(n):
            result += (
                self._scalar_ml_deriv_k(c, k)
                * np.linalg.matrix_power(N, k)
            )
        return result

    def _scalar_ml_deriv_k(self, c: complex, k: int) -> complex:
        """f^(k)(c)/k! where f = E^γ_{α,β}, via mpmath.diff."""
        import mpmath as _mp
        with _mp.workdps(30):
            c_mp = _mp.mpc(c)

            def f(z):
                return self._mp_series_at(
                    z, self.alpha, self.beta, self.gamma,
                )

            if k == 0:
                return complex(f(c_mp))
            d = _mp.diff(f, c_mp, k)
            return complex(d) / math.factorial(k)

    def _scalar_ml_complex_real(self, x: float, y: float) -> float:
        """Wrapper for scipy.differentiate (real argument only)."""
        return self._scalar_ml(complex(x, y)).real

    # ------------------------------------------------------------------
    # Method 3: contour integral (Cauchy)
    # ------------------------------------------------------------------
    def _via_contour(self, A: np.ndarray) -> np.ndarray:
        """
        Cauchy integral: f(A) = (1/2πi) ∮_C f(z)(zI − A)⁻¹ dz.

        Uses:
        * circular contour with radius = 3·(||A−cI||_2 + 0.5)
        * 192 Gauss-Legendre nodes (dense enough for machine eps)
        * high-precision f(z) via mpmath at every node
        """
        n = A.shape[0]
        if n == 1:
            return np.array([[self._scalar_ml(A[0, 0])]], dtype=complex)

        center = np.trace(A) / n
        eigs = np.linalg.eigvals(A - center * np.eye(n))
        radius = 3.0 * (float(np.max(np.abs(eigs))) + 0.5)

        n_quad = 192
        nodes, weights = np.polynomial.legendre.leggauss(n_quad)
        theta = np.pi * (nodes + 1.0)
        w = np.pi * weights
        z_pts = center + radius * np.exp(1j * theta)
        dz_dtheta = 1j * radius * np.exp(1j * theta)   # ← YENİ

        acc = np.zeros((n, n), dtype=complex)
        eye = np.eye(n, dtype=complex)

        for z, wz, dzi in zip(z_pts, w, dz_dtheta):   # ← dzi ekle
            try:
                fz = self._scalar_ml(complex(z))
            except MittagLefflerError:
                from mlfp.mlfp import _mpmath_fallback
                fz = _mpmath_fallback(
                    complex(z), self.alpha, self.beta, self.gamma,
                )
            if not np.isfinite(fz):
                raise MittagLefflerError(
                    f"Contour: non-finite f(z) at z={z:.3e}"
                )
            inv = np.linalg.inv(z * eye - A)
            acc += wz * fz * inv * dzi                  # ← × dzi

        return acc / (2j * np.pi)


# ============================================================================
# Fractional Schrödinger solver
# ============================================================================
def solve_fractional_schrodinger(
    H: ArrayLike,
    psi0: ArrayLike,
    t: float,
    *,
    alpha: float = 1.0,
    hbar: float = 1.0,
    method: str = "schur",
    tolerance: float = 1e-13,
) -> NDArray[np.complex128]:
    r"""
    Solve the (time-)fractional Schrödinger equation:

        (iħ)^α ∂_t^α ψ = H ψ,   ψ(0) = ψ₀

    whose solution is

        ψ(t) = E_{α,1}(−(iħ)⁻ᵅ H t^α) ψ₀.

    For α = 1, this reduces to the standard Schrödinger equation:

        ψ(t) = e^{−iHt/ħ} ψ₀.

    Parameters
    ----------
    H : array_like, shape (n, n)
        Hamiltonian (Hermitian).
    psi0 : array_like, shape (n,)
        Initial state.
    t : float
        Evolution time.
    alpha : float, default 1.0
        Fractional order.
    hbar : float, default 1.0
        Reduced Planck constant.
    method : str
        See :class:`MittagLefflerMatrix`.
    tolerance : float

    Returns
    -------
    ndarray, shape (n,)
        Final state ψ(t).

    Examples
    --------
    >>> H = np.array([[0, 0.5], [0.5, 0]])   # two-level system
    >>> psi0 = np.array([1.0, 0.0])
    >>> # Classical QM limit
    >>> psi_t = solve_fractional_schrodinger(H, psi0, t=1.0, alpha=1.0)
    >>> # E_{1,1}(z)=e^z, so this should match scipy.linalg.expm
    >>> expected = sla.expm(-1j * H * 1.0) @ psi0
    >>> np.allclose(psi_t, expected, rtol=1e-12)
    True
    """
    H_arr = np.asarray(H, dtype=np.complex128)
    if H_arr.ndim != 2 or H_arr.shape[0] != H_arr.shape[1]:
        raise MittagLefflerError("H must be a square matrix.")
    psi0_arr = np.asarray(psi0, dtype=np.complex128).ravel()
    if psi0_arr.shape[0] != H_arr.shape[0]:
        raise MittagLefflerError(
            f"psi0 length {psi0_arr.size} ≠ H size {H_arr.shape[0]}"
        )

    # Matrix argument: −(iħ)⁻ᵅ H t^α
    i_hbar_alpha = (1j * hbar) ** alpha
    A = (H_arr * t ** alpha) / i_hbar_alpha

    ml = MittagLefflerMatrix(
        alpha=alpha, beta=1.0, method=method, tolerance=tolerance,
    )
    U = ml(A)
    return U @ psi0_arr


# ============================================================================
# Export
# ============================================================================
__all__ = ["MittagLefflerMatrix", "solve_fractional_schrodinger"]
