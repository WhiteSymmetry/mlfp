# SPDX-License-Identifier: MIT
r"""
Fractional Bogoliubov-de Gennes (BdG) dynamics.

Standard BdG Hamiltonian (Nambu basis):

    H_BdG = [[H_0,      Δ   ],
             [Δ†,     -H_0^T]]

Fractional evolution:

    (iħ)^α ∂_t^α Ψ = H_BdG Ψ
    Ψ(t) = E_{α,1}(H_BdG · t^α / (iħ)^α) Ψ_0

For α=1: Ψ(t) = exp(-iH_BdG t/ħ) Ψ_0   (classical BdG)

Includes Kitaev-chain helper for quick 1D topological-superconductor
experiments.
"""

from __future__ import annotations

import numpy as np
import scipy.linalg as sla
from numpy.typing import ArrayLike, NDArray

from mlfp.mlfp import MittagLefflerError
from mlfp.mlfp_matrix import MittagLefflerMatrix


__all__ = [
    "build_bdg_hamiltonian",
    "kitaev_chain",
    "solve_fractional_bdg",
    "bdg_spectrum",
    "find_majorana_modes",
    "majorana_localization",
]


# ============================================================================
# BdG Hamiltonian construction
# ============================================================================
def build_bdg_hamiltonian(
    H_0: ArrayLike,
    Delta: ArrayLike | float,
) -> NDArray[np.complex128]:
    """
    Construct the 2N×2N Bogoliubov-de Gennes Hamiltonian.

    Parameters
    ----------
    H_0 : ndarray, shape (N, N), Hermitian
        Normal-state Hamiltonian.
    Delta : ndarray (N,N) or scalar
        Superconducting pairing matrix. Scalar → Δ·I (onsite s-wave).

    Returns
    -------
    H_BdG : ndarray, shape (2N, 2N), complex128

    Notes
    -----
    Nambu spinor ordering: Ψ = (c_↑, c_↓, c_↓†, -c_↑†)^T.

    For p-wave pairing (Kitaev chain) pass an **antisymmetric** Delta:
    Δ[i,i+1] = +Δ₀, Δ[i+1,i] = -Δ₀.
    """
    H_0_arr = np.asarray(H_0, dtype=np.complex128)
    if H_0_arr.ndim != 2 or H_0_arr.shape[0] != H_0_arr.shape[1]:
        raise MittagLefflerError("H_0 must be a square matrix.")
    N = H_0_arr.shape[0]

    if np.ndim(Delta) == 0:
        Delta_arr = complex(Delta) * np.eye(N, dtype=np.complex128)
    else:
        Delta_arr = np.asarray(Delta, dtype=np.complex128)
        if Delta_arr.shape != (N, N):
            raise MittagLefflerError(
                f"Delta must be scalar or ({N},{N}); got {Delta_arr.shape}"
            )

    H_BdG = np.zeros((2 * N, 2 * N), dtype=np.complex128)
    H_BdG[:N, :N] = H_0_arr
    H_BdG[N:, N:] = -H_0_arr.T.conj()
    H_BdG[:N, N:] = Delta_arr
    H_BdG[N:, :N] = Delta_arr.T.conj()
    return H_BdG


# ============================================================================
# Kitaev-chain helper
# ============================================================================
def kitaev_chain(
    N: int,
    *,
    t: float = 1.0,
    mu: float = 0.0,
    Delta_0: float = 0.5,
) -> tuple[NDArray[np.complex128], NDArray[np.complex128]]:
    r"""
    1D Kitaev chain: N sites, hopping t, chemical potential μ,
    p-wave pairing Δ₀.

    H_0[i,i]   = -μ
    H_0[i,i±1] = -t
    Δ[i,i+1]   = +Δ₀,   Δ[i+1,i] = -Δ₀  (antisymmetric p-wave)

    Topological phase when |μ| < 2|t| → 2 Majorana zero modes at ends.

    Returns
    -------
    H_0 : ndarray (N, N), complex
    Delta : ndarray (N, N), complex
    """
    H_0 = np.zeros((N, N), dtype=np.complex128)
    for i in range(N):
        H_0[i, i] = -mu
    for i in range(N - 1):
        H_0[i, i + 1] = -t
        H_0[i + 1, i] = -t

    Delta = np.zeros((N, N), dtype=np.complex128)
    for i in range(N - 1):
        Delta[i, i + 1] = +Delta_0
        Delta[i + 1, i] = -Delta_0
    return H_0, Delta


# ============================================================================
# Fractional BdG evolution
# ============================================================================
def solve_fractional_bdg(
    H_0: ArrayLike,
    Delta: ArrayLike | float,
    psi0: ArrayLike,
    t: float,
    *,
    alpha: float = 1.0,
    hbar: float = 1.0,
    method: str = "schur",
    tolerance: float = 1e-13,
) -> NDArray[np.complex128]:
    r"""
    Solve (iħ)^α ∂_t^α Ψ = H_BdG Ψ, Ψ(0) = Ψ₀.

    Ψ(t) = E_{α,1}(H_BdG · t^α / (iħ)^α) Ψ₀

    For α=1 → classical BdG.

    Parameters
    ----------
    H_0, Delta : BdG inputs
    psi0 : ndarray, shape (2N,)
    t : float
    alpha : float, default 1.0
    hbar : float, default 1.0
    method : {"auto", "eigen", "schur", "contour"}
    tolerance : float

    Returns
    -------
    Psi_t : ndarray, shape (2N,), complex128

    Notes
    -----
    For α≠1 the evolution is NOT unitary (Naber 2004); the norm
    ``||Ψ(t)||`` is bounded but not preserved.
    """
    H_BdG = build_bdg_hamiltonian(H_0, Delta)
    psi0_arr = np.asarray(psi0, dtype=np.complex128).ravel()

    if psi0_arr.shape[0] != H_BdG.shape[0]:
        raise MittagLefflerError(
            f"psi0 length {psi0_arr.size} ≠ 2N = {H_BdG.shape[0]}"
        )

    # A = H_BdG · t^α / (iħ)^α     (α=1: A = −iH_BdG t/ħ)
    i_hbar_alpha = (1j * hbar) ** alpha
    A = H_BdG * (t ** alpha) / i_hbar_alpha

    ml = MittagLefflerMatrix(
        alpha=alpha, beta=1.0, method=method, tolerance=tolerance,
    )
    return ml(A) @ psi0_arr


# ============================================================================
# BdG spectrum
# ============================================================================
def bdg_spectrum(
    H_0: ArrayLike,
    Delta: ArrayLike | float,
    *,
    tol: float = 1e-8,
) -> dict:
    """
    Compute the BdG eigenvalue spectrum and classify zero modes.
    """
    H_BdG = build_bdg_hamiltonian(H_0, Delta)
    eigs = np.sort_complex(np.linalg.eigvals(H_BdG))
    abs_eigs = np.abs(eigs)

    zero_mask = abs_eigs < tol
    n_zero = int(zero_mask.sum())
    gap = float(abs_eigs.min()) if abs_eigs.size else 0.0

    return {
        "eigenvalues": eigs,
        "gap": gap,
        "n_zero_modes": n_zero,
        "zero_modes": eigs[zero_mask],
        "has_gap": gap > tol,
    }


# ============================================================================
# Majorana detection
# ============================================================================
def find_majorana_modes(
    H_0: ArrayLike,
    Delta: ArrayLike | float,
    *,
    tol: float = 1e-8,
) -> dict:
    """
    Identify Majorana zero modes: |E| < tol and particle-hole
    self-conjugacy.

    Returns
    -------
    dict with keys:
      n_majorana      : int (each Majorana = 1 pair in Nambu, so n_zero/2)
      energies        : ndarray of zero-eigenvalues
      wavefunctions   : ndarray (n_zero, 2N) row-wise eigenvectors
      localizations   : ndarray of inverse participation ratios
      is_topological  : bool
    """
    H_BdG = build_bdg_hamiltonian(H_0, Delta)
    N = H_BdG.shape[0] // 2

    eigs, vecs = np.linalg.eig(H_BdG)
    idx = np.argsort(np.abs(eigs))
    eigs = eigs[idx]
    vecs = vecs[:, idx]

    zero_mask = np.abs(eigs) < tol
    n_zero = int(zero_mask.sum())

    if n_zero == 0:
        return {
            "n_majorana": 0,
            "energies": np.array([], dtype=complex),
            "wavefunctions": np.zeros((0, 2 * N), dtype=complex),
            "localizations": np.array([]),
            "is_topological": False,
        }

    zero_eigs = eigs[zero_mask]
    zero_vecs = vecs[:, zero_mask].T

    ipr = []
    for psi in zero_vecs:
        rho = np.abs(psi) ** 2
        rho = rho / rho.sum()
        ipr.append(float((rho ** 2).sum()))
    ipr = np.array(ipr)

    next_eig = float(np.abs(eigs[n_zero])) if n_zero < len(eigs) else np.inf

    return {
        "n_majorana": n_zero // 2,
        "energies": zero_eigs,
        "wavefunctions": zero_vecs,
        "localizations": ipr,
        "is_topological": (n_zero >= 2) and (next_eig > 10 * tol),
    }


# ============================================================================
# Spatial localization
# ============================================================================
def majorana_localization(
    H_0: ArrayLike,
    Delta: ArrayLike | float,
    *,
    tol: float = 1e-8,
) -> dict:
    """
    Site-resolved density profile of Majorana modes:
        ρ_i = |u_i|² + |v_i|²,   i = 1..N
    """
    H_BdG = build_bdg_hamiltonian(H_0, Delta)
    N = H_BdG.shape[0] // 2

    modes = find_majorana_modes(H_0, Delta, tol=tol)
    if modes["n_majorana"] == 0:
        return {
            "profiles": np.zeros((0, N)),
            "left_edges": [],
            "right_edges": [],
            "overlap": 0.0,
        }

    profiles = []
    left_edges = []
    right_edges = []
    for k, psi in enumerate(modes["wavefunctions"]):
        u = psi[:N]
        v = psi[N:]
        rho = np.abs(u) ** 2 + np.abs(v) ** 2
        rho = rho / rho.sum()
        profiles.append(rho)
        (left_edges if rho[0] > rho[-1] else right_edges).append(k)

    profiles = np.array(profiles)

    if left_edges and right_edges:
        overlap = float(
            np.sqrt(profiles[left_edges[0]])
            @ np.sqrt(profiles[right_edges[0]])
        )
    else:
        overlap = 0.0

    return {
        "profiles": profiles,
        "left_edges": left_edges,
        "right_edges": right_edges,
        "overlap": overlap,
    }
