"""
mlfp: Mittag-Leffler and Prabhakar function, MLF, 
Mittag-Leffler and Prabhakar function evaluation: 
Mittag-Leffler ve Prabhakar fonksiyonlarının hesaplanması.

- Two-parameter Mittag-Leffler function, gamma = 1:

      E_{alpha,beta}(z)

  is evaluated with Garrappa's optimal parabolic contour (OPC)
  Laplace-transform inversion, with a Taylor series branch for
  small |z|.

- Three-parameter Mittag-Leffler / Prabhakar function, gamma != 1:

      E_{alpha,beta}^gamma(z)

  is evaluated from its defining power series. The OPC residue
  treatment is not mathematically valid for gamma != 1 because the
  singularities are generally branch points rather than simple poles.

Compatible with Python 3.11-3.15.
"""

__version__ = "0.1.8"
__author__ = "Mehmet Keçeci"
__main__ = "0.1.8"

from .mlfp import (
    # Temel eğri sınıfları
    MittagLefflerError,
    _safe_z_limit,
    _find_optimal_bounded_param,
    _find_optimal_unbounded_param,
    _invert_laplace_transform_two_parameter,
    MittagLeffler,
    mittag_leffler,
    _asymptotic_expansion,
    _mpmath_fallback,
    _reduced_evaluator,
    derivative,
    _normalize_param,
    FractionalQuantumChannel,
    MatrixMittagLeffler,
    Prabhakar,
    _is_diagonal,
    _is_hermitian,
    CompositeFractionalChannel,
    FractionalAmplitudeDamping,
    two_qubit_amplitude_damping,
    local_amplitude_damping,
    amplitude_damping_channel,
    amplitude_damping_kraus,
    MatrixAnalyzer,
    _is_sparse,
    analyze_and_plot_transform,
    _eigenvalue_comparison,
    _matrix_heatmap_axis,
    _summary_heatmap,
    _print_report,
    _collect_property_names,
    DASHBOARD_PROPERTY_ORDER,
    plot_time_process_eigenvalues,
    plot_time_process_metrics,
    time_process,
    _evaluate_with_kwargs,
    SparseHermitianMatrixMittagLeffler,
    NV_hamiltonian,
    build_surface_code_hamiltonian,
    topological_invariant,
    FractionalLindbladianProcess,
    fast_diagonal_process,
    vec,
    unvec,
    build_lindbladian,
    kitaev_chain_hamiltonian,
    kitaev_winding,


)

from mlfp.mlfp_matrix import (
    MittagLefflerMatrix,
    solve_fractional_schrodinger,
)

from mlfp.mlfp_bdg import (
    build_bdg_hamiltonian,
    kitaev_chain,
    solve_fractional_bdg,
    bdg_spectrum,
    find_majorana_modes,
    majorana_localization,
)

__all__ = [
    "MittagLefflerError",
    "_safe_z_limit",
    "_find_optimal_bounded_param",
    "_find_optimal_unbounded_param",
    "_invert_laplace_transform_two_parameter",
    "MittagLeffler",
    "mittag_leffler",
    "_asymptotic_expansion",
    "_mpmath_fallback",
    "_reduced_evaluator",
    "derivative",
    "_normalize_param",
    "MittagLefflerMatrix",
    "solve_fractional_schrodinger",
    "build_bdg_hamiltonian",
    "kitaev_chain",
    "solve_fractional_bdg",
    "bdg_spectrum",
    "find_majorana_modes",
    "majorana_localization",
    "FractionalQuantumChannel",
    "MatrixMittagLeffler",
    "Prabhakar",
    "_is_diagonal",
    "_is_hermitian",
    "CompositeFractionalChannel",
    "FractionalAmplitudeDamping",
    "two_qubit_amplitude_damping",
    "local_amplitude_damping",
    "amplitude_damping_channel",
    "amplitude_damping_kraus",
    "MatrixAnalyzer",
    "_is_sparse",
    "analyze_and_plot_transform",
    "_eigenvalue_comparison",
    "_matrix_heatmap_axis",
    "_summary_heatmap",
    "_print_report",
    "_collect_property_names",
    "DASHBOARD_PROPERTY_ORDER",
    "plot_time_process_eigenvalues",
    "plot_time_process_metrics",
    "time_process",
    "_evaluate_with_kwargs",
    "SparseHermitianMatrixMittagLeffler",
    "NV_hamiltonian",
    "kitaev_chain",
    "build_surface_code_hamiltonian",
    "topological_invariant",
    "FractionalLindbladianProcess",
    "fast_diagonal_process",
    "vec",
    "unvec",
    "build_lindbladian",
    "kitaev_chain_hamiltonian",
    "kitaev_winding",
]
