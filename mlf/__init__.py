"""
mlf: Mittag-Leffler and Prabhakar function, MLF, 
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

Compatible with Python 3.15+.
"""

__version__ = "0.1.0"
__author__ = "Mehmet Keçeci"
__main__ = "0.1.0"

from .mlf import (
    # Temel eğri sınıfları
    MittagLefflerError,
    _safe_z_limit,
    _find_optimal_bounded_param,
    _find_optimal_unbounded_param,
    _invert_laplace_transform_two_parameter,
    MittagLeffler,
    mittag_leffler,

)

__all__ = [
    "MittagLefflerError",
    "_safe_z_limit",
    "_find_optimal_bounded_param",
    "_find_optimal_unbounded_param",
    "_invert_laplace_transform_two_parameter",
    "MittagLeffler",
    vmittag_leffler",

]
