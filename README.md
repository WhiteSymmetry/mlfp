# mlf
mlf: Mittag-Leffler and Prabhakar function, MLF, Mittag-Leffler and Prabhakar function evaluation: Mittag-Leffler ve Prabhakar fonksiyonlarının hesaplanması.

---

Parameter validation tests
=======================================================
  [OK]   negative alpha → MittagLefflerError
  [OK]   zero alpha → MittagLefflerError
  [OK]   nan alpha → MittagLefflerError
  [OK]   infinite beta → MittagLefflerError
  [OK]   zero gamma → MittagLefflerError
  [OK]   tolerance = 1 → MittagLefflerError
  [OK]   tolerance = 0 → MittagLefflerError
  [OK]   negative series_threshold → MittagLefflerError
  [OK]   zero max_series_terms → MittagLefflerError
  [OK]   complex z to __call__ → MittagLefflerError

========================================================
Safe-limit overflow detection
====================================================================================================
  [OK  ] α=0.25, |z|=10.3 → MittagLefflerError (The Mittag-Leffler value is too large to represent in f…)
  [OK  ] α=0.5 , |z|=49 → MittagLefflerError (The Mittag-Leffler value is too large to represent in f…)
  [OK  ] α=1.0 , |z|=1.28e+03 → MittagLefflerError (The Mittag-Leffler value is too large to represent in f…)
  [OK  ] α=2.0 , |z|=9.07e+05 → MittagLefflerError (The Mittag-Leffler value is too large to represent in f…)

============================================================================
Numerical comparison:  local  vs  mpmath (ground truth)  vs  pymittagleffler
  local tolerance      : 1.0e-13
  mpmath digits        : dynamic (>=50, up to 600)
  pass criterion       : relative error < 1e-09
=============================================================================

angle = 0.000000 rad  (0.0°)
----------------------------------------------------------------------------------------------------
  [OK  ] α=0.25, β=0.5  | n= 13 | vs mpmath: max_abs=1.11e+190, max_rel=7.23e-14 | vs pymtl: max_rel=7.19e-14
  [OK  ] α=0.25, β=1.0  | n= 13 | vs mpmath: max_abs=5.22e+188, max_rel=1.38e-14 | vs pymtl: max_rel=9.21e-15
  [OK  ] α=0.25, β=1.5  | n= 13 | vs mpmath: max_abs=2.39e+187, max_rel=1.38e-14 | vs pymtl: max_rel=6.88e-15
  [OK  ] α=0.25, β=2.0  | n= 13 | vs mpmath: max_abs=1.10e+186, max_rel=1.36e-14 | vs pymtl: max_rel=7.22e-15
  [OK  ] α=0.5 , β=0.5  | n= 32 | vs mpmath: max_abs=1.25e+45, max_rel=3.36e-14 | vs pymtl: max_rel=3.36e-14
  [OK  ] α=0.5 , β=1.0  | n= 32 | vs mpmath: max_abs=1.12e+44, max_rel=7.37e-15 | vs pymtl: max_rel=7.51e-15
  [OK  ] α=0.5 , β=1.5  | n= 32 | vs mpmath: max_abs=1.14e+46, max_rel=4.79e-15 | vs pymtl: max_rel=4.79e-15
  [OK  ] α=0.5 , β=2.0  | n= 32 | vs mpmath: max_abs=7.14e+44, max_rel=4.81e-15 | vs pymtl: max_rel=1.20e-15
  [OK  ] α=1.0 , β=0.5  | n= 32 | vs mpmath: max_abs=1.16e-10, max_rel=1.56e-14 | vs pymtl: max_rel=1.54e-14
  [OK  ] α=1.0 , β=1.0  | n= 32 | vs mpmath: max_abs=7.99e-15, max_rel=2.49e-15 | vs pymtl: max_rel=2.49e-15
  [OK  ] α=1.0 , β=1.5  | n= 32 | vs mpmath: max_abs=7.28e-12, max_rel=1.49e-15 | vs pymtl: max_rel=1.55e-15
  [OK  ] α=1.0 , β=2.0  | n= 32 | vs mpmath: max_abs=1.82e-12, max_rel=1.68e-14 | vs pymtl: max_rel=1.67e-14
  [OK  ] α=1.5 , β=0.5  | n= 32 | vs mpmath: max_abs=2.27e-13, max_rel=4.83e-14 | vs pymtl: max_rel=4.82e-14
  [OK  ] α=1.5 , β=1.0  | n= 32 | vs mpmath: max_abs=8.53e-14, max_rel=1.76e-14 | vs pymtl: max_rel=1.75e-14
  [OK  ] α=1.5 , β=1.5  | n= 32 | vs mpmath: max_abs=3.55e-14, max_rel=5.79e-15 | vs pymtl: max_rel=5.65e-15
  [OK  ] α=1.5 , β=2.0  | n= 32 | vs mpmath: max_abs=1.78e-14, max_rel=2.29e-15 | vs pymtl: max_rel=2.12e-15
  [OK  ] α=2.0 , β=0.5  | n= 32 | vs mpmath: max_abs=4.18e-14, max_rel=2.68e-14 | vs pymtl: max_rel=2.67e-14
  [OK  ] α=2.0 , β=1.0  | n= 32 | vs mpmath: max_abs=2.22e-14, max_rel=1.35e-14 | vs pymtl: max_rel=1.35e-14
  [OK  ] α=2.0 , β=1.5  | n= 32 | vs mpmath: max_abs=6.69e-15, max_rel=4.45e-15 | vs pymtl: max_rel=4.30e-15
  [OK  ] α=2.0 , β=2.0  | n= 32 | vs mpmath: max_abs=4.22e-15, max_rel=2.78e-15 | vs pymtl: max_rel=2.78e-15

angle = 0.785398 rad  (45.0°)
----------------------------------------------------------------------------------------------------
  [OK  ] α=0.25, β=0.5  | n= 13 | vs mpmath: max_abs=3.40e-14, max_rel=3.40e-14 | vs pymtl: max_rel=3.41e-14
  [OK  ] α=0.25, β=1.0  | n= 13 | vs mpmath: max_abs=3.42e-14, max_rel=2.57e-14 | vs pymtl: max_rel=2.62e-14
  [OK  ] α=0.25, β=1.5  | n= 13 | vs mpmath: max_abs=4.15e-13, max_rel=2.86e-13 | vs pymtl: max_rel=2.84e-13
  [OK  ] α=0.25, β=2.0  | n= 13 | vs mpmath: max_abs=2.43e-15, max_rel=1.95e-15 | vs pymtl: max_rel=1.93e-15
  [OK  ] α=0.5 , β=0.5  | n= 32 | vs mpmath: max_abs=3.10e-13, max_rel=6.63e-14 | vs pymtl: max_rel=6.64e-14
  [OK  ] α=0.5 , β=1.0  | n= 32 | vs mpmath: max_abs=2.65e-14, max_rel=1.89e-14 | vs pymtl: max_rel=1.86e-14
  [OK  ] α=0.5 , β=1.5  | n= 32 | vs mpmath: max_abs=5.12e-14, max_rel=4.27e-14 | vs pymtl: max_rel=4.22e-14
  [OK  ] α=0.5 , β=2.0  | n= 32 | vs mpmath: max_abs=4.50e-15, max_rel=3.67e-15 | vs pymtl: max_rel=4.03e-15
  [OK  ] α=1.0 , β=0.5  | n= 32 | vs mpmath: max_abs=2.57e-12, max_rel=4.69e-14 | vs pymtl: max_rel=4.71e-14
  [OK  ] α=1.0 , β=1.0  | n= 32 | vs mpmath: max_abs=4.58e-13, max_rel=1.21e-14 | vs pymtl: max_rel=1.22e-14
  [OK  ] α=1.0 , β=1.5  | n= 32 | vs mpmath: max_abs=2.34e-13, max_rel=4.29e-15 | vs pymtl: max_rel=4.56e-15
  [OK  ] α=1.0 , β=2.0  | n= 32 | vs mpmath: max_abs=5.68e-14, max_rel=2.58e-14 | vs pymtl: max_rel=2.56e-14
  [OK  ] α=1.5 , β=0.5  | n= 32 | vs mpmath: max_abs=1.62e-13, max_rel=3.19e-14 | vs pymtl: max_rel=3.19e-14
  [OK  ] α=1.5 , β=1.0  | n= 32 | vs mpmath: max_abs=7.42e-14, max_rel=9.07e-15 | vs pymtl: max_rel=9.25e-15
  [OK  ] α=1.5 , β=1.5  | n= 32 | vs mpmath: max_abs=2.84e-14, max_rel=3.46e-15 | vs pymtl: max_rel=3.42e-15
  [OK  ] α=1.5 , β=2.0  | n= 32 | vs mpmath: max_abs=1.29e-14, max_rel=2.05e-15 | vs pymtl: max_rel=1.90e-15
  [OK  ] α=2.0 , β=0.5  | n= 32 | vs mpmath: max_abs=5.78e-14, max_rel=5.78e-14 | vs pymtl: max_rel=5.95e-14
  [OK  ] α=2.0 , β=1.0  | n= 32 | vs mpmath: max_abs=1.12e-14, max_rel=9.77e-15 | vs pymtl: max_rel=9.76e-15
  [OK  ] α=2.0 , β=1.5  | n= 32 | vs mpmath: max_abs=5.57e-15, max_rel=4.59e-15 | vs pymtl: max_rel=4.77e-15
  [OK  ] α=2.0 , β=2.0  | n= 32 | vs mpmath: max_abs=2.07e-15, max_rel=1.97e-15 | vs pymtl: max_rel=1.75e-15

angle = 1.570796 rad  (90.0°)
----------------------------------------------------------------------------------------------------
  [OK  ] α=0.25, β=0.5  | n= 13 | vs mpmath: max_abs=2.56e-14, max_rel=2.56e-14 | vs pymtl: max_rel=2.57e-14
  [OK  ] α=0.25, β=1.0  | n= 13 | vs mpmath: max_abs=7.66e-15, max_rel=7.66e-15 | vs pymtl: max_rel=7.21e-15
  [OK  ] α=0.25, β=1.5  | n= 13 | vs mpmath: max_abs=1.64e-13, max_rel=1.55e-13 | vs pymtl: max_rel=1.53e-13
  [OK  ] α=0.25, β=2.0  | n= 13 | vs mpmath: max_abs=1.28e-15, max_rel=1.28e-15 | vs pymtl: max_rel=1.26e-15
  [OK  ] α=0.5 , β=0.5  | n= 32 | vs mpmath: max_abs=2.87e-14, max_rel=2.87e-14 | vs pymtl: max_rel=2.96e-14
  [OK  ] α=0.5 , β=1.0  | n= 32 | vs mpmath: max_abs=1.80e-14, max_rel=1.80e-14 | vs pymtl: max_rel=1.85e-14
  [OK  ] α=0.5 , β=1.5  | n= 32 | vs mpmath: max_abs=1.04e-14, max_rel=1.04e-14 | vs pymtl: max_rel=9.12e-15
  [OK  ] α=0.5 , β=2.0  | n= 32 | vs mpmath: max_abs=2.59e-15, max_rel=2.59e-15 | vs pymtl: max_rel=2.69e-15
  [OK  ] α=1.0 , β=0.5  | n= 32 | vs mpmath: max_abs=3.29e-14, max_rel=2.95e-14 | vs pymtl: max_rel=2.97e-14
  [OK  ] α=1.0 , β=1.0  | n= 32 | vs mpmath: max_abs=1.39e-14, max_rel=1.39e-14 | vs pymtl: max_rel=1.39e-14
  [OK  ] α=1.0 , β=1.5  | n= 32 | vs mpmath: max_abs=6.91e-15, max_rel=6.51e-15 | vs pymtl: max_rel=6.69e-15
  [OK  ] α=1.0 , β=2.0  | n= 32 | vs mpmath: max_abs=2.84e-14, max_rel=2.84e-14 | vs pymtl: max_rel=2.84e-14
  [OK  ] α=1.5 , β=0.5  | n= 32 | vs mpmath: max_abs=7.39e-14, max_rel=6.67e-14 | vs pymtl: max_rel=6.68e-14
  [OK  ] α=1.5 , β=1.0  | n= 32 | vs mpmath: max_abs=2.46e-14, max_rel=2.08e-14 | vs pymtl: max_rel=2.06e-14
  [OK  ] α=1.5 , β=1.5  | n= 32 | vs mpmath: max_abs=9.22e-15, max_rel=7.74e-15 | vs pymtl: max_rel=7.60e-15
  [OK  ] α=1.5 , β=2.0  | n= 32 | vs mpmath: max_abs=3.94e-15, max_rel=3.91e-15 | vs pymtl: max_rel=3.84e-15
  [OK  ] α=2.0 , β=0.5  | n= 32 | vs mpmath: max_abs=4.00e-14, max_rel=3.14e-14 | vs pymtl: max_rel=3.10e-14
  [OK  ] α=2.0 , β=1.0  | n= 32 | vs mpmath: max_abs=1.77e-14, max_rel=1.60e-14 | vs pymtl: max_rel=1.60e-14
  [OK  ] α=2.0 , β=1.5  | n= 32 | vs mpmath: max_abs=4.25e-15, max_rel=3.67e-15 | vs pymtl: max_rel=4.48e-15
  [OK  ] α=2.0 , β=2.0  | n= 32 | vs mpmath: max_abs=2.00e-15, max_rel=1.99e-15 | vs pymtl: max_rel=1.90e-15

angle = 2.356194 rad  (135.0°)
----------------------------------------------------------------------------------------------------
  [OK  ] α=0.25, β=0.5  | n= 13 | vs mpmath: max_abs=2.23e-14, max_rel=2.23e-14 | vs pymtl: max_rel=2.25e-14
  [OK  ] α=0.25, β=1.0  | n= 13 | vs mpmath: max_abs=6.42e-15, max_rel=6.42e-15 | vs pymtl: max_rel=5.63e-15
  [OK  ] α=0.25, β=1.5  | n= 13 | vs mpmath: max_abs=1.20e-13, max_rel=1.20e-13 | vs pymtl: max_rel=1.19e-13
  [OK  ] α=0.25, β=2.0  | n= 13 | vs mpmath: max_abs=1.69e-15, max_rel=1.69e-15 | vs pymtl: max_rel=1.39e-15
  [OK  ] α=0.5 , β=0.5  | n= 32 | vs mpmath: max_abs=2.34e-14, max_rel=2.34e-14 | vs pymtl: max_rel=2.34e-14
  [OK  ] α=0.5 , β=1.0  | n= 32 | vs mpmath: max_abs=3.74e-15, max_rel=3.74e-15 | vs pymtl: max_rel=3.93e-15
  [OK  ] α=0.5 , β=1.5  | n= 32 | vs mpmath: max_abs=4.13e-14, max_rel=4.13e-14 | vs pymtl: max_rel=3.60e-14
  [OK  ] α=0.5 , β=2.0  | n= 32 | vs mpmath: max_abs=2.22e-15, max_rel=2.22e-15 | vs pymtl: max_rel=1.96e-15
  [OK  ] α=1.0 , β=0.5  | n= 32 | vs mpmath: max_abs=5.39e-14, max_rel=5.39e-14 | vs pymtl: max_rel=5.41e-14
  [OK  ] α=1.0 , β=1.0  | n= 32 | vs mpmath: max_abs=2.61e-14, max_rel=2.61e-14 | vs pymtl: max_rel=2.61e-14
  [OK  ] α=1.0 , β=1.5  | n= 32 | vs mpmath: max_abs=8.89e-15, max_rel=8.89e-15 | vs pymtl: max_rel=8.55e-15
  [OK  ] α=1.0 , β=2.0  | n= 32 | vs mpmath: max_abs=9.14e-15, max_rel=9.14e-15 | vs pymtl: max_rel=9.14e-15
  [OK  ] α=1.5 , β=0.5  | n= 32 | vs mpmath: max_abs=5.09e-14, max_rel=5.09e-14 | vs pymtl: max_rel=4.92e-14
  [OK  ] α=1.5 , β=1.0  | n= 32 | vs mpmath: max_abs=2.30e-14, max_rel=2.30e-14 | vs pymtl: max_rel=2.27e-14
  [OK  ] α=1.5 , β=1.5  | n= 32 | vs mpmath: max_abs=8.24e-15, max_rel=8.24e-15 | vs pymtl: max_rel=8.31e-15
  [OK  ] α=1.5 , β=2.0  | n= 32 | vs mpmath: max_abs=3.18e-15, max_rel=3.18e-15 | vs pymtl: max_rel=3.03e-15
  [OK  ] α=2.0 , β=0.5  | n= 32 | vs mpmath: max_abs=3.34e-14, max_rel=3.34e-14 | vs pymtl: max_rel=3.34e-14
  [OK  ] α=2.0 , β=1.0  | n= 32 | vs mpmath: max_abs=1.90e-14, max_rel=1.90e-14 | vs pymtl: max_rel=1.91e-14
  [OK  ] α=2.0 , β=1.5  | n= 32 | vs mpmath: max_abs=4.92e-15, max_rel=4.92e-15 | vs pymtl: max_rel=4.91e-15
  [OK  ] α=2.0 , β=2.0  | n= 32 | vs mpmath: max_abs=1.89e-15, max_rel=1.89e-15 | vs pymtl: max_rel=2.00e-15

angle = 3.141593 rad  (180.0°)
----------------------------------------------------------------------------------------------------
  [OK  ] α=0.25, β=0.5  | n= 13 | vs mpmath: max_abs=2.08e-14, max_rel=2.08e-14 | vs pymtl: max_rel=2.09e-14
  [OK  ] α=0.25, β=1.0  | n= 13 | vs mpmath: max_abs=6.99e-15, max_rel=6.99e-15 | vs pymtl: max_rel=6.44e-15
  [OK  ] α=0.25, β=1.5  | n= 13 | vs mpmath: max_abs=1.11e-13, max_rel=1.11e-13 | vs pymtl: max_rel=1.10e-13
  [OK  ] α=0.25, β=2.0  | n= 13 | vs mpmath: max_abs=1.11e-15, max_rel=1.11e-15 | vs pymtl: max_rel=1.78e-15
  [OK  ] α=0.5 , β=0.5  | n= 32 | vs mpmath: max_abs=2.22e-14, max_rel=2.22e-14 | vs pymtl: max_rel=2.23e-14
  [OK  ] α=0.5 , β=1.0  | n= 32 | vs mpmath: max_abs=3.70e-15, max_rel=3.70e-15 | vs pymtl: max_rel=3.81e-15
  [OK  ] α=0.5 , β=1.5  | n= 32 | vs mpmath: max_abs=4.13e-14, max_rel=4.13e-14 | vs pymtl: max_rel=3.61e-14
  [OK  ] α=0.5 , β=2.0  | n= 32 | vs mpmath: max_abs=2.55e-15, max_rel=2.55e-15 | vs pymtl: max_rel=2.33e-15
  [OK  ] α=1.0 , β=0.5  | n= 32 | vs mpmath: max_abs=2.61e-14, max_rel=2.61e-14 | vs pymtl: max_rel=2.54e-14
  [OK  ] α=1.0 , β=1.0  | n= 32 | vs mpmath: max_abs=6.63e-15, max_rel=6.63e-15 | vs pymtl: max_rel=6.63e-15
  [OK  ] α=1.0 , β=1.5  | n= 32 | vs mpmath: max_abs=8.88e-15, max_rel=8.88e-15 | vs pymtl: max_rel=1.01e-14
  [OK  ] α=1.0 , β=2.0  | n= 32 | vs mpmath: max_abs=6.48e-14, max_rel=6.48e-14 | vs pymtl: max_rel=6.48e-14
  [OK  ] α=1.5 , β=0.5  | n= 32 | vs mpmath: max_abs=6.80e-14, max_rel=6.80e-14 | vs pymtl: max_rel=6.79e-14
  [OK  ] α=1.5 , β=1.0  | n= 32 | vs mpmath: max_abs=1.55e-14, max_rel=1.55e-14 | vs pymtl: max_rel=1.60e-14
  [OK  ] α=1.5 , β=1.5  | n= 32 | vs mpmath: max_abs=5.44e-15, max_rel=5.44e-15 | vs pymtl: max_rel=5.45e-15
  [OK  ] α=1.5 , β=2.0  | n= 32 | vs mpmath: max_abs=3.27e-15, max_rel=3.27e-15 | vs pymtl: max_rel=3.37e-15
  [OK  ] α=2.0 , β=0.5  | n= 32 | vs mpmath: max_abs=5.46e-14, max_rel=5.46e-14 | vs pymtl: max_rel=5.46e-14
  [OK  ] α=2.0 , β=1.0  | n= 32 | vs mpmath: max_abs=2.10e-14, max_rel=2.10e-14 | vs pymtl: max_rel=2.11e-14
  [OK  ] α=2.0 , β=1.5  | n= 32 | vs mpmath: max_abs=5.79e-15, max_rel=5.79e-15 | vs pymtl: max_rel=5.91e-15
  [OK  ] α=2.0 , β=2.0  | n= 32 | vs mpmath: max_abs=2.31e-15, max_rel=2.31e-15 | vs pymtl: max_rel=2.42e-15

angle = 3.926991 rad  (225.0°)
----------------------------------------------------------------------------------------------------
  [OK  ] α=0.25, β=0.5  | n= 13 | vs mpmath: max_abs=2.11e-14, max_rel=2.11e-14 | vs pymtl: max_rel=2.13e-14
  [OK  ] α=0.25, β=1.0  | n= 13 | vs mpmath: max_abs=5.88e-15, max_rel=5.88e-15 | vs pymtl: max_rel=5.08e-15
  [OK  ] α=0.25, β=1.5  | n= 13 | vs mpmath: max_abs=1.21e-13, max_rel=1.21e-13 | vs pymtl: max_rel=1.19e-13
  [OK  ] α=0.25, β=2.0  | n= 13 | vs mpmath: max_abs=1.08e-15, max_rel=1.08e-15 | vs pymtl: max_rel=1.12e-15
  [OK  ] α=0.5 , β=0.5  | n= 32 | vs mpmath: max_abs=2.50e-14, max_rel=2.50e-14 | vs pymtl: max_rel=2.53e-14
  [OK  ] α=0.5 , β=1.0  | n= 32 | vs mpmath: max_abs=4.33e-15, max_rel=4.33e-15 | vs pymtl: max_rel=4.34e-15
  [OK  ] α=0.5 , β=1.5  | n= 32 | vs mpmath: max_abs=4.14e-14, max_rel=4.14e-14 | vs pymtl: max_rel=3.64e-14
  [OK  ] α=0.5 , β=2.0  | n= 32 | vs mpmath: max_abs=2.52e-15, max_rel=2.52e-15 | vs pymtl: max_rel=2.37e-15
  [OK  ] α=1.0 , β=0.5  | n= 32 | vs mpmath: max_abs=5.00e-14, max_rel=5.00e-14 | vs pymtl: max_rel=5.01e-14
  [OK  ] α=1.0 , β=1.0  | n= 32 | vs mpmath: max_abs=2.64e-14, max_rel=2.64e-14 | vs pymtl: max_rel=2.64e-14
  [OK  ] α=1.0 , β=1.5  | n= 32 | vs mpmath: max_abs=7.90e-15, max_rel=7.90e-15 | vs pymtl: max_rel=7.57e-15
  [OK  ] α=1.0 , β=2.0  | n= 32 | vs mpmath: max_abs=9.13e-15, max_rel=9.13e-15 | vs pymtl: max_rel=9.13e-15
  [OK  ] α=1.5 , β=0.5  | n= 32 | vs mpmath: max_abs=6.39e-14, max_rel=6.39e-14 | vs pymtl: max_rel=6.39e-14
  [OK  ] α=1.5 , β=1.0  | n= 32 | vs mpmath: max_abs=1.63e-14, max_rel=1.63e-14 | vs pymtl: max_rel=1.64e-14
  [OK  ] α=1.5 , β=1.5  | n= 32 | vs mpmath: max_abs=6.86e-15, max_rel=6.86e-15 | vs pymtl: max_rel=6.96e-15
  [OK  ] α=1.5 , β=2.0  | n= 32 | vs mpmath: max_abs=3.22e-15, max_rel=3.22e-15 | vs pymtl: max_rel=3.33e-15
  [OK  ] α=2.0 , β=0.5  | n= 32 | vs mpmath: max_abs=3.13e-14, max_rel=3.13e-14 | vs pymtl: max_rel=3.14e-14
  [OK  ] α=2.0 , β=1.0  | n= 32 | vs mpmath: max_abs=1.91e-14, max_rel=1.91e-14 | vs pymtl: max_rel=1.93e-14
  [OK  ] α=2.0 , β=1.5  | n= 32 | vs mpmath: max_abs=4.45e-15, max_rel=4.45e-15 | vs pymtl: max_rel=4.56e-15
  [OK  ] α=2.0 , β=2.0  | n= 32 | vs mpmath: max_abs=2.01e-15, max_rel=2.01e-15 | vs pymtl: max_rel=2.01e-15

=====================================================================
Series vs OPC consistency (gamma = 1) — small |z|, both methods valid
=====================================================================
  [OK  ] α=0.5 , β=1.0  | max_abs=2.01e-13, max_rel=5.13e-14
  [OK  ] α=1.0 , β=1.5  | max_abs=3.65e-14, max_rel=1.18e-14
  [OK  ] α=1.5 , β=0.5  | max_abs=8.20e-14, max_rel=5.60e-14
  [OK  ] α=2.0 , β=2.0  | max_abs=3.11e-15, max_rel=2.66e-15

=====================================================================
Sample values (angle = π/2, mpmath reference uses adaptive precision)
=====================================================================
  alpha   beta          r         local Re           ref Re         local Im           ref Im
---------------------------------------------------------------------------------------------
  0.500  1.000     0.0000   1.00000000e+00   1.00000000e+00   0.00000000e+00   0.00000000e+00
  0.500  1.000     1.5484   9.09447953e-02   9.09447953e-02   4.67700971e-01   4.67700971e-01
  0.500  1.000     3.0968   6.84087098e-05   6.84087098e-05   1.93863251e-01   1.93863251e-01
  0.500  1.000     5.0323   1.00501737e-11   1.00482567e-11   1.14474501e-01   1.14474501e-01
  0.500  1.000     6.5806   1.77442800e-15  -5.28709605e-18   8.67610345e-02   8.67610345e-02
  0.500  1.000     8.5161   1.79246827e-15  -4.14355680e-18   6.67160627e-02   6.67160627e-02
  0.500  1.000    10.0645   1.04123319e-15  -3.48464925e-18   5.63382053e-02   5.63382053e-02
  0.500  1.000    12.0000   1.13372614e-15  -2.90940947e-18   4.71807787e-02   4.71807787e-02
---------------------------------------------------------------------------------------------
  1.000  1.500     0.0000   1.12837917e+00   1.12837917e+00   0.00000000e+00   0.00000000e+00
  1.000  1.500     1.5484   5.09688802e-01   5.09688802e-01   8.75068324e-01   8.75068324e-01
  1.000  1.500     3.0968  -4.07036496e-01  -4.07036496e-01   5.92423924e-01   5.92423924e-01
  1.000  1.500     5.0323  -2.10109182e-01  -2.10109182e-01  -2.88870382e-01  -2.88870382e-01
  1.000  1.500     6.5806   3.38254207e-01   3.38254207e-01  -9.82928006e-02  -9.82928006e-02
  1.000  1.500     8.5161   3.84081301e-02   3.84081301e-02   4.05701628e-01   4.05701628e-01
  1.000  1.500    10.0645  -3.14568714e-01  -3.14568714e-01   1.01423543e-01   1.01423543e-01
  1.000  1.500    12.0000   6.08109247e-02   6.08109247e-02  -2.34994626e-01  -2.34994626e-01
---------------------------------------------------------------------------------------------
  2.000  2.000     0.0000   1.00000000e+00   1.00000000e+00   0.00000000e+00   0.00000000e+00
  2.000  2.000     1.5484   9.80036649e-01   9.80036649e-01   2.57328179e-01   2.57328179e-01
  2.000  2.000     3.0968   9.20336545e-01   9.20336545e-01   5.10243686e-01   5.10243686e-01
  2.000  2.000     5.0323   7.90734428e-01   7.90734428e-01   8.13505739e-01   8.13505739e-01
  2.000  2.000     6.5806   6.44280740e-01   6.44280740e-01   1.04054046e+00   1.04054046e+00
  2.000  2.000     8.5161   4.10062960e-01   4.10062960e-01   1.29792936e+00   1.29792936e+00
  2.000  2.000    10.0645   1.83988068e-01   1.83988068e-01   1.47772066e+00   1.47772066e+00
  2.000  2.000    12.0000  -1.43335456e-01  -1.43335456e-01   1.66334926e+00   1.66334926e+00
---------------------------------------------------------------------------------------------

====================
Overall: PASS ✅
====================

---


