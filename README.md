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
  
angle = 0.785398 rad  (45.0°)

----------------------------------------------------------------------------------------------------

  [OK  ] α=0.25, β=0.5  | n= 13 | vs mpmath: max_abs=3.40e-14, max_rel=3.40e-14 | vs pymtl: max_rel=3.41e-14
  

angle = 1.570796 rad  (90.0°)

----------------------------------------------------------------------------------------------------

  [OK  ] α=0.25, β=0.5  | n= 13 | vs mpmath: max_abs=2.56e-14, max_rel=2.56e-14 | vs pymtl: max_rel=2.57e-14

angle = 2.356194 rad  (135.0°)

----------------------------------------------------------------------------------------------------

  [OK  ] α=0.25, β=0.5  | n= 13 | vs mpmath: max_abs=2.23e-14, max_rel=2.23e-14 | vs pymtl: max_rel=2.25e-14

angle = 3.141593 rad  (180.0°)

----------------------------------------------------------------------------------------------------

  [OK  ] α=0.25, β=0.5  | n= 13 | vs mpmath: max_abs=2.08e-14, max_rel=2.08e-14 | vs pymtl: max_rel=2.09e-14
 
angle = 3.926991 rad  (225.0°)

----------------------------------------------------------------------------------------------------

  [OK  ] α=0.25, β=0.5  | n= 13 | vs mpmath: max_abs=2.11e-14, max_rel=2.11e-14 | vs pymtl: max_rel=2.13e-14
  
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
 
================

Overall: PASS ✅

================

---


