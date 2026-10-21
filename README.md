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

Numerical comparison:  mlf  vs  mpmath (ground truth)  vs  pymittagleffler

  mlf tolerance      : 1.0e-13
  
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

  alpha   beta          r         mlf Re           ref Re         mlf Im           ref Im
  
---------------------------------------------------------------------------------------------

  0.500  1.000     0.0000   1.00000000e+00   1.00000000e+00   0.00000000e+00   0.00000000e+00
 
================

Overall: PASS ✅

================

---
# Examples

```python
import numpy as np
from mlf import mittag_leffler

alpha, beta = 1.0, 2.0
z = np.linspace(0.0, 2.0, 64)
result = mittag_leffler(z, alpha, beta)
print(result)

```

```python
import numpy as np
from mlf import MittagLeffler, mittag_leffler

# 1. Two-parameter Mittag-Leffler function
ml = MittagLeffler(alpha=0.5, beta=1.0)

print(ml(0.0))       # 1.0
print(ml(0.5))       # 1.9523604891825568
print(ml(1.0))       # 5.008980080762283
print(ml(2.0))       # 108.94090438997797

# 2. Array evaluation
z = np.array([0.0, 0.25, 0.5, 1.0, 2.0, 5.0])
print(ml(z))

# 3. Small arguments now use the reliable Taylor-series branch
ml07 = MittagLeffler(alpha=0.7, beta=1.0)
print(ml07(0.001))   # 1.0011013529...
print(ml07(0.01))    # 1.0110864355...
print(ml07(1.0))     # 3.7041461454...

# 4. Prabhakar / three-parameter function
prabhakar = MittagLeffler(alpha=0.5, beta=1.0, gamma=0.7)
print(prabhakar(0.1))  # 1.0853656185...
print(prabhakar(0.5))  # 1.6154153961...
print(prabhakar(1.0))  # 3.3153696560...
print(prabhakar(2.0))  # 47.137926541...

# 5. Functional interface
print(mittag_leffler(1.0, alpha=0.5, beta=1.0))
print(mittag_leffler(0.5, alpha=0.5, beta=1.0, gamma=0.7))

```

```python
import numpy as np
import mpmath as mp
from pymittagleffler import mittag_leffler as pymittag_leffler
from mlf import MittagLeffler

mp.mp.dps = 40

def mpmath_reference(z, alpha, beta, terms=2000):
    z = mp.mpc(z)
    total = mp.mpc(0)
    for k in range(terms):
        term = z**k / mp.gamma(alpha * k + beta)
        total += term
        if abs(term) < mp.mpf(10) ** -30:
            break
    return complex(total)

alpha, beta = 1.0, 1.5
r = np.linspace(0.0, 12.0, 256)
z = r * np.exp(1j * np.pi / 2)

mlf = MittagLeffler(
    alpha=alpha,
    beta=beta,
    gamma=1.0,
    tolerance=1.0e-13,
    series_threshold=0.1,
).evaluate_complex(z)

reference = pymittag_leffler(z, alpha, beta)
mp_reference = np.array(
    [mpmath_reference(value, alpha, beta) for value in z],
    dtype=np.complex128,
)

difference = np.abs(mlf - reference)
jmax = int(np.argmax(difference))

print(f"Maximum difference at r = {r[jmax]:.6f}")
print(f"z                 = {z[jmax]}")
print(f"mlf             = {mlf[jmax]}")
print(f"pymittagleffler   = {reference[jmax]}")
print(f"mpmath reference  = {mp_reference[jmax]}")
print(f"|mlf - pymtl|   = {difference[jmax]:.6e}")
print(f"|mlf - mpmath|  = {abs(mlf[jmax] - mp_reference[jmax]):.6e}")
print(f"|pymtl  - mpmath| = {abs(reference[jmax] - mp_reference[jmax]):.6e}")

```

```python
from mlf import mittag_leffler
z = np.linspace(0.4, 2.0) + 3j * np.linspace(0.5, 1.5)
ml = mittag_leffler(z, alpha=0.5, beta=2.5)
print(z)
print(ml)
```

---




