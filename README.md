# mlfp
mlfp: Mittag-Leffler and Prabhakar function, MLF, Mittag-Leffler and Prabhakar function evaluation: Mittag-Leffler ve Prabhakar fonksiyonlarının hesaplanması.

---

# mlfp <img src="https://github.com/WhiteSymmetry/mlfp/blob/main/docs/logo.png" alt="mlfp" align="right" height="140"/>

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

Numerical comparison:  mlfp  vs  mpmath (ground truth)  vs  pymittagleffler

  mlfp tolerance      : 1.0e-13
  
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

  alpha   beta          r         mlfp Re           ref Re         mlfp Im           ref Im
  
---------------------------------------------------------------------------------------------

  0.500  1.000     0.0000   1.00000000e+00   1.00000000e+00   0.00000000e+00   0.00000000e+00
 
================

Overall: PASS ✅

================

---
# Examples

```python
import numpy as np
from mlfp import mittag_leffler

alpha, beta = 1.0, 2.0
z = np.linspace(0.0, 2.0, 64)
result = mittag_leffler(z, alpha, beta)
print(result)

```

```python
import numpy as np
from mlfp import MittagLeffler, mittag_leffler

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
from mlfp import MittagLeffler

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

mlfp = MittagLeffler(
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

difference = np.abs(mlfp - reference)
jmax = int(np.argmax(difference))

print(f"Maximum difference at r = {r[jmax]:.6f}")
print(f"z                 = {z[jmax]}")
print(f"mlfp             = {mlfp[jmax]}")
print(f"pymittagleffler   = {reference[jmax]}")
print(f"mpmath reference  = {mp_reference[jmax]}")
print(f"|mlfp - pymtl|   = {difference[jmax]:.6e}")
print(f"|mlfp - mpmath|  = {abs(mlfp[jmax] - mp_reference[jmax]):.6e}")
print(f"|pymtl  - mpmath| = {abs(reference[jmax] - mp_reference[jmax]):.6e}")

```

```python
import numpy as np
from mlfp import mittag_leffler
z = np.linspace(0.4, 2.0) + 3j * np.linspace(0.5, 1.5)
ml = mittag_leffler(z, alpha=0.5, beta=2.5)
print(z)
print(ml)
```

---

# mlfp — Mittag-Leffler ve Prabhakar Fonksiyonları için Python Modülü

> **PyPI:** [https://pypi.org/project/mlfp/](https://pypi.org/project/mlfp/)
> **Lisans:** MIT · **Python:** 3.11+ / 3.15 · **Bağımlılıklar:** `numpy`, `scipy`

`mlfp` (*Mittag-Leffler Functions for Python*), iki parametreli **Mittag-Leffler** fonksiyonu $E_{\alpha,\beta}(z)$ ile üç parametreli **Prabhakar** fonksiyonu $E_{\alpha,\beta}^{\gamma}(z)$'yi float64 duyarlığında güvenli ve doğru biçimde hesaplayan açık kaynaklı bir kütüphanedir.

---

## İçindekiler

1. [Kurulum](#kurulum)
2. [Hızlı Başlangıç](#hızlı-başlangıç)
3. [Matematiksel Özellikler](#matematiksel-özellikler)
4. [Sayısal Güvenlik Özellikleri](#sayısal-güvenlik-özellikleri)
5. [API Özellikleri](#api-özellikleri)
6. [Yazılım Mühendisliği Özellikleri](#yazılım-mühendisliği-özellikleri)
7. [Doğrulama](#doğrulama)
8. [Mevcut Kütüphanelerle Karşılaştırma](#mevcut-kütüphanelerle-karşılaştırma)
9. [Sınırlamalar](#sınırlamalar)
10. [Kaynaklar](#kaynaklar)

---

## Kurulum

```bash
pip install mlfp
```

Geliştirme sürümü için:

```bash
pip install "mlfp[dev]"
```

---

## Hızlı Başlangıç

```python
from mlfp import MittagLeffler, mittag_leffler
import numpy as np

# ─── Temel kullanım ───
ml = MittagLeffler(alpha=0.5, beta=1.0)
print(ml(1.0))                  # E_{0.5, 1}(1.0) → 5.008980080762283

# ─── Prabhakar fonksiyonu (γ ≠ 1) ───
prabhakar = MittagLeffler(alpha=0.5, beta=1.0, gamma=0.7)
print(prabhakar(0.5))           # E_{0.5, 1}^{0.7}(0.5) → 1.6154153961...

# ─── Vektörleştirilmiş hesap ───
z = np.array([0.1, 0.5, 1.0], dtype=complex)
print(ml.evaluate_complex(z))   # NDArray[np.complex128]

# ─── Fonksiyon tabanlı API ───
value = mittag_leffler(1.0, alpha=0.5, beta=1.0)
print(value) 
```

---

## Matematiksel Özellikler

### Desteklenen Fonksiyonlar

**İki parametreli Mittag-Leffler:**

$$
E_{\alpha,\beta}(z) = \sum_{k=0}^{\infty} \frac{z^{k}}{\Gamma(\alpha k + \beta)},
\qquad \alpha > 0,\ \beta \in \mathbb{R}
$$

**Üç parametreli Prabhakar:**

$$
E_{\alpha,\beta}^{\gamma}(z) = \frac{1}{\Gamma(\gamma)}
  \sum_{k=0}^{\infty}
  \frac{(\gamma)_{k}\, z^{k}}{k!\,\Gamma(\alpha k + \beta)},
\qquad \alpha > 0,\ \gamma > 0
$$

burada $(\gamma)_{k} = \Gamma(\gamma + k) / \Gamma(\gamma)$ Pochhammer sembolüdür.

| Fonksiyon | Hesaplama Yöntemi |
|-----------|-------------------|
| İki parametreli Mittag-Leffler | Hibrit OPC + Taylor serisi |
| Üç parametreli Prabhakar | Taylor serisi |

### Sayısal Yöntemler

- **Garrappa'nın Optimal Parabolik Kontur (OPC)** yöntemi — ters Laplace dönüşümü üzerinden yüksek hassasiyetli hesap.
- **Hibrit geçiş** — $\lvert z \rvert \le \tau$ (varsayılan $\tau = 0.1$) için Taylor serisi, aksi halde OPC.
- **Prabhakar için doğru ayrım** — $\gamma \neq 1$ durumunda OPC geçersizdir (dallanma noktaları), bu nedenle seri kullanılır.
- **Karmaşık düzlemde tam destek** — hem gerçek hem sanal argümanlar.

### Özel Durumlar

Aşağıdaki kapalı formlar otomatik olarak doğru hesaplanır:

$$
\begin{aligned}
E_{1,1}(z) &= e^{z} \\
E_{2,1}(z) &= \cosh\!\left(\sqrt{z}\right) \\
E_{1/2,1}(z) &= e^{z^{2}}\,\mathrm{erfc}(-z) \\
E_{\alpha,\beta}(0) &= \frac{1}{\Gamma(\beta)} \\
E_{\alpha,\beta}^{\gamma}(0) &= \frac{1}{\Gamma(\beta)} \quad (\gamma\text{'dan bağımsız})
\end{aligned}
$$

---

## Sayısal Güvenlik Özellikleri

### Float64 Taşma Koruması

- **Dinamik güvenli $\lvert z \rvert$ sınırı** — `_safe_z_limit(α)` fonksiyonu $\alpha$'ya bağlı olarak

$$
\lvert z \rvert_{\max} \approx 0.9 \cdot \bigl(\log(\text{float64}_{\max})\bigr)^{\alpha}
$$

  değerini hesaplar.

- **Rezidü log-büyüklük denetimi** — $\exp(t\,\Re(s))$ çağrısından **önce** taşma riski kontrol edilir.
- **Açıklayıcı hata** — taşma durumunda sessiz `nan` yerine `MittagLefflerError` fırlatılır:

```text
The Mittag-Leffler value is too large to represent in float64:
max log|residue| ≈ 1523.4 exceeds log(float64_max) ≈ 709.8.
|z| = 12.3 exceeds the safe range for α = 0.25 (safe limit ≈ 4.65).
Reduce |z|, or use mpmath for higher-precision arithmetic.
```

### Katastrofik İptal Tespiti

- Küçük $\alpha$ ve büyük $\lvert z \rvert$ kombinasyonlarında tek terimler $\sim e^{\lvert z \rvert^{1/\alpha}}$ mertebesine ulaşırken sonuç $\mathcal{O}(1)$ kalır.
- Modül bu durumlarda test referansı için **adaptif `mpmath.dps`** önerir.
- Kullanıcıya güvenli aralık hakkında bilgi verir.

### Parametre Doğrulaması

`__post_init__` içinde otomatik doğrulama:

| Parametre | Kural |
|-----------|-------|
| `alpha` | $\alpha > 0$, sonlu |
| `beta` | Sonlu reel |
| `gamma` | $\gamma > 0$, sonlu |
| `tolerance` | $0 < \varepsilon < 1$ |
| `series_threshold` | $\ge 0$ |
| `max_series_terms` | Pozitif tam sayı |

---

## API Özellikleri

### İçe Aktarma

```python
from mlfp import MittagLeffler, mittag_leffler, MittagLefflerError
```

### Sınıf Tabanlı API (Tekrar Kullanılabilir)

```python
ml = MittagLeffler(alpha=0.5, beta=1.0, tolerance=1e-14)
print(ml(1.0))                    # gerçek argüman → float
print(ml.evaluate_complex(1.0j))  # karmaşık argüman → complex
```

### Fonksiyon Tabanlı API (Tek Seferlik)

```python
value = mittag_leffler(1.0, alpha=0.5, beta=1.0)
print(value)
```

### Çift Çıktı Modu

| Metot | Girdi | Çıktı | Kullanım |
|-------|-------|-------|----------|
| `__call__` / `evaluate` | Gerçek $z$ | `float64` | Gerçek hesaplamalar |
| `evaluate_complex` | Karmaşık $z$ | `complex128` | Tam karmaşık sonuç |

Karmaşık $z$ gerçek API'ye verilirse **açıklayıcı hata** fırlatılır:

```text
For complex z, use evaluate_complex() to obtain the full complex-valued result.
```

### Vektörleştirme

- NumPy dizileri doğrudan desteklenir.
- Skaler/dizi tutarlılığı garanti edilir.
- Boş girdiler güvenle işlenir.

```python
import numpy as np
z = np.array([0.1, 0.5, 1.0, 2.0], dtype=complex)
result = ml.evaluate_complex(z)   # NDArray[np.complex128]
print(z)
print(result)
```

### Tip İpuçları

- Tam tip ipuçları — `ArrayLike`, `NDArray[np.float64]`, `NDArray[np.complex128]`.
- `from __future__ import annotations` — ileri uyumluluk.
- `ClassVar` — `slots=True` dataclass'ta sabit doğru tanımlı.
- IDE/mypy dostu.

### Hata Yönetimi

```python
class MittagLefflerError(ValueError):
    """Geçersiz parametre veya argümanlar için fırlatılır."""
```

Tüm hatalar tek bir istisna sınıfından türetilir → `try/except` ile kolay yakalama.

---

## Yazılım Mühendisliği Özellikleri

### Performans

- Vektörleştirilmiş hesaplama — OPC integrandı NumPy ile bir kerede değerlendirilir.
- Otomatik düğüm sayısı — $\text{nodes} \le 200$ hedefine göre dinamik ayar.
- Hibrit geçiş — küçük $\lvert z \rvert$ için seri (OPC'den çok daha hızlı).
- Yeniden kullanılabilir sınıf — parametre doğrulaması bir kez yapılır.

### Kod Kalitesi

- `slots=True` dataclass — bellek verimliliği.
- `Final` / `ClassVar` — sabitlerin doğru tanımı.
- Kapsamlı docstring'ler — tüm public API belgelenmiş.
- SPDX lisans başlığı — MIT lisansı otomatik tanınır.
- `__repr__` — hata ayıklama için okunabilir gösterim.

```python
>>> MittagLeffler(alpha=0.5, beta=1.0, gamma=1.0)
MittagLeffler(alpha=0.5, beta=1.0, gamma=1.0, tolerance=1e-15,
              series_threshold=0.1, max_series_terms=100000)
```

### Python 3.15 Uyumluluğu

- `from __future__ import annotations` — yeni tip sözdizimi.
- PEP 604 (`X | Y`) — modern birleşim tipleri.
- `typing.override` — override dekoratörü.
- NumPy 2.x — yeni array API'si.

---

## Doğrulama

### Kapsamlı Test Altyapısı

Modül, **37 pytest testi** ile doğrulanmıştır:

| Test Grubu | Kapsam | Sayı |
|------------|--------|------|
| Özel değerler | $e^{z}$, $\cosh\sqrt{z}$, $\mathrm{erfc}$ | 5 |
| Prabhakar | $\gamma$ indirgemesi, Pochhammer | 6 |
| Sınır koşulları | $z=0$, negatif eksen, taşma | 5 |
| Parametre doğrulama | Geçersiz girdiler | 10 |
| Vektörleştirme | Dizi–skaler tutarlılığı | 3 |
| Seri–OPC | İç tutarlılık | 4 |
| Rekürans | Fonksiyonel ilişki | 3 |
| **Toplam** | | **36 PASS + 1 SKIP** |

### Adaptif mpmath Referansı

Referans değerler, gereken ondalık basamak sayısı

$$
\text{dps} \ge \frac{\lvert z \rvert^{1/\alpha}}{\ln 10} + 40
$$

ile dinamik olarak belirlenerek `mpmath` ile hesaplanır. Bu, küçük $\alpha$ ve büyük $\lvert z \rvert$ durumlarında oluşan **katastrofik iptali** engeller.

### Doğrulanmış Sayısal Sonuçlar

- Tüm geçerli kombinasyonlarda **göreli hata** $< 10^{-13}$.
- Örnek: $\alpha=1.0$, $\beta=1.5$, $\theta=0$ için $\max_{\text{rel}} = 1.49 \times 10^{-15}$.
- Örnek: $\alpha=0.25$, $\beta=0.5$, $\theta=0$ için $\max_{\text{rel}} = 7.23 \times 10^{-14}$.

### Doğrulama Araçları

| Dosya | İşlev |
|-------|-------|
| `test_special_cases.py` | 37 pytest testi |
| `run_comparison.py` | Kapsamlı konsol karşılaştırması |
| `plot_comparison.py` | Değer, mutlak hata, göreli hata panelleri |
| `_ml_helpers.py` | Adaptif referans + hata istatistikleri |

---

## Mevcut Kütüphanelerle Karşılaştırma

| Özellik | `mlfp` | `pymittagleffler` | `khinsen/ml` | `numfracpy` |
|---------|:------:|:-----------------:|:------------:|:-----------:|
| İki parametreli ML | ✅ | ✅ | ✅ | ✅ |
| Prabhakar ($\gamma \neq 1$) | ✅ | ❌ | ❌ | Kısmi |
| Float64 taşma güvenliği | ✅ | ❌ | ❌ | ❌ |
| Açıklayıcı hata mesajları | ✅ | ❌ | ❌ | ❌ |
| Adaptif test referansı | ✅ | ❌ | ❌ | ❌ |
| Çift API (gerçek/karmaşık) | ✅ | ❌ | ✅ | ✅ |
| Vektörleştirme | ✅ | ✅ | ❌ | ✅ |
| NumPy uyumlu tipler | ✅ | ✅ | ❌ | Kısmi |
| Modern Python 3.15 | ✅ | ✅ | ❌ | ❌ |

---

## Sınırlamalar

| Sınırlama | Açıklama | Gelecek Çalışma |
|-----------|----------|-----------------|
| $\alpha = 0$ | Matematiksel olarak $E_{0,\beta}(z) = 1/(\Gamma(\beta)(1-z))$ iyi tanımlı; mevcut OPC uygulaması $\alpha > 0$ varsayar | Özel bir `alpha == 0` dalı eklenebilir |
| Prabhakar büyük $\lvert z \rvert$ | $\gamma \neq 1$ ve büyük $\lvert z \rvert$ için seri yakınsaması yavaşlayabilir | Asimptotik açılım entegrasyonu |
| Matris argüman | Yalnızca skaler/karmaşık $z$ desteklenir | Matris değerli ML fonksiyonu eklenebilir |
| Yüksek duyarlık | float64 ile sınırlı | İsteğe bağlı `mpmath` arka ucu |

---

## Kaynaklar

- **PyPI:** [https://pypi.org/project/mlfp/](https://pypi.org/project/mlfp/)
- **Lisans:** MIT
- **Ana referanslar:**
  - Garrappa, R. (2015). *Numerical evaluation of two and three parameter Mittag-Leffler functions.* SIAM Journal on Numerical Analysis, 53(3), 1350–1369. [doi:10.1137/140971191](https://doi.org/10.1137/140971191)
  - Gorenflo, R., Kilbas, A. A., Mainardi, F., & Rogosin, S. V. (2020). *Mittag-Leffler Functions, Related Topics and Applications.* Springer. [doi:10.1007/978-3-662-61550-8](https://doi.org/10.1007/978-3-662-61550-8)
  - Prabhakar, T. R. (1971). *A singular integral equation with a generalized Mittag-Leffler function in the kernel.* Yokohama Mathematical Journal, 19, 7–15.
  - Mainardi, F., & Garrappa, R. (2017). *On complete monotonicity and asymptotic behaviour of the Prabhakar function.* Fractional Calculus and Applied Analysis, 20(6), 1369–1386. [doi:10.1515/fca-2017-0073](https://doi.org/10.1515/fca-2017-0073)

---

## Lisans

MIT Lisansı altında dağıtılmaktadır. Ayrıntılar için `LICENSE` dosyasına bakın.

---

## Atıf

Bu kütüphaneyi akademik bir çalışmada kullanıyorsanız, aşağıdaki biçimde atıf yapabilirsiniz:

```bibtex
@article{kececi_2026_23605556,
  author       = {Keçeci, Mehmet},
  title        = {Kesirli Kuantum Dinamiği Altında Stratum Mimarisi:
                   NV-Weyl-Majorana Karma Sistemlerinde Topolojik
                   Koruma, Kuantum Hata Düzeltme ve Işınlanma
                  },
  month        = oct,
  year         = 2026,
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.23605556},
  url          = {https://doi.org/10.5281/zenodo.23605556},
}
@article{kececi_2026_23550352,
  author       = {Keçeci, Mehmet},
  title        = {Mittag-Leffler ve Prabhakar Fonksiyonlarının
                   Python'da Güvenli ve Doğru Hesaplanması
                  },
  month        = oct,
  year         = 2026,
  publisher    = {Open Science Articles (OSAs), Zenodo},
  doi          = {10.5281/zenodo.23550352},
  url          = {https://doi.org/10.5281/zenodo.23550352},
}
@software{kececi_2026_23549172,
  author       = {Keçeci, Mehmet},
  title        = {mlfp},
  month        = oct,
  year         = 2026,
  publisher    = {Open Science Articles (OSAs), Zenodo},
  version      = {0.1.1},
  doi          = {10.5281/zenodo.23549172},
  url          = {https://doi.org/10.5281/zenodo.23549172},
}
@software{mlfp2026,
  title  = {mlfp: Mittag-Leffler and Prabhakar Functions for Python},
  author = {Keçeci, Mehmet},
  year   = {2026},
  url    = {https://github.com/WhiteSymmetry/mlfp},
  note   = {MIT License}
}
@software{mlfp2026p,
  title  = {mlfp: Mittag-Leffler and Prabhakar Functions for Python},
  author = {Keçeci, Mehmet},
  year   = {2026},
  url    = {https://pypi.org/project/mlfp/},
  note   = {MIT License}
}
@software{mlfp2026a,
  title  = {mlfp: Mittag-Leffler and Prabhakar Functions for Python},
  author = {Keçeci, Mehmet},
  year   = {2026},
  url    = {https://anaconda.org/channels/bilgi/packages/mlfp},
  note   = {MIT License}
}
@software{mlfp2026d,
  title  = {mlfp: Mittag-Leffler and Prabhakar Functions for Python},
  author = {Keçeci, Mehmet},
  year   = {2026},
  url    = {https://prefix.dev/channels/bilgi/packages/mlfp},
  note   = {MIT License}
}

```

---
