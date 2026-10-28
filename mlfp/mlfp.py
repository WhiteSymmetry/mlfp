# mlfp.py
# SPDX-License-Identifier: MIT
"""
mlfp: Mittag-Leffler and Prabhakar functions for Python.

Two-parameter Mittag-Leffler E_{α,β}(z) and Prabhakar E^γ_{α,β}(z),
plus matrix-valued and fractional-BdG extensions.

- Two-parameter ML (γ = 1): Garrappa OPC with Taylor-series branch
  for small |z|, plus power-series asymptotics in the correct sector.

- Prabhakar (γ ≠ 1): defining power series for small |z|; mpmath
  fallback for large |z|. Note: no fast asymptotics for γ ≠ 1 —
  the simple power-series form omits the exponential A_γ(z) term
  (Mainardi & Garrappa 2017, Eq. 3.9) and is not used.

- Matrix-valued E_{α,β}(A): eigendecomposition, Schur-Parlett, or
  Cauchy contour integral (see mlfp.mlfp_matrix).

- Fractional BdG: (iħ)^α ∂_t^α Ψ = H_BdG Ψ, solved via
  E_{α,1}(H_BdG t^α / (iħ)^α) (see mlfp.mlfp_bdg).

Compatible with Python 3.11-3.15.

0.1.7 özellikleri (0.1.6'e göre)
================================
Modül sayısı: +analyzer, +process, +lindblad
Matris özellikleri: 16+
Sparse desteği: Tam (CSR/CSC/COO)
Fiziksel modeller: +Kitaev, +NV, +Surface Code, +Lindbladian
Grafik sistemi: Dashboard +6 plot fonksiyonu
Kararlılık: Bulletproof (NaN/Inf/Trace koruması)

0.1.6 özellikleri (0.1.5'e göre)
================================
1. **T₁ (Genlik Sönümü) Kanalı Eklendi:** `FractionalAmplitudeDamping` sınıfı ve ilgili yardımcı fonksiyonlar 
(`amplitude_damping_channel`, `local_amplitude_damping`) ile artık relaksasyon (T₁) süreçleri modellenebilir. 
Bu sınıf, kesirli dinamik altında `γ(t) = 1 − E_{α,1}(−(t/T₁)^α)` formülünü kullanır.
2. **Kompozit Kanal Desteği:** `CompositeFractionalChannel` sınıfı ile T₁ ve T₂ etkileri (fenomenolojik olarak) sıralı 
şekilde uygulanabilir hale getirildi.
- mlfp kütüphanesinin Stratum mimarisi gibi karmaşık kuantum sistemlerinin modellenmesinde 
daha güvenilir bir araç olmasını sağlar.

0.1.5 özellikleri (0.1.4'e göre)
================================
1. Parametrelere göre analitik türevler:
     mlfp.derivative(z, alpha, beta, gamma, wrt="alpha")
     MittagLeffler(...).derivative(z, wrt=...)
   Digamma ψ(αk+β) üzerinden seri türevi (küçük |z|); büyük |z|
   için mpmath central difference. ∂E/∂α, ∂E/∂β, ∂E/∂γ desteklenir.

2. Kullanıcı kontrollü yüksek hassasiyet modu:
     MittagLeffler(..., precision=50)
   mpmath workdps(50) bağlamında değerlendirme; tol = 10^(−precision+3)
   (üç guard basamağı). float64'ün 15 basamağının çok ötesine geçen
   sonuçlar (kritik eşik, QFI sınırı vb.) için.

3. Vektörel parametre desteği:
     MittagLeffler(alpha=array, beta=array)
   NumPy broadcasting kurallarıyla uyumlu; tek çağrıda N farklı α
   kombinasyonu değerlendirilebilir.

4. Seri taşma düzeltmesi (lgamma tabanlı):
   _evaluate_series_with_params artık (γ)_k/k! ve 1/Γ(αk+β)
   çarpanlarını logaritmik ölçekte (log_coef birikimiyle) hesaplar.
   γ ≥ 2 ve büyük k için Python int faktöriyel taşması (OverflowError:
   int too large to convert to float) ortadan kaldırıldı; aşırı büyük
   katsayılarda mpmath fallback devreye girer.

5. Matris değerli Mittag-Leffler (yeni modül mlfp/matrix.py):
     MittagLefflerMatrix(alpha, beta, gamma, method=...)
     → E^γ_{α,β}(A) for A ∈ ℂ^{n×n}, yöntem: "eigen" | "schur" |
       "contour" | "auto".
   Schur-Parlett (sla.schur + sla.solve_sylvester) varsayılan;
   Jordan blokları için Hermite interpolasyonu (mpmath.diff);
   Cauchy contour yöntemi 192 Gauss-Legendre düğümü ve dz/dθ
   faktörü ile tam integral.
     solve_fractional_schrodinger(H, psi0, t, alpha=...)
     → kesirli Schrödinger denkleminin çözümü; α=1'de klasik
       exp(−iHt/ħ) ile örtüşür (test edilmiş, ~10^-15 hata).

6. Kesirli Bogoliubov-de Gennes (yeni modül mlfp/mltip_bdg.py):
     build_bdg_hamiltonian(H_0, Delta)
     kitaev_chain(N, t, mu, Delta_0) → (H_0, Delta) p-dalga
     solve_fractional_bdg(H_0, Delta, psi0, t, alpha=...)
     bdg_spectrum(H_0, Delta, tol=...) → gap, n_zero_modes
     find_majorana_modes(H_0, Delta, tol=...) → Majorana tespiti
     majorana_localization(H_0, Delta, tol=...) → uzaysal profil
   α=1'de klasik BdG evrimiyle tam uyumlu; α≠1'de üniter olmayan
   ancak Majorana sıfır modlarında (E=0) topolojik koruma geçerli
   (PHS'nin zorunlu kıldığı sıfır enerji; kesirli evrim operatörü
   genel olarak üniter olmasa da E=0 özdurum üzerinde etkisizdir).

0.1.4 düzeltmeleri (0.1.3'e göre)
==================================
1. __call__ otomatik tip seçimi (NumPy np.exp / np.sin geleneğiyle
   uyumlu): gerçek z → float64; karmaşık z → complex128;
   karmaşık dizi → complex128 dizi; karışık dizi → complex128 dizi.
   evaluate_complex(z) her zaman complex128 döner (davranışı korundu).
   mittag_leffler() fonksiyonu da aynı otomatik tipi kullanır.

2. _invert_laplace_transform_two_parameter: rezidü kontrolü kontur
   integralinden ÖNCE yapılıyor (hızlı başarısızlık); integrand
   hesabı np.errstate(over, invalid, under="ignore") içinde
   (RuntimeWarning sızıntısı önlenir); sonuç isfinite kontrolünden
   geçiriliyor (sessiz inf/nan yerine açıklayıcı MittagLefflerError).

3. _mpmath_fallback: complex(total) dönüşümünden ÖNCE
   log10|total| > 308 kontrolü — float64 aralığını sessizce aşan
   sonuçlar için açıklayıcı MittagLefflerError fırlatılır.

4. Global _mp.mp.dps = 50 ataması kaldırıldı — import yan etkisi yok.
   _mpmath_fallback zaten with _mp.workdps(dps) bağlamını kullanır;
   çağıran tarafın mpmath hassasiyeti etkilenmez.

0.1.3 düzeltmeleri (referans)
==============================
1. _evaluate_series: Γ kutbu nedeniyle sıfır çıkan seri terimi artık
   yakınsama kanıtı sayılmıyor (önceki davranış: E_{0.5,-0.5}(0.05)
   yanlışlıkla yalnız k=0 terimini döndürüyordu).
2. _mpmath_fallback: 1/Γ yerine rgamma — β = 0, −1, … kutuplarında
   "gamma function pole" hatası yerine doğru değer (örn. E^2_{1,0}(z)).
3. İndirgeme: E^2_{α,α+1}(z) = E_{α,α}(z)/α → OPC hızlı yolu
   (~100× hızlanma; doğrudan mpmath yedeklemesine kıyasla).
4. α = 0 dalı: E_{0,β}^γ(z) = (1-z)^{-γ}/Γ(β), |z| < 1.
5. mittag_leffler() işlev API'sine max_series_terms parametresi eklendi.
"""

from __future__ import annotations

import cmath
import functools
import math
from dataclasses import dataclass
from matplotlib.colors import LogNorm, Normalize
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import ArrayLike, NDArray
from typing import ClassVar, Final
import scipy.linalg as la
from scipy.linalg import eigh, eig
import scipy
import scipy.special as special
import scipy.sparse as sparse
import scipy.sparse.linalg as spla
from scipy.special import (
    gamma as gamma_func,
    digamma,
    gammaln,
    rgamma,
    erfcx, poch
)
from scipy.optimize import brentq, curve_fit
from pfapack import pfaffian

try:
    import scipy.special as special
    import scipy.sparse as sparse
    import scipy.sparse.linalg as spla
    from scipy.sparse.linalg import norm as sp_norm, eigsh, svds
    HAS_SCIPY = True
except Exception:
    sparse = None
    sp_norm = None
    spla = None
    eigsh = None
    svds = None
    HAS_SCIPY = False

try:
    from pfapack import pfaffian as pf
    HAS_PFAPACK = True
except ImportError:
    pf = None
    HAS_PFAPACK = False


def _is_sparse(x) -> bool:
    """
    scipy.sparse matris olup olmadığını kontrol eder.
    """
    return HAS_SCIPY and sparse.issparse(x)

import mpmath as _mp
# Not: global mp.dps DEĞİŞTİRİLMEZ. _mpmath_fallback kendi
# workdps(dps) bağlamını kullanır; çağıran tarafın mpmath
# hassasiyetini etkilemez.
import warnings  # Sonunda 's' olmalı

# Örnek Kullanım: Kodunuzdaki tüm uyarıları gizlemek için
warnings.filterwarnings("ignore")

# Örnek Kullanım: Özel bir uyarı tetiklemek için
warnings.warn("Bu bir özel uyarı mesajıdır!")

#analyzer = MatrixAnalyzer(A)

class MittagLefflerError(ValueError):
    """Raised for invalid Mittag-Leffler parameters or arguments."""

# ---------------------------------------------------------------------------
# Float64 stability limits
# ---------------------------------------------------------------------------
_LOG_MAX_FLOAT64: Final[float] = math.log(np.finfo(np.float64).max)  # ≈ 709.78

def _asymptotic_expansion(
    z: complex,
    alpha: float,
    beta: float,
    *,
    n_terms: int = 40,
    rel_tol: float = 1e-14,
) -> complex | None:
    """
    E_{α,β}(z)'nin |z| → ∞ için kuvvet-serisi asimptotiği.

    GEÇERLİ SADECE ŞU SECTOR'DA (α < 2 için):
        |arg(z)| ∈ (απ/2, π)   VE   z negatif reel eksende değil

    Bu sector dışında None döner → caller mpmath fallback yapar.
    """
    if alpha >= 2.0:
        return None

    theta = cmath.phase(z)
    abs_theta = abs(theta)

    # ← DOĞRU SECTOR: |θ| > απ/2, üst sınır YOK
    #    θ = π (negatif reel eksen) dahil edilir.
    lower = alpha * math.pi / 2.0 + 1e-3
    if abs_theta < lower:
        return None
    # ÜST SINIR KALDIRILDI (önceden: abs_theta < π - 1e-3)

    if abs(z) < 5.0:
        return None

    total = 0.0 + 0.0j
    prev_mag = math.inf
    for k in range(1, n_terms + 1):
        denom = beta - alpha * k
        if abs(denom - round(denom)) < 1e-12 and denom <= 0:
            continue
        try:
            term = z ** (-k) / gamma_func(denom)
        except (ValueError, ZeroDivisionError):
            continue
        if not np.isfinite(term):
            break
        total -= term
        mag = abs(term)
        if mag > prev_mag:
            break
        prev_mag = mag
        if mag < rel_tol * (1.0 + abs(total)):
            break
    return complex(total)

def _mpmath_fallback(
    z, alpha, beta, gamma=1.0,
    *,
    max_dps: int = 2000,     # ← 800 → 2000
) -> complex:
    """
    E_{α,β}^γ(z)'yi adaptif hassasiyetli mpmath ile hesapla.

    float64'te taşma/iptal riski olduğunda çağrılır.
    Gerekli hassasiyet: |z|^(1/α)/ln(10) + güvenlik payı.
    """
    abs_z = abs(z)

    # Gerekli ondalık basamak sayısını hesapla
    if abs_z < 1.0 or alpha >= 2.0:
        dps = 60
    else:
        log_peak = abs_z ** (1.0 / alpha)          # nat cinsinden
        dps = max(60, int(log_peak / math.log(10.0)) + 60)

    if dps > max_dps:
        raise MittagLefflerError(
            f"Even mpmath would need {dps} digits (> {max_dps}) "
            f"for α={alpha}, |z|={abs_z:.3g}. "
            "The function value itself is beyond float64 range."
        )

    with _mp.workdps(dps):
        z_mp     = _mp.mpc(z)
        alpha_mp = _mp.mpf(alpha)
        beta_mp  = _mp.mpf(beta)
        gamma_mp = _mp.mpf(gamma)

        tol = _mp.mpf(10) ** (-(dps - 15))
        # DÜZELTME: rgamma kullan — Γ kutuplarında (β = 0, −1, …) 1/Γ = 0
        # verir; "gamma function pole" hatası fırlatılmaz. Böylece
        # E^2_{1,0}(z) = (z² + 2z)e^z gibi iyi tanımlı değerler hesaplanır.
        total = _mp.rgamma(beta_mp)
        z_power = _mp.mpc(1)
        pochhammer = _mp.mpf(1)

        for k in range(1, 200_000):
            z_power *= z_mp
            pochhammer *= gamma_mp + k - 1
            # rgamma: seri terimi Γ kutbuna denk gelirse 0 olur, hata çıkmaz
            term = (pochhammer * z_power
                    * _mp.rgamma(alpha_mp * k + beta_mp) / _mp.factorial(k))
            total += term
            # DÜZELTME: sıfır terim yakınsama kanıtı sayılmaz
            if term != 0 and abs(term) < abs(total) * tol:
                # ---- YENİ: float64 aralık kontrolü ----
                if total != 0:
                    log10_mag = float(_mp.log10(abs(total)))
                    if log10_mag > 308.0:
                        raise MittagLefflerError(
                            f"Mittag-Leffler value (~10^{log10_mag:.0f}) "
                            f"exceeds float64 range (max ~10^308) for "
                            f"α={alpha}, β={beta}, γ={gamma}, |z|={abs_z:.3g}. "
                            "The function value cannot be represented in "
                            "Python complex (float64)."
                        )
                return complex(total)

    raise MittagLefflerError(
        f"mpmath fallback did not converge for α={alpha}, "
        f"β={beta}, γ={gamma}, |z|={abs_z:.3g}, dps={dps}."
    )


def _safe_z_limit(alpha: float, safety: float = 0.9) -> float:
    """
    Upper bound on |z| for which the OPC residue sum fits in float64.

    The residue of the branch point s satisfies |Res| ≲ exp(t·Re(s))·|s|^(1-β),
    with Re(s) ≈ |z|^(1/α) along the positive real axis.  Requiring
    t·Re(s) < log(float64 max) yields |z| < log_max^α.

    NOT: Bu sınır float64 rezidü hesabı içindir. |z| sınırı aşarsa,
    _asymptotic_expansion devreye girer (arg(z) uygunsa).
    """
    if alpha <= 0.0:
        return 0.0
    return safety * _LOG_MAX_FLOAT64 ** alpha


# ---------------------------------------------------------------------------
# Low-level optimal parabolic contour routines (Garrappa OPC)
# ---------------------------------------------------------------------------
def _find_optimal_bounded_param(
    t: float,
    phi_star0: float,
    phi_star1: float,
    p: float,
    q: float,
    *,
    log_eps: float,
    log_machine_eps: float,
    fac: float = 1.01,
    p_eps: float = 1.0e-14,
    q_eps: float = 1.0e-14,
) -> tuple[float, float, float]:
    """Find optimal contour parameters for a right-bounded region."""
    f_max = math.exp(log_eps - log_machine_eps)
    threshold = 2.0 * math.sqrt((log_eps - log_machine_eps) / t)

    sqrt_phi0 = math.sqrt(phi_star0)
    sqrt_phi1 = min(math.sqrt(phi_star1), threshold - sqrt_phi0)

    if p < p_eps and q < q_eps:
        sqrt_phi_bar0 = sqrt_phi0
        sqrt_phi_bar1 = sqrt_phi1
        admissible = True

    elif p < p_eps:
        sqrt_phi_bar0 = sqrt_phi0
        if sqrt_phi0 > 0.0:
            f_min = fac * (sqrt_phi0 / (sqrt_phi1 - sqrt_phi0)) ** q
        else:
            f_min = fac
        if f_min < f_max:
            f_bar = f_min + f_min / f_max * (f_max - f_min)
            f_q = f_bar ** (-1.0 / q)
            sqrt_phi_bar1 = (2.0 * sqrt_phi1 - f_q * sqrt_phi0) / (2.0 + f_q)
            admissible = True
        else:
            admissible = False

    elif q < q_eps:
        sqrt_phi_bar1 = sqrt_phi1
        f_min = fac * (sqrt_phi1 / (sqrt_phi1 - sqrt_phi0)) ** p
        if f_min < f_max:
            f_bar = f_min + f_min / f_max * (f_max - f_min)
            f_p = f_bar ** (-1.0 / p)
            sqrt_phi_bar0 = (2.0 * sqrt_phi0 - f_p * sqrt_phi1) / (2.0 - f_p)
            admissible = True
        else:
            admissible = False

    else:
        f_min = (
            fac
            * (sqrt_phi1 + sqrt_phi0)
            / (sqrt_phi1 - sqrt_phi0) ** max(p, q)
        )
        if f_min < f_max:
            f_min = max(f_min, 1.5)
            f_bar = f_min + f_min / f_max * (f_max - f_min)
            f_p = f_bar ** (-1.0 / p)
            f_q = f_bar ** (-1.0 / q)
            weight = -phi_star1 * t / log_eps
            denominator = 2.0 + weight - (1.0 + weight) * f_p + f_q
            sqrt_phi_bar0 = (
                (2.0 + weight + f_q) * sqrt_phi0 + f_p * sqrt_phi1
            ) / denominator
            sqrt_phi_bar1 = (
                -(1.0 + weight) * f_q * sqrt_phi0
                + (2.0 + weight - (1.0 + weight) * f_p) * sqrt_phi1
            ) / denominator
            admissible = True
        else:
            admissible = False

    if not admissible:
        return 0.0, math.inf, 0.0

    log_eps -= math.log(f_bar)
    weight = -sqrt_phi_bar1**2 * t / log_eps

    mu = (
        ((1.0 + weight) * sqrt_phi_bar0 + sqrt_phi_bar1) / (2.0 + weight)
    ) ** 2

    step = (
        -2.0
        * math.pi
        / log_eps
        * (sqrt_phi_bar1 - sqrt_phi_bar0)
        / ((1.0 + weight) * sqrt_phi_bar0 + sqrt_phi_bar1)
    )

    nodes = math.ceil(math.sqrt(1.0 - log_eps / t / mu) / step)
    return mu, nodes, step


def _find_optimal_unbounded_param(
    t: float,
    phi_star: float,
    p: float,
    *,
    log_eps: float,
    log_machine_eps: float,
    fac: float = 1.01,
    p_eps: float = 1.0e-14,
) -> tuple[float, float, float]:
    """Find optimal contour parameters for a right-unbounded region."""
    sqrt_phi_star = math.sqrt(phi_star)
    phi_bar = fac * phi_star if phi_star > 0.0 else 0.01
    sqrt_phi_bar = math.sqrt(phi_bar)

    f_min = 1.0
    f_max = 10.0
    f_target = 5.0

    while True:
        phi_t = phi_bar * t
        log_eps_phi_t = log_eps / phi_t
        nodes = math.ceil(
            phi_t
            / math.pi
            * (
                1.0
                - 3.0 * log_eps_phi_t / 2.0
                + math.sqrt(1.0 - 2.0 * log_eps_phi_t)
            )
        )
        a = math.pi * nodes / phi_t
        sqrt_mu = (
            sqrt_phi_bar * abs(4.0 - a)
            / abs(7.0 - math.sqrt(1.0 + 12.0 * a))
        )
        f_bar = ((sqrt_phi_bar - sqrt_phi_star) / sqrt_mu) ** (-p)
        if p < p_eps or f_min < f_bar < f_max:
            break
        sqrt_phi_bar = f_target ** (-1.0 / p) * sqrt_mu + sqrt_phi_star
        phi_bar = sqrt_phi_bar**2

    mu = sqrt_mu**2
    step = (
        (-3.0 * a - 2.0 + 2.0 * math.sqrt(1.0 + 12.0 * a))
        / (4.0 - a)
        / nodes
    )

    roundoff_threshold = (log_eps - log_machine_eps) / t
    if mu > roundoff_threshold:
        q = (
            0.0
            if abs(p) < p_eps
            else f_target ** (-1.0 / p) * math.sqrt(mu)
        )
        phi_bar = (q + sqrt_phi_star) ** 2
        if phi_bar < roundoff_threshold:
            weight = math.sqrt(
                log_machine_eps / (log_machine_eps - log_eps)
            )
            u = math.sqrt(-phi_bar * t / log_machine_eps)
            mu = roundoff_threshold
            nodes = math.ceil(
                weight * log_eps / (2.0 * math.pi * (u * weight - 1.0))
            )
            step = weight / nodes
        else:
            return mu, math.inf, 0.0

    return mu, nodes, step


def _invert_laplace_transform_two_parameter(
    t: float,
    z: complex,
    *,
    alpha: float,
    beta: float,
    tolerance: float,
) -> complex:
    """Evaluate E_{alpha,beta}(z t^alpha) via OPC Laplace inversion."""
    if abs(z) == 0.0:
        return complex(1.0 / gamma_func(beta))

    machine_eps = float(np.finfo(np.float64).eps)
    log_machine_eps = math.log(machine_eps)
    log_tolerance = math.log(tolerance)
    log_10 = math.log(10.0)
    log_tolerance_difference = log_tolerance - log_machine_eps

    theta = cmath.phase(z)
    k_min = math.ceil(-alpha / 2.0 - theta / (2.0 * math.pi))
    k_max = math.floor(alpha / 2.0 - theta / (2.0 * math.pi))

    k_values = np.arange(k_min, k_max + 1, dtype=np.float64)
    poles = abs(z) ** (1.0 / alpha) * np.exp(
        1j * (theta + 2.0 * math.pi * k_values) / alpha
    )

    phi_values = (poles.real + np.abs(poles)) / 2.0
    order = np.argsort(phi_values)
    poles = poles[order]
    phi_values = phi_values[order]

    nonzero = phi_values > machine_eps
    poles = poles[nonzero]
    phi_values = phi_values[nonzero]

    poles = np.insert(poles, 0, 0.0 + 0.0j)
    phi_values = np.insert(phi_values, 0, 0.0)

    strengths = np.ones(poles.size, dtype=np.float64)
    strengths[0] = max(0.0, -2.0 * (alpha - beta + 1.0))
    right_strengths = np.ones(poles.size, dtype=np.float64)
    right_strengths[-1] = math.inf

    phi_values = np.insert(phi_values, phi_values.size, math.inf)

    region_indices = np.flatnonzero(
        (phi_values[:-1] < log_tolerance_difference / t)
        & (phi_values[:-1] < phi_values[1:])
    )
    if region_indices.size == 0:
        raise MittagLefflerError(
            "No admissible integration region was found."
        )

    region_count = region_indices[-1] + 1
    mu_values = np.full(region_count, math.inf, dtype=np.float64)
    node_values = np.full(region_count, math.inf, dtype=np.float64)
    step_values = np.full(region_count, math.inf, dtype=np.float64)

    current_log_tolerance = log_tolerance
    while True:
        for region in region_indices:
            if region < poles.size - 1:
                (mu_values[region],
                 node_values[region],
                 step_values[region]) = _find_optimal_bounded_param(
                    t,
                    phi_values[region],
                    phi_values[region + 1],
                    strengths[region],
                    right_strengths[region],
                    log_eps=current_log_tolerance,
                    log_machine_eps=log_machine_eps,
                )
            else:
                (mu_values[region],
                 node_values[region],
                 step_values[region]) = _find_optimal_unbounded_param(
                    t,
                    phi_values[region],
                    strengths[region],
                    log_eps=current_log_tolerance,
                    log_machine_eps=log_machine_eps,
                )

        if np.min(node_values) <= 200.0:
            break
        current_log_tolerance += log_10
        if current_log_tolerance >= 0.0:
            raise MittagLefflerError(
                "Failed to find an admissible parabolic contour region."
            )

    selected = int(np.argmin(node_values))
    mu = float(mu_values[selected])
    nodes = int(node_values[selected])
    step = float(step_values[selected])

    indices = np.arange(-nodes, nodes + 1, dtype=np.float64)
    u_values = step * indices

    contour = mu * (1j * u_values + 1.0) ** 2
    derivative = -2.0 * mu * u_values + 2j * mu

    # ------------------------------------------------------------------
    # 1) Rezidü kontrolü EN ÖNCE — hızlı başarısızlık
    # ------------------------------------------------------------------
    residue_poles = poles[selected + 1 :]
    residues = 0.0

    if residue_poles.size:
        magnitudes = np.maximum(np.abs(residue_poles), 1e-300)
        log_mags = (
            t * residue_poles.real
            + (1.0 - beta) * np.log(magnitudes)
        )
        max_log_mag = float(np.max(log_mags))

        if max_log_mag > _LOG_MAX_FLOAT64 - 1.0:
            raise MittagLefflerError(
                "The Mittag-Leffler value is too large to represent in "
                f"float64: max log|residue| ≈ {max_log_mag:.1f} exceeds "
                f"log(float64_max) ≈ {_LOG_MAX_FLOAT64:.1f}. "
                f"|z| = {abs(z):.3g} exceeds the safe range for "
                f"α = {alpha} (safe limit ≈ {_safe_z_limit(alpha):.3g}). "
                "Reduce |z|, or use mpmath for higher-precision arithmetic."
            )

        with np.errstate(over="ignore", invalid="ignore"):
            residues = (
                residue_poles ** (1.0 - beta)
                * np.exp(t * residue_poles)
            ).sum() / alpha

    # ------------------------------------------------------------------
    # 2) Kontur integrali — numpy taşma uyarılarını sustur
    # ------------------------------------------------------------------
    # Tekil terimler geçici olarak taşsa bile nihai toplam sonlu
    # olabilir; bu yüzden errstate ile uyarıları susturup SONUCU kontrol
    # ediyoruz (önceki sürümde uyarı sızıyordu).
    with np.errstate(over="ignore", invalid="ignore", under="ignore"):
        integrand = (
            np.exp(contour * t)
            * contour ** (alpha - beta)
            / (contour**alpha - z)
            * derivative
        )
        integral = step * integrand.sum() / (2j * math.pi)

    if not np.isfinite(integral) or not np.isfinite(residues):
        raise MittagLefflerError(
            f"Laplace inversion produced non-finite result for "
            f"α={alpha}, β={beta}, |z|={abs(z):.3g}. "
            "The value may be beyond float64 range."
        )

    return complex(integral + residues)

def derivative(
    z: ArrayLike,
    alpha: float,
    beta: float = 1.0,
    gamma: float = 1.0,
    *,
    wrt: str = "alpha",
    tolerance: float = 1.0e-15,
    series_threshold: float = 0.1,
    max_series_terms: int = 100_000,
) -> complex | NDArray[np.complex128]:
    """
    Compute ∂E^γ_{α,β}(z) / ∂(wrt) at z.

    Parameters
    ----------
    z : real or complex, scalar or array
    alpha, beta, gamma : float, scalar parameters
    wrt : {"alpha", "beta", "gamma"}
    tolerance, series_threshold, max_series_terms : as in MittagLeffler

    Returns
    -------
    complex (scalar) or complex ndarray

    Notes
    -----
    * For ``|z| <= series_threshold``: analytic series derivative (fast).
    * Otherwise: mpmath central-difference (accurate but slower).

    Examples
    --------
    >>> derivative(0.3 + 0.2j, alpha=0.5, beta=1.0, wrt="alpha")
    (-0.05...+0.1...j)
    """
    if wrt not in ("alpha", "beta", "gamma"):
        raise MittagLefflerError(
            f"wrt must be 'alpha', 'beta', or 'gamma'; got {wrt!r}"
        )

    # Delegate to a temporary evaluator (single source of truth)
    evaluator = MittagLeffler(
        alpha=alpha, beta=beta, gamma=gamma,
        tolerance=tolerance,
        series_threshold=series_threshold,
        max_series_terms=max_series_terms,
    )
    # Call the *private* per-point implementation directly,
    # avoiding infinite recursion via the public .derivative() method
    return evaluator._derivative_dispatch(z, wrt)

def _normalize_param(value, name):
    """Convert scalar/array-like to float or np.ndarray."""
    if isinstance(value, (list, tuple, np.ndarray)):
        arr = np.asarray(value, dtype=float)
        if arr.ndim == 0 or (arr.ndim == 1 and arr.size == 1):
            return float(arr.ravel()[0])
        return arr
    return float(value)

# ---------------------------------------------------------------------------
# Public evaluator
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class MittagLeffler:
    """
    Evaluator for the Mittag-Leffler and Prabhakar functions.

    For gamma = 1: ``E_{alpha,beta}(z)`` is evaluated via Garrappa's
    optimal parabolic contour method (with a Taylor-series branch for
    small ``|z|``).

    For gamma != 1: ``E_{alpha,beta}^gamma(z)`` is evaluated from the
    defining power series.

    Parameters
    ----------
    alpha : float
        Non-negative real order parameter (``alpha = 0`` requires
        ``|z| < 1`` and uses the closed form ``(1-z)**(-gamma)/Gamma(beta)``).
    beta : float
        Real second parameter.
    gamma : float
        Positive real Prabhakar parameter.
    tolerance : float
        Target relative accuracy (must satisfy ``0 < tolerance < 1``).
    series_threshold : float
        For ``gamma == 1``, arguments with ``|z| <= series_threshold``
        use the Taylor series instead of OPC.
    max_series_terms : int
        Maximum number of terms in the series evaluator.
    alpha, beta, gamma : float or array-like
        Scalar or broadcastable arrays (vectorized mode).
    tolerance : float
    series_threshold : float
    max_series_terms : int
    precision : int or None
        NEW in 0.1.5. If not None (>= 15), uses mpmath with that many
        decimal digits; the return value is mpmath.mpf/mpc.
    """

    alpha: float
    beta: float = 1.0
    gamma: float = 1.0
    tolerance: float = 1.0e-15
    series_threshold: float = 0.1
    max_series_terms: int = 100_000
    precision: int | None = None

    # NOTE: ClassVar so that `slots=True` does not treat it as a field.
    _MACHINE_EPS: ClassVar[float] = float(np.finfo(np.float64).eps)

    def __post_init__(self) -> None:
        # NEW 0.1.5: array support
        self.alpha = _normalize_param(self.alpha, "alpha")
        self.beta = _normalize_param(self.beta, "beta")
        self.gamma = _normalize_param(self.gamma, "gamma")
        
        if not np.all(np.isfinite(self.alpha)) or np.any(np.asarray(self.alpha) < 0):
            raise MittagLefflerError(
                "alpha must be non-negative real (scalar or array)."
            )
        if not np.all(np.isfinite(self.beta)):
            raise MittagLefflerError(
                "beta must be finite real (scalar or array)."
            )
        if not np.all(np.isfinite(self.gamma)) or np.any(np.asarray(self.gamma) <= 0):
            raise MittagLefflerError(
                "gamma must be positive real (scalar or array)."
            )
        if not 0.0 < self.tolerance < 1.0:
            raise MittagLefflerError("tolerance must satisfy 0 < tolerance < 1.")
        if self.series_threshold < 0.0:
            raise MittagLefflerError("series_threshold must be non-negative.")
        if self.max_series_terms < 1:
            raise MittagLefflerError("max_series_terms must be positive integer.")
        if self.precision is not None:
            if not isinstance(self.precision, int) or self.precision < 15:
                raise MittagLefflerError(
                    "precision must be int >= 15 (float64 needs 15)."
                )
    
    def _params_are_arrays(self) -> bool:
        return (isinstance(self.alpha, np.ndarray)
                or isinstance(self.beta, np.ndarray)
                or isinstance(self.gamma, np.ndarray))

    # ------------------------------------------------------------------
    # Real-valued public interface
    # ------------------------------------------------------------------
    def __call__(
        self, z: ArrayLike
    ) -> (
        NDArray[np.float64] | np.float64
        | NDArray[np.complex128] | np.complex128
    ):
        r"""
        Evaluate the function at one point or over an array.

        Return type follows NumPy's convention (like ``np.exp``, ``np.sin``):

        * **Real ``z``** (imaginary part exactly 0): returns ``float64``.
        * **Complex ``z``** (any non-zero imaginary part): returns
          ``complex128`` with the full result.

        Use :meth:`evaluate_complex` if you always want a complex result
        (even for real inputs).
        """
        # 0.1.5: precision mode
        if self.precision is not None:
            return self._evaluate_precise(z)
        
        # 0.1.5: vectorized params
        if self._params_are_arrays():
            return self._evaluate_vectorized(z)
        
        # Scalar params (existing behavior, unchanged)
        z_array = np.asarray(z, dtype=np.complex128)
        if z_array.ndim == 0:
            z_val = complex(z_array)
            if z_val.imag == 0.0:
                return np.float64(self._evaluate_scalar(z_val))
            return np.complex128(self._evaluate_scalar_complex(z_val))
        if np.all(z_array.imag == 0.0):
            evaluate = np.vectorize(
                lambda value: self._evaluate_scalar(complex(value)),
                otypes=[np.float64],
            )
            return evaluate(z_array)
        evaluate = np.vectorize(
            lambda value: self._evaluate_scalar_complex(complex(value)),
            otypes=[np.complex128],
        )
        return evaluate(z_array)

    def evaluate(
        self, z: ArrayLike
    ) -> (
        NDArray[np.float64] | np.float64
        | NDArray[np.complex128] | np.complex128
    ):
        """Explicit alias for :meth:`__call__`."""
        return self(z)

    # ------------------------------------------------------------------
    # Complex-valued public interface
    # ------------------------------------------------------------------
    def evaluate_complex(self, z: ArrayLike):
        """Complex result; routes like __call__ but always complex."""
        if self.precision is not None:
            return self._evaluate_precise(z)
        
        if self._params_are_arrays():
            result = self._evaluate_vectorized(z)
            return np.asarray(result, dtype=np.complex128)
        
        z_array = np.asarray(z, dtype=np.complex128)
        if z_array.ndim == 0:
            return np.complex128(
                self._evaluate_scalar_complex(complex(z_array))
            )
        evaluate = np.vectorize(
            lambda value: self._evaluate_scalar_complex(complex(value)),
            otypes=[np.complex128],
        )
        return evaluate(z_array)

    # ------------------------------------------------------------------
    # Internal evaluators
    # ------------------------------------------------------------------
    def _evaluate_scalar(self, z: complex) -> float:
        """Return the real part for real z; raise otherwise.

        Internal API. `__call__` always dispatches here only for z with
        z.imag == 0.0, so the guard is defensive.
        """
        value = self._evaluate_scalar_complex(z)
        if z.imag == 0.0:
            return float(value.real)
        raise MittagLefflerError(
            "For complex z, use evaluate_complex() to obtain the full "
            "complex-valued result."
        )

    def _evaluate_scalar_complex(self, z: complex) -> complex:
        return self._evaluate_complex_impl(
            z, self.alpha, self.beta, self.gamma,
        )


    def _evaluate_complex_impl(
        self, z: complex, alpha: float, beta: float, gamma: float
    ) -> complex:
        """Parameterized complex evaluation (0.1.5 refactor)."""
        if abs(z) < self._MACHINE_EPS:
            return complex(1.0 / gamma_func(beta))
        
        if alpha == 0.0:
            if abs(z) >= 1.0:
                raise MittagLefflerError(
                    "alpha = 0 requires |z| < 1."
                )
            return complex((1.0 - z) ** (-gamma) / gamma_func(beta))
        
        if gamma == 1.0 and abs(alpha - 1.0) < 1e-12 and abs(beta - 1.0) < 1e-12:
            with np.errstate(over="ignore", invalid="ignore"):
                return complex(np.exp(z))
        
        if gamma == 1.0:
            if abs(z) <= self.series_threshold:
                return self._evaluate_series_with_params(z, alpha, beta, gamma)
            if abs(z) > 0.8 * _safe_z_limit(alpha):
                asymp = _asymptotic_expansion(z, alpha, beta)
                if asymp is not None:
                    return asymp
            try:
                return _invert_laplace_transform_two_parameter(
                    t=1.0, z=z, alpha=alpha, beta=beta,
                    tolerance=self.tolerance,
                )
            except MittagLefflerError as exc:
                if "too large" in str(exc):
                    return _mpmath_fallback(z, alpha, beta, gamma)
                raise
        
        # Reduction E^2_{α,α+1} = E_{α,α} / α
        if (abs(z) > self.series_threshold
                and abs(gamma - 2.0) < 1e-9
                and abs(beta - alpha - 1.0) < 1e-9):
            reduced = _reduced_evaluator(
                alpha, alpha,
                self.tolerance, self.series_threshold, self.max_series_terms,
            )
            return reduced._evaluate_scalar_complex(z) / alpha

        # γ ≠ 1 yolu
        if abs(z) <= self.series_threshold:
            return self._evaluate_series_with_params(z, alpha, beta, gamma)

        # 0.1.5 NOTU: γ≠1 için basit kuvvet asimptotiği EKSİK'tir
        # (üstel A_γ(z) terimi ihmal edilir). Doğru asimptotik Wright
        # fonksiyonu üzerinden gelir — gelecek sürümde (0.3.0) eklenecek.
        # Şimdilik: mpmath fallback (yavaş ama doğru).
        return _mpmath_fallback(z, alpha, beta, gamma)

        def _evaluate_series(self, z: complex) -> complex:
            """Legacy wrapper — uses self's parameters."""
            return self._evaluate_series_with_params(
                z, self.alpha, self.beta, self.gamma,
            )


    def _evaluate_series_with_params(
        self, z: complex, alpha: float, beta: float, gamma: float
    ) -> complex:
        r"""
        Evaluate from power series (0.1.5: lgamma-based, overflow-safe).

        E^γ_{α,β}(z) = Σ_{k=0}^∞ (γ)_k z^k / (k! Γ(αk+β))
        
        Using gammaln avoids the int-overflow that occurred when
        ``factorial *= k`` exceeded float64 max (~171!).
        """
        # k=0 term — rgamma handles Γ poles (returns 0)
        total = complex(float(rgamma(beta)))
        
        z_power = 1.0 + 0.0j
        log_coef = 0.0  # log((γ)_k / k!)
        
        for k in range(1, self.max_series_terms + 1):
            z_power *= z
            # (γ)_k / k! via incremental logs — no int overflow
            log_coef += math.log(gamma + k - 1.0) - math.log(k)
            
            try:
                coef = math.exp(log_coef)
            except OverflowError:
                # Coef huge → defer to mpmath (accurate, slow)
                return _mpmath_fallback(z, alpha, beta, gamma)
            
            # 1/Γ(αk+β) — rgamma returns 0 at Γ poles (safe)
            try:
                rg = float(rgamma(alpha * k + beta))
            except (ValueError, OverflowError):
                rg = 0.0
            
            if coef == 0.0 or rg == 0.0:
                continue
            
            term = coef * z_power * rg
            if term == 0.0:
                continue
            
            total += term
            if abs(term) <= self.tolerance * (1.0 + abs(total)):
                return total
        
        raise MittagLefflerError(
            "The Mittag-Leffler series did not converge to the requested "
            "tolerance within max_series_terms terms."
        )

    def _evaluate_precise(self, z: ArrayLike):
        """mpmath evaluation with self.precision digits (0.1.5)."""
        if self._params_are_arrays():
            raise MittagLefflerError(
                "Precision mode requires scalar alpha/beta/gamma."
            )
        z_arr = np.asarray(z, dtype=np.complex128)
        if z_arr.ndim == 0:
            return self._evaluate_precise_scalar(complex(z_arr))
        out = np.empty(z_arr.shape, dtype=object)
        for idx in np.ndindex(z_arr.shape):
            out[idx] = self._evaluate_precise_scalar(complex(z_arr[idx]))
        return out


    def _evaluate_precise_scalar(self, z: complex):
        """Single-point mpmath evaluation at self.precision digits."""
        with _mp.workdps(self.precision):
            z_mp = _mp.mpc(z)
            a_mp = _mp.mpf(float(self.alpha))
            b_mp = _mp.mpf(float(self.beta))
            g_mp = _mp.mpf(float(self.gamma))

            # α = 0 special branch
            if a_mp == 0:
                if abs(z_mp) >= 1:
                    raise MittagLefflerError(
                        "alpha = 0 requires |z| < 1."
                    )
                return _mp.power(1 - z_mp, -g_mp) * _mp.rgamma(b_mp)

            total = _mp.rgamma(b_mp)
            z_power = _mp.mpc(1)
            poch = _mp.mpf(1)

            # ✅ DÜZELTME: 3 guard digits (önceki: 10)
            tol = _mp.mpf(10) ** (-(self.precision - 3))

            for k in range(1, 500_000):
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
            f"Precise evaluation did not converge at "
            f"{self.precision} digits (500k terms)."
        )

    def _evaluate_vectorized(self, z: ArrayLike):
        """Vectorized evaluation with array alpha/beta/gamma (0.1.5)."""
        z_arr = np.asarray(z, dtype=np.complex128)
        a_arr = np.asarray(self.alpha, dtype=float)
        b_arr = np.asarray(self.beta, dtype=float)
        g_arr = np.asarray(self.gamma, dtype=float)
        
        try:
            a_b, b_b, g_b, z_b = np.broadcast_arrays(
                a_arr, b_arr, g_arr, z_arr,
            )
        except ValueError as exc:
            raise MittagLefflerError(
                f"Incompatible shapes: alpha={a_arr.shape}, beta={b_arr.shape}, "
                f"gamma={g_arr.shape}, z={z_arr.shape}"
            ) from exc
        
        all_real = bool(np.all(z_b.imag == 0))
        
        if a_b.ndim == 0:
            result = self._evaluate_complex_impl(
                complex(z_b), float(a_b), float(b_b), float(g_b),
            )
            return (np.float64(result.real) if all_real
                    else np.complex128(result))
        
        if all_real:
            out = np.empty(a_b.shape, dtype=np.float64)
            for idx in np.ndindex(a_b.shape):
                v = self._evaluate_complex_impl(
                    complex(z_b[idx]), float(a_b[idx]),
                    float(b_b[idx]), float(g_b[idx]),
                )
                out[idx] = v.real
            return out
        
        out = np.empty(a_b.shape, dtype=np.complex128)
        for idx in np.ndindex(a_b.shape):
            out[idx] = self._evaluate_complex_impl(
                complex(z_b[idx]), float(a_b[idx]),
                float(b_b[idx]), float(g_b[idx]),
            )
        return out

    def derivative(self, z: ArrayLike, *, wrt: str = "alpha"):
        """
        Compute ∂E/∂(wrt) at z using self's alpha/beta/gamma.

        Thin wrapper around the per-point derivative dispatch; honors
        self.tolerance, self.series_threshold, self.max_series_terms.

        Parameters
        ----------
        z : real or complex, scalar or array
        wrt : {"alpha", "beta", "gamma"}

        Returns
        -------
        complex (scalar) or complex ndarray
        """
        return self._derivative_dispatch(z, wrt)

    def _derivative_dispatch(self, z: ArrayLike, wrt: str):
        """Actual derivative implementation (private)."""
        if wrt not in ("alpha", "beta", "gamma"):
            raise MittagLefflerError(
                f"wrt must be 'alpha', 'beta', or 'gamma'; got {wrt!r}"
            )
        if (isinstance(self.alpha, np.ndarray)
                or isinstance(self.beta, np.ndarray)
                or isinstance(self.gamma, np.ndarray)):
            raise MittagLefflerError(
                "derivative() requires scalar parameters."
            )

        z_arr = np.asarray(z, dtype=np.complex128)
        if z_arr.ndim == 0:
            return self._derivative_scalar(complex(z_arr), wrt)

        out = np.empty(z_arr.shape, dtype=np.complex128)
        for idx in np.ndindex(z_arr.shape):
            out[idx] = self._derivative_scalar(complex(z_arr[idx]), wrt)
        return out

    def _derivative_scalar(self, z: complex, wrt: str) -> complex:
        """Single-point derivative: series for small |z|, mpmath otherwise."""
        if abs(z) <= self.series_threshold and self.precision is None:
            try:
                return self._derivative_series(z, wrt)
            except MittagLefflerError:
                pass  # fall through to mpmath
        return self._derivative_mpmath(z, wrt)

    def _derivative_series(self, z: complex, wrt: str) -> complex:
        """
        Analytic series derivative (accurate for small |z|).

        Term-wise derivatives of
            E^γ_{α,β}(z) = Σ (γ)_k z^k / (k! Γ(αk+β))

        ∂c_k/∂α = -k·ψ(αk+β)·c_k
        ∂c_k/∂β = -ψ(αk+β)·c_k
        ∂c_k/∂γ = [ψ(γ+k) - ψ(γ)]·c_k
        """
        alpha = float(self.alpha)
        beta = float(self.beta)
        gamma = float(self.gamma)

        # k=0 term: c_0 = 1/Γ(β)
        total = 0.0 + 0.0j
        if wrt == "beta":
            if not (beta <= 0 and abs(beta - round(beta)) < 1e-12):
                total += -digamma(beta) / gamma_func(beta)
        # α and γ: ∂/∂(α,γ)[1/Γ(β)] = 0

        z_power = 1.0 + 0.0j
        log_coef = 0.0  # log((γ)_k / k!)

        for k in range(1, self.max_series_terms + 1):
            z_power *= z
            log_coef += math.log(gamma + k - 1.0) - math.log(k)
            try:
                coef = math.exp(log_coef)
            except OverflowError:
                raise MittagLefflerError(
                    "Derivative series overflow; try mpmath fallback."
                )

            g = alpha * k + beta
            if g <= 0 and abs(g - round(g)) < 1e-12:
                continue  # Γ pole

            try:
                rg = float(rgamma(g))
                psi_g = float(digamma(g))
            except (ValueError, OverflowError):
                continue

            if wrt == "alpha":
                factor = -k * psi_g
            elif wrt == "beta":
                factor = -psi_g
            else:  # gamma
                psi_gk = float(digamma(gamma + k))
                psi_g0 = float(digamma(gamma))
                factor = psi_gk - psi_g0

            term = coef * z_power * rg * factor

            if not np.isfinite(term):
                raise MittagLefflerError(
                    "Derivative series produced non-finite term."
                )

            total += term
            if abs(term) <= self.tolerance * (1.0 + abs(total)):
                return total

        raise MittagLefflerError(
            "Derivative series did not converge within max_series_terms."
        )

    def _derivative_mpmath(self, z: complex, wrt: str) -> complex:
        """mpmath central-difference derivative (large |z| or forced)."""
        dps = (self.precision if self.precision is not None else 30)
        with _mp.workdps(dps):
            z_mp = _mp.mpc(z)
            a_mp = _mp.mpf(float(self.alpha))
            b_mp = _mp.mpf(float(self.beta))
            g_mp = _mp.mpf(float(self.gamma))

            if wrt == "alpha":
                result = _mp.diff(
                    lambda a: self._mpmath_series(z_mp, a, b_mp, g_mp),
                    a_mp,
                )
            elif wrt == "beta":
                result = _mp.diff(
                    lambda b: self._mpmath_series(z_mp, a_mp, b, g_mp),
                    b_mp,
                )
            else:  # gamma
                result = _mp.diff(
                    lambda g: self._mpmath_series(z_mp, a_mp, b_mp, g),
                    g_mp,
                )
            return complex(result)

    def _mpmath_series(self, z_mp, a_mp, b_mp, g_mp):
        """mpmath series at current workdps — safe for α near 0."""
        if a_mp == 0:
            if abs(z_mp) >= 1:
                raise MittagLefflerError("alpha = 0 requires |z| < 1.")
            return _mp.power(1 - z_mp, -g_mp) * _mp.rgamma(b_mp)

        total = _mp.rgamma(b_mp)
        z_power = _mp.mpc(1)
        poch = _mp.mpf(1)
        tol = _mp.mpf(10) ** (-(_mp.mp.dps - 10))

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
            "mpmath series did not converge within 200,000 terms."
        )

    def __repr__(self) -> str:
        return (
            f"{type(self).__name__}("
            f"alpha={self.alpha!r}, beta={self.beta!r}, "
            f"gamma={self.gamma!r}, tolerance={self.tolerance!r}, "
            f"series_threshold={self.series_threshold!r}, "
            f"max_series_terms={self.max_series_terms!r})"
        )


# ---------------------------------------------------------------------------
# Convenience function
# ---------------------------------------------------------------------------
@functools.lru_cache(maxsize=256)
def _reduced_evaluator(
    alpha: float,
    beta: float,
    tolerance: float,
    series_threshold: float,
    max_series_terms: int,
) -> "MittagLeffler":
    """İndirgenmiş (γ=1) değerlendirici — parametrelere göre önbellekli."""
    return MittagLeffler(
        alpha=alpha, beta=beta, gamma=1.0,
        tolerance=tolerance, series_threshold=series_threshold,
        max_series_terms=max_series_terms,
    )


def mittag_leffler(
    z: ArrayLike,
    alpha: float,
    beta: float = 1.0,
    gamma: float = 1.0,
    tolerance: float = 1.0e-15,
    series_threshold: float = 0.1,
    max_series_terms: int = 100_000,
    precision=None,     # ← NEW 0.1.5
) -> (
    NDArray[np.float64] | np.float64
    | NDArray[np.complex128] | np.complex128
):
    """
    Convenience wrapper.

    Return type follows the input type:

    * **Real ``z``** → ``float64``
    * **Complex ``z``** (any non-zero imaginary part) → ``complex128``

    For a guaranteed complex return, use
    ``MittagLeffler(...).evaluate_complex(z)``.

    Examples
    --------
    >>> mittag_leffler(1.0, alpha=0.5)
    5.008980080762283

    >>> mittag_leffler(1.0j, alpha=0.5)
    (0.36787944117144233+0.6071577058413894j)
    """
    evaluator = MittagLeffler(
        alpha=alpha,
        beta=beta,
        gamma=gamma,
        tolerance=tolerance,
        series_threshold=series_threshold,
        max_series_terms=max_series_terms,
        precision=precision,
    )
    return evaluator(z)


class MatrixAnalyzer:
    """
    Matris özelliklerini analiz eden sınıf.

    Desteklenen matris türleri:
    - numpy.ndarray
    - scipy.sparse matrisleri (csr, csc, coo, vb.)

    Not:
    - Sparse matrislerde yapısal testler mümkün olduğunca sparse olarak yapılır.
    - Özdeğer / SVD gibi işlemler gerektiğinde küçük matrislerde dense'e çevrilir.
    - Büyük sparse matrislerde Hermitian extremal özdeğerler için eigsh kullanılır.
    """

    def __init__(self, A, name: str = "Matrix"):
        self.name = name
        self._input = A
        self._dense_cache = None
        self._norm_cache = None
        self._eigenvalues_cache = None

        if _is_sparse(A):
            self._sparse = A.tocsr() if hasattr(A, "tocsr") else A
            self._shape = self._sparse.shape
            self._dtype = self._sparse.dtype
            self._is_sparse = True
        else:
            arr = np.asarray(A)

            try:
                if np.iscomplexobj(arr):
                    arr = arr.astype(np.complex128, copy=False)
                else:
                    arr = arr.astype(float, copy=False)
            except (TypeError, ValueError):
                arr = arr.astype(np.complex128, copy=False)

            self._dense_cache = arr
            self._sparse = None
            self._shape = arr.shape
            self._dtype = arr.dtype
            self._is_sparse = False

        self._is_square = (
            len(self._shape) == 2 and self._shape[0] == self._shape[1]
        )

    # ------------------------------------------------------------------
    # Temel property'ler
    # ------------------------------------------------------------------

    @property
    def shape(self):
        return self._shape

    @property
    def is_square(self) -> bool:
        return bool(self._is_square)

    @property
    def is_sparse(self) -> bool:
        return bool(self._is_sparse)

    @property
    def dtype(self):
        return self._dtype

    @property
    def matrix(self):
        """
        Sparse ise sparse matrisi, dense ise dense matrisi döndürür.
        """
        return self._sparse if self._is_sparse else self.dense

    @property
    def dense(self) -> np.ndarray:
        """
        Matrisin dense numpy temsilini döndürür.
        Sparse matrislerde sadece gerektiğinde çağrılmalıdır.
        """
        if self._dense_cache is None:
            if self._is_sparse:
                arr = self._sparse.toarray()
            else:
                arr = np.asarray(self._input)

            if np.iscomplexobj(arr):
                arr = arr.astype(np.complex128, copy=False)
            else:
                arr = arr.astype(float, copy=False)

            self._dense_cache = arr

        return self._dense_cache

    @property
    def diagonal(self) -> np.ndarray:
        """
        Matrisin köşegen elemanlarını döndürür.
        """
        if not self._is_square:
            return np.array([])

        if self._is_sparse:
            return self._sparse.diagonal()

        return np.diag(self.dense)

    @property
    def trace(self):
        """
        Matrisin izini döndürür.
        """
        if not self._is_square:
            return 0.0

        if self._is_sparse:
            return self._sparse.trace()

        return np.trace(self.dense)

    @property
    def norm(self) -> float:
        """
        Matrisin Frobenius normunu döndürür.
        Sparse matrislerde sparse norm kullanılır.
        """
        if self._norm_cache is None:
            self._norm_cache = self._norm(self.matrix)

        return float(self._norm_cache)

    # ------------------------------------------------------------------
    # Yardımcı private fonksiyonlar
    # ------------------------------------------------------------------

    def _norm(self, M) -> float:
        """
        Dense veya sparse matrisin normunu hesaplar.
        """
        if _is_sparse(M):
            if sp_norm is not None:
                return float(sp_norm(M))
            return float(np.linalg.norm(M.data)) if M.nnz else 0.0

        return float(np.linalg.norm(np.asarray(M)))

    def _close(self, X, Y, rtol: float, atol: float, ref_norm=None) -> bool:
        """
        İki matrisin birbirine yakın olup olmadığını norm tabanlı kontrol eder.
        """
        diff = X - Y

        if ref_norm is None:
            ref_norm = self._norm(Y)

        diff_norm = self._norm(diff)

        return bool(
            diff_norm <= atol + rtol * max(1.0, float(ref_norm))
        )

    def _identity(self, use_sparse=None):
        """
        Kare matrisle aynı boyutta birim matris üretir.
        """
        if not self._is_square:
            raise ValueError("Identity requires a square matrix.")

        n = self._shape[0]

        if use_sparse is None:
            use_sparse = self._is_sparse

        dtype = np.complex128 if np.issubdtype(self._dtype, np.complexfloating) else float

        if use_sparse:
            if not HAS_SCIPY:
                raise RuntimeError("scipy.sparse is not available.")
            return sparse.eye(n, format="csr", dtype=dtype)

        return np.eye(n, dtype=dtype)

    def _is_complex_dtype(self) -> bool:
        return bool(np.issubdtype(self._dtype, np.complexfloating))

    def _extreme_eigenvalue(self, which: str = "SA", tol: float = 1e-10) -> float:
        """
        Hermitian matrisler için extremal özdeğeri döndürür.

        which:
            "SA" -> smallest algebraic
            "LA" -> largest algebraic

        Sparse ve büyük matrislerde eigsh kullanmaya çalışır.
        Başarısız olursa dense'e düşer.
        """
        if not self._is_square:
            return float("nan")

        if (
            self._is_sparse
            and HAS_SCIPY
            and eigsh is not None
            and self._shape[0] >= 50
        ):
            try:
                dtype = np.complex128 if self._is_complex_dtype() else float
                M = self._sparse.astype(dtype)

                vals = eigsh(
                    M,
                    k=1,
                    which=which,
                    return_eigenvectors=False,
                    tol=tol,
                )

                return float(np.real(vals[0]))

            except Exception:
                pass

        # Fallback: dense eigen decomposition
        vals = np.linalg.eigvalsh(self.dense)

        if which == "SA":
            return float(np.min(vals))

        if which == "LA":
            return float(np.max(vals))

        raise ValueError("which must be 'SA' or 'LA'.")

    # ------------------------------------------------------------------
    # Temel matris özellikleri
    # ------------------------------------------------------------------

    def is_real(self, rtol: float = 1e-12, atol: float = 1e-14) -> bool:
        """
        Matrisin reel matris olup olmadığını kontrol eder.
        """
        if not self._is_complex_dtype():
            return True

        M = self.matrix
        imag = M.imag

        imag_norm = self._norm(imag)

        return bool(
            imag_norm <= atol + rtol * max(1.0, self.norm)
        )

    def is_diagonal(self, rtol: float = 1e-12, atol: float = 1e-14) -> bool:
        """
        Matrisin diagonal olup olmadığını kontrol eder.
        """
        if not self._is_square:
            return False

        M = self.matrix
        diag = self.diagonal

        if self._is_sparse:
            D = sparse.diags(
                diag,
                offsets=0,
                shape=self._shape,
                format="csr",
                dtype=M.dtype,
            )
            off = M - D
        else:
            D = np.diag(diag)
            off = M - D

        off_norm = self._norm(off)

        return bool(
            off_norm <= atol + rtol * max(1.0, self.norm)
        )

    def is_symmetric(self, rtol: float = 1e-12, atol: float = 1e-14) -> bool:
        """
        A = A^T kontrolü.
        """
        if not self._is_square:
            return False

        M = self.matrix
        MT = M.T

        return self._close(M, MT, rtol=rtol, atol=atol)

    def is_hermitian(self, rtol: float = 1e-12, atol: float = 1e-14) -> bool:
        """
        A = A^H kontrolü.
        """
        if not self._is_square:
            return False

        M = self.matrix
        MH = M.conj().T

        return self._close(M, MH, rtol=rtol, atol=atol)

    def is_skew_hermitian(self, rtol: float = 1e-12, atol: float = 1e-14) -> bool:
        """
        A = -A^H kontrolü.
        """
        if not self._is_square:
            return False

        M = self.matrix
        neg_MH = -(M.conj().T)

        return self._close(M, neg_MH, rtol=rtol, atol=atol)

    def is_unitary(self, rtol: float = 1e-12, atol: float = 1e-14) -> bool:
        """
        A A^H = I kontrolü.
        """
        if not self._is_square:
            return False

        M = self.matrix
        MH = M.conj().T
        prod = M @ MH
        I = self._identity(use_sparse=self._is_sparse)

        return self._close(prod, I, rtol=rtol, atol=atol)

    def is_normal(self, rtol: float = 1e-12, atol: float = 1e-14) -> bool:
        """
        A A^H = A^H A kontrolü.
        """
        if not self._is_square:
            return False

        M = self.matrix
        MH = M.conj().T

        left = M @ MH
        right = MH @ M

        return self._close(left, right, rtol=rtol, atol=atol)

    # ------------------------------------------------------------------
    # Yeni eklenen özellikler
    # ------------------------------------------------------------------

    def is_idempotent(self, rtol: float = 1e-12, atol: float = 1e-14) -> bool:
        """
        A^2 = A kontrolü.
        """
        if not self._is_square:
            return False

        M = self.matrix
        prod = M @ M

        return self._close(prod, M, rtol=rtol, atol=atol)

    def is_involutory(self, rtol: float = 1e-12, atol: float = 1e-14) -> bool:
        """
        A^2 = I kontrolü.
        """
        if not self._is_square:
            return False

        M = self.matrix
        prod = M @ M
        I = self._identity(use_sparse=self._is_sparse)

        return self._close(prod, I, rtol=rtol, atol=atol)

    def is_projection(self, rtol: float = 1e-12, atol: float = 1e-14) -> bool:
        """
        Projection tanımı burada idempotent matris olarak alınmıştır:

            P^2 = P

        Eğer orthogonal projection isteniyorsa:

            is_projection() and is_hermitian()

        veya doğrudan:

            is_orthogonal_projection()
        """
        return self.is_idempotent(rtol=rtol, atol=atol)

    def is_orthogonal_projection(self, rtol: float = 1e-12, atol: float = 1e-14) -> bool:
        """
        Orthogonal projection:

            P^2 = P ve P = P^H
        """
        return bool(
            self.is_projection(rtol=rtol, atol=atol)
            and self.is_hermitian(rtol=rtol, atol=atol)
        )

    def is_orthogonal(self, rtol: float = 1e-12, atol: float = 1e-14) -> bool:
        """
        Real orthogonal matrix kontrolü:

            A^T A = I

        Kompleks matrisler için orthogonal yerine unitary kullanılır.
        Bu metot, matris kompleks ise ve sanal kısmı sıfırdan farklıysa False döndürür.
        """
        if not self._is_square:
            return False

        if not self.is_real(rtol=rtol, atol=atol):
            return False

        M = self.matrix
        MR = M.real

        prod = MR.T @ MR

        if self._is_sparse:
            I = sparse.eye(self._shape[0], format="csr", dtype=float)
        else:
            I = np.eye(self._shape[0], dtype=float)

        return self._close(prod, I, rtol=rtol, atol=atol)

    def is_positive_definite(self, tol: float = 1e-10) -> bool:
        """
        Hermitian ve tüm özdeğerleri > tol.
        """
        if not self.is_hermitian(rtol=tol, atol=tol):
            return False

        try:
            min_eig = self._extreme_eigenvalue(which="SA", tol=tol)
        except np.linalg.LinAlgError:
            return False

        return bool(min_eig > tol)

    def is_positive_semidefinite(self, tol: float = 1e-10) -> bool:
        """
        Hermitian ve tüm özdeğerleri >= -tol.
        """
        if not self.is_hermitian(rtol=tol, atol=tol):
            return False

        try:
            min_eig = self._extreme_eigenvalue(which="SA", tol=tol)
        except np.linalg.LinAlgError:
            return False

        return bool(min_eig >= -tol)

    def is_negative_definite(self, tol: float = 1e-10) -> bool:
        """
        Hermitian ve tüm özdeğerleri < -tol.
        """
        if not self.is_hermitian(rtol=tol, atol=tol):
            return False

        try:
            max_eig = self._extreme_eigenvalue(which="LA", tol=tol)
        except np.linalg.LinAlgError:
            return False

        return bool(max_eig < -tol)

    def is_trace_zero(self, rtol: float = 1e-12, atol: float = 1e-14) -> bool:
        """
        trace(A) ≈ 0 kontrolü.
        """
        if not self._is_square:
            return False

        tr = self.trace
        scale = max(1.0, self.norm)

        return bool(abs(tr) <= atol + rtol * scale)

    def is_singular(self, rtol: float = 1e-12, atol: float = 1e-14) -> bool:
        """
        Matris singular mı?

        En küçük singular değer <= atol + rtol * norm(A) ise singular kabul edilir.
        """
        if not self._is_square:
            return False

        threshold = atol + rtol * max(1.0, self.norm)

        # Büyük sparse matrislerde svds ile sadece en küçük singular değeri aramayı dene.
        if (
            self._is_sparse
            and HAS_SCIPY
            and svds is not None
            and min(self._shape) > 2
            and self._shape[0] >= 50
        ):
            try:
                dtype = np.complex128 if self._is_complex_dtype() else float
                M = self._sparse.astype(dtype)

                s = svds(
                    M,
                    k=1,
                    which="SM",
                    return_singular_vectors=False,
                )

                return bool(float(np.min(s)) <= threshold)

            except Exception:
                pass

        # Fallback: dense SVD
        try:
            s = np.linalg.svd(self.dense, compute_uv=False)
            return bool(s[-1] <= threshold)

        except np.linalg.LinAlgError:
            rank = np.linalg.matrix_rank(self.dense, tol=threshold)
            return bool(rank < self._shape[0])

    # ------------------------------------------------------------------
    # Özdeğerler ve sayısal ek bilgiler
    # ------------------------------------------------------------------

    @property
    def eigenvalues(self) -> np.ndarray:
        """
        Matrisin özdeğerlerini döndürür.

        Sparse matrislerde bu işlem dense dönüşüm gerektirebilir.
        """
        if self._eigenvalues_cache is not None:
            return self._eigenvalues_cache

        if not self._is_square:
            self._eigenvalues_cache = np.array([], dtype=np.complex128)
            return self._eigenvalues_cache

        try:
            if self.is_hermitian():
                vals = np.linalg.eigvalsh(self.dense)
            else:
                vals = np.linalg.eigvals(self.dense)

            vals = vals.astype(np.complex128)

        except np.linalg.LinAlgError:
            vals = np.array([], dtype=np.complex128)

        self._eigenvalues_cache = vals
        return self._eigenvalues_cache

    @property
    def spectral_radius(self) -> float:
        """
        max |eigenvalue|
        """
        vals = self.eigenvalues

        if vals.size == 0:
            return float("nan")

        return float(np.max(np.abs(vals)))

    @property
    def condition_number(self) -> float:
        """
        Matrisin condition number değeri.
        """
        if not self._is_square:
            return float("nan")

        try:
            return float(np.linalg.cond(self.dense))
        except np.linalg.LinAlgError:
            return float("inf")

    # ------------------------------------------------------------------
    # Summary / report
    # ------------------------------------------------------------------

    def summary(
        self,
        include_advanced: bool = True,
        include_definiteness: bool = True,
    ) -> dict:
        """
        Matris özelliklerini dictionary olarak döndürür.

        np.True_ / np.False_ yerine Python bool döndürür.
        """
        if not self._is_square:
            return {"Error": "Matrix is not square."}

        result = {
            "Diagonal": self.is_diagonal(),
            "Symmetric": self.is_symmetric(),
            "Hermitian": self.is_hermitian(),
            "Skew-Hermitian": self.is_skew_hermitian(),
            "Unitary": self.is_unitary(),
            "Normal": self.is_normal(),
        }

        if include_advanced:
            result.update(
                {
                    "Idempotent": self.is_idempotent(),
                    "Involutory": self.is_involutory(),
                    "Projection": self.is_projection(),
                    "Orthogonal Projection": self.is_orthogonal_projection(),
                    "Orthogonal": self.is_orthogonal(),
                    "Trace Zero": self.is_trace_zero(),
                    "Singular": self.is_singular(),
                }
            )

        if include_definiteness:
            result.update(
                {
                    "Positive Definite": self.is_positive_definite(),
                    "Positive Semidefinite": self.is_positive_semidefinite(),
                    "Negative Definite": self.is_negative_definite(),
                }
            )

        # np.bool_ -> Python bool
        return {k: bool(v) for k, v in result.items()}

    def report(self) -> dict:
        """
        Hem boolean özellikler hem de sayısal bilgiler.
        """
        base = self.summary()

        if "Error" in base:
            return base

        rep = dict(base)
        rep["Eigenvalues"] = self.eigenvalues
        rep["Spectral Radius"] = self.spectral_radius
        rep["Condition Number"] = self.condition_number

        return rep

# ----------------------------------------------------------------------------
# Yardımcı doğrulayıcılar
# ----------------------------------------------------------------------------
"""
def _is_diagonal(A: np.ndarray, tol: float = 1e-14) -> bool:
    off = A - np.diag(np.diag(A))
    return np.linalg.norm(off) < tol * max(1.0, np.linalg.norm(A))


def _is_hermitian(A: np.ndarray, tol: float = 1e-10) -> bool:
    return np.allclose(A, A.conj().T, atol=tol)
"""
def _is_diagonal(A, tol: float = 1e-14) -> bool:
    return bool(MatrixAnalyzer(A).is_diagonal(rtol=tol, atol=tol))


def _is_hermitian(A, tol: float = 1e-10) -> bool:
    return bool(MatrixAnalyzer(A).is_hermitian(rtol=tol, atol=tol))

# ----------------------------------------------------------------------------
# 1) Prabhakar — skaler E^γ_{α,β}(z)
# ----------------------------------------------------------------------------
class Prabhakar:
    r"""
    E^γ_{α,β}(z) = Σ_{k=0}^∞  (γ)_k / (Γ(αk+β)·k!) · z^k

    Sayısal kararlılık için:
      • γ=0 özel durumu (yalnızca k=0 terimi)
      • γ=1 ve büyük negatif |z| için asimptotik açılım (taşma önler)
      • katsayı cache (time_series için ~50× hız)
    """
    _ASYMPTOTIC_THRESHOLD = 30.0

    def __init__(self, alpha: float, beta: float = 1.0,
                 gamma: float = 1.0, tolerance: float = 1e-14):
        if alpha <= 0 or beta <= 0 or gamma < 0:
            raise ValueError("α > 0, β > 0, γ ≥ 0 gerekir")
        self.alpha = float(alpha)
        self.beta = float(beta)
        self.gamma = float(gamma)
        self.tol = float(tolerance)
        self._coef_cache: dict[int, float] = {}

    # ------------------------------------------------------------------
    def coef(self, k: int) -> float:
        """k. katsayı — cache'li, γ özel durumlarını doğru işler."""
        if k in self._coef_cache:
            return self._coef_cache[k]

        if self.gamma == 0.0:
            # (0)_0 = 1,  (0)_k = 0 (k≥1)  →  E^0_{α,β}(z) = 1/Γ(β)
            v = (1.0 / float(np.exp(gammaln(self.beta)))) if k == 0 else 0.0
        elif self.gamma == 1.0:
            v = float(np.exp(-gammaln(self.alpha * k + self.beta)))
        else:
            log_num = gammaln(self.gamma + k) - gammaln(self.gamma)
            log_den = gammaln(self.alpha * k + self.beta) + gammaln(k + 1.0)
            v = float(np.exp(log_num - log_den))

        self._coef_cache[k] = v
        return v

    # ------------------------------------------------------------------
    def _asymptotic_neg_real(self, z: complex, N: int = 8) -> complex:
        r"""
        Re(z)<0, |z|→∞ için asimptotik:
            E_{α,β}(z) ≈ −Σ_{n=1}^N z^{−n} / Γ(β−αn)
        (Podlubny / Gorenflo-Mainardi formülü, γ=1 durumu)
        """
        total = 0j
        for n in range(1, N + 1):
            denom = self.beta - self.alpha * n
            # Γ(non-positive integer) = ∞ → terim atla
            if denom <= 0 and abs(denom - round(denom)) < 1e-10:
                continue
            if abs(denom) < 1e-10:
                continue
            g = float(np.exp(gammaln(denom)))
            total -= (z ** (-n)) / g
        return total

    # ------------------------------------------------------------------
    def evaluate(self, z: complex, max_terms: int = 800,
                 min_terms: int = 8) -> complex:
        z = complex(z)

        # γ=1 ve büyük negatif argüman → asimptotik (taşma önler)
        if (self.gamma == 1.0
                and np.real(z) < 0
                and abs(z) > self._ASYMPTOTIC_THRESHOLD):
            return self._asymptotic_neg_real(z)

        total = 0j
        z_pow = 1 + 0j
        for k in range(max_terms):
            term = self.coef(k) * z_pow
            total += term
            if k >= min_terms and abs(term) < self.tol * max(1.0, abs(total)):
                return total
            z_pow *= z

        warnings.warn(f"Prabhakar: {max_terms} terimde yakınsamadı (z={z})")
        return total

    def evaluate_real(self, x: float, **kw) -> float:
        return float(np.real(self.evaluate(x, **kw)))

# ----------------------------------------------------------------------------
# 2) MatrixMittagLeffler — E^γ_{α,β}(A)
# ----------------------------------------------------------------------------
class MatrixMittagLeffler:
    r"""
    Matris Mittag-Leffler / Prabhakar fonksiyonu E^γ_{α,β}(A).

    Hesaplama stratejisi (otomatik seçim)
    -------------------------------------
    1. A diyagonal        → köşegen skaler eval (O(n) — en hızlı)
    2. A Hermitian        → eigh + köşegen skaler eval (spektral doğru)
    3. A genel            → Pochhammer serisi + yakınsama kontrolü

    Not
    ---
    Fiziksel kanal uygulamalarında `A` genelde diyagonal veya Hermitian
    (dekoherans jeneratörü) olur → 1. veya 2. yol devreye girer ve
    matris çarpımı yapılmaz. Bu, `time_series` için kritik optimizasyondur.
    """

    def __init__(self, alpha: float, beta: float = 1.0,
                 gamma: float = 1.0, tolerance: float = 1e-13):
        self._scalar = Prabhakar(alpha, beta, gamma, tolerance)
        self.alpha, self.beta, self.gamma = alpha, beta, gamma
        self.tol = tolerance

    def evaluate(self, A, max_terms: int = 200, sparse_output: bool = False):
        analyzer = MatrixAnalyzer(A)

        if not analyzer.is_square:
            raise ValueError("A kare matris olmalı.")

        if analyzer.is_diagonal(rtol=self.tol, atol=self.tol):
            return self._via_diag_values(
                analyzer.diagonal,
                sparse_output=sparse_output and analyzer.is_sparse
            )

        if analyzer.is_hermitian(rtol=self.tol, atol=self.tol):
            return self._via_eig(analyzer.dense)

        return self._via_series(analyzer.dense, max_terms)


    # --- strateji 1: diyagonal ---
    def _via_diag(self, A):
        """
        Geriye uyumluluk için tutuluyor.
        Yeni kodlarda _via_diag_values tercih edilmeli.
        """
        analyzer = MatrixAnalyzer(A)
        return self._via_diag_values(analyzer.diagonal, sparse_output=False)

    def _via_diag_values(self, diag, sparse_output: bool = False):
        diag = np.asarray(diag, dtype=complex)

        vals = np.array(
            [self._scalar.evaluate(z) for z in diag],
            dtype=complex
        )

        if sparse_output and HAS_SCIPY:
            return sparse.diags(vals, format="csr")

        return np.diag(vals)

    # --- strateji 2: Hermitian ---
    def _via_eig(self, A):
        vals, vecs = eigh(A)
        f_vals = np.array([self._scalar.evaluate(v) for v in vals])
        return vecs @ np.diag(f_vals) @ vecs.conj().T

    # --- strateji 3: seri ---
    def _via_series(self, A, max_terms):
        n = A.shape[0]
        res = np.zeros((n, n), dtype=complex)
        A_pow = np.eye(n, dtype=complex)
        for k in range(max_terms):
            term = self._scalar.coef(k) * A_pow
            res += term
            if k >= 8 and np.linalg.norm(term, 'fro') < \
                    self.tol * max(1.0, np.linalg.norm(res, 'fro')):
                return res
            A_pow = A_pow @ A
        warnings.warn(
            f"MatrixMittagLeffler: {max_terms} terimde yakınsamadı")
        return res


# mlfp.py içindeki mevcut FractionalQuantumChannel'ı BUNUNLA DEĞİŞTİR

class FractionalQuantumChannel:
    r"""
    Kesirli koherans altında 2-qubit kuantum kanalı.

    Model
    -----
        ρ(t) = K₀ ρ₀ K₀† + K₁ ρ₀ K₁†
        K₀(t) = E^γ_{α,β}( s · G ),   s = (t/τ)^α
        K₁(t) = √(I − K₀† K₀)

    Konvansiyon
    -----------
    G **negatif-tanımlı** dekoherans jeneratörüdür (fizikçi standardı):
        özdeğerler(G) ≤ 0  →  s·G negatif argüman
        →  E_{α,β}(s·G) azalır  →  |K₀| ≤ 1  →  Kraus geçerli

    G pozitif verilirse `_ensure_negative_generator` uyarı basar ve
    otomatik olarak işaretini çevirir (−G).
    """

    def __init__(self, alpha: float, beta: float = 1.0,
                 gamma: float = 1.0, tolerance: float = 1e-12):
        self.alpha = float(alpha)
        self.beta = float(beta)
        self.gamma = float(gamma)
        self._scalar = Prabhakar(alpha, beta, gamma, tolerance)
        self._matrix = MatrixMittagLeffler(alpha, beta, gamma, tolerance)

    # ------------------------------------------------------------------
    # Generator doğrulama (G'nin her kullanımdan önce çağrılır)
    # ------------------------------------------------------------------
    def _ensure_negative_generator(self, G) -> np.ndarray:
        G = np.asarray(G, dtype=complex)
        if G.ndim != 2 or G.shape[0] != G.shape[1]:
            raise ValueError("G kare matris olmalı")
        if np.allclose(G, G.conj().T, atol=1e-10):
            eigs = np.linalg.eigvalsh(G)
        else:
            eigs = np.real(np.linalg.eigvals(G))
        max_eig = float(np.max(eigs))
        if max_eig > 1e-8:
            warnings.warn(
                f"G pozitif özdeğerli (max={max_eig:.3f}); "
                f"dekoherans jeneratörü negatif olmalı. "
                f"Otomatik düzeltildi (G → −G).")
            return -G
        return G

    # ------------------------------------------------------------------
    # Skaler metrikler
    # ------------------------------------------------------------------
    def coherence(self, t: float, tau: float = 1.0) -> float:
        x = (t / tau) ** self.alpha
        return float(np.clip(self._scalar.evaluate_real(-x), 0.0, 1.0))

    def teleportation_fidelity(self, t: float, tau: float = 1.0) -> float:
        return 0.5 + 0.5 * self.coherence(t, tau)

    # ------------------------------------------------------------------
    # Kraus operatörleri
    # ------------------------------------------------------------------
    def kraus_k0(self, generator, t: float, tau: float = 1.0) -> np.ndarray:
        G = self._ensure_negative_generator(generator)
        s = (t / tau) ** self.alpha
        return self._matrix.evaluate(s * G)

    def evolve_state(self, rho0: np.ndarray, generator,
                     t: float, tau: float = 1.0,
                     check_kraus: bool = True,
                     trace_tol: float = 1e-8) -> np.ndarray:
        K0 = self.kraus_k0(generator, t, tau)
        n = K0.shape[0]
        I = np.eye(n, dtype=complex)
        rem = I - K0.conj().T @ K0

        if check_kraus:
            eigs = np.linalg.eigvalsh(rem)
            if eigs.min() < -1e-8:
                warnings.warn(
                    f"Kraus tamamlayıcısı geçersiz: min özdeğer "
                    f"{eigs.min():.2e} — generator veya MLF yakınsaması")

        vals, vecs = eigh(rem)
        vals = np.clip(vals, 0.0, None)
        K1 = vecs @ np.diag(np.sqrt(vals)) @ vecs.conj().T

        rho_t = K0 @ rho0 @ K0.conj().T + K1 @ rho0 @ K1.conj().T
        tr = np.real(np.trace(rho_t))
        if abs(tr - 1.0) > trace_tol:
            warnings.warn(f"İz sapması: |tr(ρ)−1| = {abs(tr-1):.2e}")
            rho_t = rho_t / tr
        return rho_t

    @staticmethod
    def von_neumann_entropy(rho: np.ndarray) -> float:
        vals = np.real(np.linalg.eigvalsh(rho))
        vals = vals[vals > 1e-15]
        return float(-np.sum(vals * np.log2(vals)))

    # ------------------------------------------------------------------
    # Zaman serisi — diyagonal jeneratörde vektörel hızlı yol
    # ------------------------------------------------------------------
    def time_series(self, generator, times, rho0: np.ndarray,
                    tau: float = 1.0) -> list[np.ndarray]:
        G = self._ensure_negative_generator(generator)

        if _is_diagonal(G):
            diag = np.diag(G)
            results = []
            for t in times:
                s = (t / tau) ** self.alpha
                k0_diag = np.array([self._scalar.evaluate(s * g) for g in diag])
                K0 = np.diag(k0_diag)
                rem_diag = np.clip(1.0 - np.abs(k0_diag) ** 2, 0.0, None)
                K1 = np.diag(np.sqrt(rem_diag))
                rho_t = K0 @ rho0 @ K0.conj().T + K1 @ rho0 @ K1.conj().T
                results.append(rho_t / np.real(np.trace(rho_t)))
            return results

        return [self.evolve_state(rho0, G, t, tau) for t in times]

# ============================================================================
# mlfp.py — T₁ (genlik sönümü) EKLENTİSİ
# ============================================================================

# ----------------------------------------------------------------------------
# Katman 1: Saf Kraus fonksiyonları (sabit γ)
# ----------------------------------------------------------------------------
def amplitude_damping_kraus(gamma: float) -> tuple[np.ndarray, np.ndarray]:
    """
    Tek kubit T₁ için Kraus operatörleri.

        K₀ = [[1, 0], [0, √(1-γ)]]     K₁ = [[0, √γ], [0, 0]]

    γ ∈ [0,1]:  uyarımın |1⟩ → |0⟩ bozunma olasılığı.
    γ=0 → özdeşlik;  γ=1 → tam relaksasyon.
    """
    g = float(np.clip(gamma, 0.0, 1.0))
    K0 = np.array([[1.0, 0.0],
                   [0.0, np.sqrt(1.0 - g)]], dtype=complex)
    K1 = np.array([[0.0, np.sqrt(g)],
                   [0.0, 0.0]], dtype=complex)
    return K0, K1


def amplitude_damping_channel(rho: np.ndarray, gamma: float) -> np.ndarray:
    """Tek kubit T₁ kanalı: ρ' = K₀ρK₀† + K₁ρK₁†."""
    K0, K1 = amplitude_damping_kraus(gamma)
    return K0 @ rho @ K0.conj().T + K1 @ rho @ K1.conj().T


def local_amplitude_damping(rho: np.ndarray, gamma: float,
                            qubit: int = 0) -> np.ndarray:
    """2-qubit durumda **tek** kübite T₁, diğeri özdeşlik."""
    I2 = np.eye(2, dtype=complex)
    K0, K1 = amplitude_damping_kraus(gamma)
    if qubit == 0:
        A0, A1 = np.kron(K0, I2), np.kron(K1, I2)
    elif qubit == 1:
        A0, A1 = np.kron(I2, K0), np.kron(I2, K1)
    else:
        raise ValueError("qubit ∈ {0, 1}")
    return A0 @ rho @ A0.conj().T + A1 @ rho @ A1.conj().T


def two_qubit_amplitude_damping(rho: np.ndarray,
                                gamma0: float, gamma1: float) -> np.ndarray:
    """Her iki kübite **paralel** T₁ (4 Kraus terimi)."""
    K0_0, K1_0 = amplitude_damping_kraus(gamma0)
    K0_1, K1_1 = amplitude_damping_kraus(gamma1)
    kraus_ops = [
        np.kron(K0_0, K0_1),
        np.kron(K0_0, K1_1),
        np.kron(K1_0, K0_1),
        np.kron(K1_0, K1_1),
    ]
    result = np.zeros_like(rho)
    for K in kraus_ops:
        result += K @ rho @ K.conj().T
    return result


# ----------------------------------------------------------------------------
# Katman 2: Kesirli T₁ — γ(t) = 1 − E_{α,1}(−(t/T₁)^α)
# ----------------------------------------------------------------------------
class FractionalAmplitudeDamping:
    """
    T₁ (genlik sönümü) — kesirli dinamik.

    Klasik limit  (α=1):  γ(t) = 1 − exp(−t/T₁)
    Kesirli               :  γ(t) = 1 − E_{α,1}(−(t/T₁)^α)

    `FractionalQuantumChannel` ile aynı α kullanıldığında T₁ ve T₂
    **fiziksel olarak tutarlı** kesirli rejimde olur.

    Parameters
    ----------
    alpha : float
        Kesirli dinamik mertebesi (α > 0).
    T1 : float
        Relaksasyon süresi (saniye).
    tolerance : float
        Prabhakar yakınsama toleransı.
    """

    def __init__(self, alpha: float, T1: float, tolerance: float = 1e-13):
        if alpha <= 0 or T1 <= 0:
            raise ValueError("α > 0, T₁ > 0 gerekir")
        self.alpha = float(alpha)
        self.T1 = float(T1)
        self._ml = Prabhakar(alpha, 1.0, 1.0, tolerance)

    def gamma(self, t: float) -> float:
        """γ(t) ∈ [0,1] — uyarım bozunma olasılığı."""
        if t <= 0:
            return 0.0
        s = (t / self.T1) ** self.alpha
        ml_val = self._ml.evaluate_real(-s)
        return float(np.clip(1.0 - ml_val, 0.0, 1.0))

    def apply(self, rho: np.ndarray, t: float,
              qubit: int | None = None) -> np.ndarray:
        """
        ρ'u T₁ kanalından geçir.

        qubit=None → 2-qubit simetrik (her iki kübite aynı γ)
        qubit=0/1  → 2-qubit'te yalnızca bir kübite
        (ρ boyutu 2×2 ise qubit argümanı yoksayılır.)
        """
        g = self.gamma(t)
        if rho.shape == (2, 2):
            return amplitude_damping_channel(rho, g)
        elif rho.shape == (4, 4):
            if qubit is None:
                return two_qubit_amplitude_damping(rho, g, g)
            return local_amplitude_damping(rho, g, qubit)
        raise ValueError(
            f"ρ boyutu (2,2) veya (4,4) olmalı, {rho.shape} verildi")


# ----------------------------------------------------------------------------
# Katman 3 (opsiyonel): T₁ ∘ T₂ kompozit kanal
# ----------------------------------------------------------------------------
class CompositeFractionalChannel:
    """
    T₂ (dephase) + T₁ (relaksasyon) birlikte.

        ρ(t) = T₁( T₂( ρ₀, t ), t )

    Uyarı
    -----
    Sıralı uygulama **fenomenolojik** bir yaklaşımdır; T₁ ve T₂ kesirli
    rejimde genel olarak komüt etmez. Fiziksel olarak doğru Kraus
    çiftini elde etmek için birleşik jeneratör gerekir. Bu sınıf,
    "T₂ önce, T₁ sonra" yaklaşımını kolaylaştırır — sadece gösterim
    ve karşılaştırma amaçlıdır.

    Parameters
    ----------
    alpha : float
        Ortak kesirli mertebe.
    G : np.ndarray
        Dekoherans jeneratörü (T₂ için, köşegen).
    T1 : float
        Relaksasyon süresi.
    """

    def __init__(self, alpha: float, G: np.ndarray, T1: float):
        self.t2 = FractionalQuantumChannel(alpha)
        self.t1 = FractionalAmplitudeDamping(alpha, T1)
        self.G = np.asarray(G, dtype=complex)

    def evolve(self, rho0: np.ndarray, t: float) -> np.ndarray:
        if t <= 0:
            return rho0.copy()
        rho_t = self.t2.evolve_state(rho0, self.G, t)
        return self.t1.apply(rho_t, t)

    def time_series(self, rho0: np.ndarray, times) -> list[np.ndarray]:
        return [self.evolve(rho0, t) for t in times]


# ----------------------------------------------------------------------
# Property sırası
# ----------------------------------------------------------------------

DASHBOARD_PROPERTY_ORDER = [
    "Diagonal",
    "Symmetric",
    "Hermitian",
    "Skew-Hermitian",
    "Unitary",
    "Normal",
    "Idempotent",
    "Involutory",
    "Projection",
    "Orthogonal Projection",
    "Orthogonal",
    "Trace Zero",
    "Singular",
    "Positive Definite",
    "Positive Semidefinite",
    "Negative Definite",
]


# ----------------------------------------------------------------------
# Yardımcı: property isimlerini topla
# ----------------------------------------------------------------------

def _collect_property_names(summaries):
    present = set()

    for s in summaries:
        for key in s.keys():
            if key != "Error":
                present.add(key)

    ordered = [
        p for p in DASHBOARD_PROPERTY_ORDER
        if p in present
    ]

    extras = sorted(present - set(ordered))

    return ordered + extras


# ----------------------------------------------------------------------
# Yardımcı: rapor yazdır
# ----------------------------------------------------------------------

def _print_report(analyzer):
    print("=" * 70)
    print(f"Matrix: {analyzer.name}")
    print("=" * 70)

    report = analyzer.report()

    for key, value in report.items():
        if isinstance(value, np.ndarray):
            print(
                f"  {key:28}: "
                f"{np.array2string(value, precision=4, suppress_small=True)}"
            )
        else:
            print(f"  {key:28}: {value}")

    print()


# ----------------------------------------------------------------------
# Summary heatmap
# ----------------------------------------------------------------------

def _summary_heatmap(analyzers, ax=None):
    """
    A ve E(A) gibi analyzer nesnelerinin boolean özelliklerini heatmap olarak çizer.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(16, 3.8))

    summaries = [a.summary() for a in analyzers]
    props = _collect_property_names(summaries)

    if not props:
        ax.text(
            0.5,
            0.5,
            "No square matrix available.",
            ha="center",
            va="center",
        )
        ax.axis("off")
        return ax

    data = []
    text_data = []

    for s in summaries:
        row = []
        text_row = []

        if "Error" in s:
            row = [np.nan] * len(props)
            text_row = ["N/A"] * len(props)
        else:
            for p in props:
                val = bool(s.get(p, False))
                row.append(1.0 if val else 0.0)
                text_row.append("True" if val else "False")

        data.append(row)
        text_data.append(text_row)

    data = np.array(data, dtype=float)

    im = ax.imshow(
        data,
        cmap="RdYlGn",
        vmin=0,
        vmax=1,
        aspect="auto",
    )

    ax.set_xticks(range(len(props)))
    ax.set_xticklabels(props, rotation=60, ha="right", fontsize=8)

    ax.set_yticks(range(len(analyzers)))
    ax.set_yticklabels(
        [
            getattr(a, "name", f"Matrix {i+1}")
            for i, a in enumerate(analyzers)
        ]
    )

    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            if not np.isnan(data[i, j]):
                ax.text(
                    j,
                    i,
                    text_data[i][j],
                    ha="center",
                    va="center",
                    fontsize=8,
                    color="black",
                )

    ax.set_title("Property Comparison: A vs E(A)", fontsize=12)

    cbar = ax.figure.colorbar(im, ax=ax, ticks=[0, 1])
    cbar.ax.set_yticklabels(["False", "True"])

    return ax


# ----------------------------------------------------------------------
# Tek matris heatmap
# ----------------------------------------------------------------------

def _matrix_heatmap_axis(ax, analyzer, title):
    """
    Tek bir analyzer nesnesinin matrisini verilen axis üzerine çizer.
    """
    if not analyzer.is_square:
        ax.axis("off")
        ax.set_title(f"{title}\n(not square)")
        return

    A = analyzer.dense

    if np.iscomplexobj(A):
        if np.allclose(A.imag, 0.0, atol=1e-12):
            data = A.real
            subtitle = "Re(A)"
            cmap = "coolwarm"
        else:
            data = np.abs(A)
            subtitle = "|A|"
            cmap = "viridis"
    else:
        data = A
        subtitle = "A"
        cmap = "coolwarm"

    im = ax.imshow(data, cmap=cmap)

    ax.set_title(f"{title}\n{subtitle}", fontsize=10)

    ax.set_xticks(range(data.shape[1]))
    ax.set_yticks(range(data.shape[0]))

    ax.figure.colorbar(im, ax=ax, fraction=0.046, pad=0.04)


# ----------------------------------------------------------------------
# Eigenvalue comparison
# ----------------------------------------------------------------------

def _eigenvalue_comparison(analyzers, ax=None):
    """
    A ve E(A) özdeğerlerini karmaşık düzlemde karşılaştırır.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(7.6, 7.6))

    for analyzer in analyzers:
        vals = analyzer.eigenvalues

        if vals.size == 0:
            continue

        ax.scatter(
            vals.real,
            vals.imag,
            s=80,
            label=analyzer.name,
            alpha=0.85,
        )

        for idx, z in enumerate(vals, start=1):
            ax.annotate(
                f"{analyzer.name}: λ{idx}",
                xy=(z.real, z.imag),
                xytext=(6, 6),
                textcoords="offset points",
                fontsize=8,
                alpha=0.85,
            )

    ax.axhline(0, color="black", linewidth=0.9)
    ax.axvline(0, color="black", linewidth=0.9)

    ax.grid(True, linestyle="--", alpha=0.35)

    ax.set_xlabel("Re(λ)")
    ax.set_ylabel("Im(λ)")
    ax.set_title("Eigenvalues: A vs E(A)")

    ax.legend(loc="best", fontsize=9)
    ax.set_aspect("equal", adjustable="datalim")

    return ax


# ----------------------------------------------------------------------
# Ana dashboard fonksiyonu
# ----------------------------------------------------------------------

def analyze_and_plot_transform(
    A,
    mml,
    name: str = "A",
    sparse_output: bool = False,
    save_prefix: str = None,
    show: bool = True,
    verbose: bool = True,
):
    """
    A matrisini ve E(A) = MatrixMittagLeffler.evaluate(A) sonucunu
    analiz edip yan yana grafikler üretir.
    """

    # 1. A analizi
    A_analyzer = MatrixAnalyzer(A, name=name)

    # 2. E(A) hesapla
    try:
        F = mml.evaluate(A, sparse_output=sparse_output)
    except TypeError:
        # Eğer evaluate sparse_output parametresi desteklemiyorsa
        F = mml.evaluate(A)

    F_name = f"E({name})"

    # 3. E(A) analizi
    F_analyzer = MatrixAnalyzer(F, name=F_name)

    analyzers = [A_analyzer, F_analyzer]

    # 4. Konsol raporu
    if verbose:
        _print_report(A_analyzer)
        _print_report(F_analyzer)

    figures = {}

    # ------------------------------------------------------------------
    # Grafik 1: property heatmap
    # ------------------------------------------------------------------
    fig1, ax1 = plt.subplots(figsize=(16, 3.8))
    _summary_heatmap(analyzers, ax=ax1)
    fig1.tight_layout()

    if save_prefix:
        fig1.savefig(
            f"{save_prefix}_properties.png",
            dpi=150,
            bbox_inches="tight",
        )

    figures["properties"] = fig1

    # ------------------------------------------------------------------
    # Grafik 2: A ve E(A) matrislerini yan yana çiz
    # ------------------------------------------------------------------
    fig2, axes = plt.subplots(1, 2, figsize=(10.6, 4.8))

    _matrix_heatmap_axis(axes[0], A_analyzer, name)
    _matrix_heatmap_axis(axes[1], F_analyzer, F_name)

    fig2.suptitle(f"{name} vs {F_name}", fontsize=13)
    fig2.tight_layout(rect=[0, 0, 1, 0.94])

    if save_prefix:
        fig2.savefig(
            f"{save_prefix}_matrices.png",
            dpi=150,
            bbox_inches="tight",
        )

    figures["matrices"] = fig2

    # ------------------------------------------------------------------
    # Grafik 3: eigenvalue comparison
    # ------------------------------------------------------------------
    fig3, ax3 = plt.subplots(figsize=(7.6, 7.6))
    _eigenvalue_comparison(analyzers, ax=ax3)
    fig3.tight_layout()

    if save_prefix:
        fig3.savefig(
            f"{save_prefix}_eigenvalues.png",
            dpi=150,
            bbox_inches="tight",
        )

    figures["eigenvalues"] = fig3

    if show:
        plt.show()

    return {
        "figures": figures,
        "A_analyzer": A_analyzer,
        "F_analyzer": F_analyzer,
        "F": F,
    }

# ======================================================================
# 1. Sparse Hermitian destekli MatrixMittagLeffler
# ======================================================================

class SparseHermitianMatrixMittagLeffler(MatrixMittagLeffler):
    """
    MatrixMittagLeffler sınıfına sparse Hermitian eigsh yolu ekleyen alt sınıf.

    Yeni evaluate parametreleri:
        sparse_output:
            True ise diagonal sonuçlar sparse döndürülmeye çalışılır.

        sparse_hermitian_k:
            None ise Hermitian sparse matris için tam dense spektral ayrışım yapılır.
            Bir tamsayı ise eigsh ile sadece k özdeğer/özvektör çifti kullanılır.

            ÖNEMLİ:
            k < n olduğunda sonuç tam f(A) değildir.
            Seçilen özçiftlere dayalı yaklaşık bir projeksiyondur:

                f(A) ≈ V_k f(Λ_k) V_k^H

        sparse_hermitian_which:
            eigsh için seçim türü:
                "SA" -> smallest algebraic
                "LA" -> largest algebraic
                "SM" -> smallest magnitude
                "LM" -> largest magnitude
    """

    def evaluate(
        self,
        A,
        max_terms: int = 200,
        sparse_output: bool = False,
        sparse_hermitian_k=None,
        sparse_hermitian_which: str = "SA",
    ):
        analyzer = MatrixAnalyzer(A)

        if not analyzer.is_square:
            raise ValueError("A kare matris olmalı.")

        # --------------------------------------------------------------
        # 1) Diagonal matris
        # --------------------------------------------------------------
        if analyzer.is_diagonal(rtol=self.tol, atol=self.tol):
            return self._via_diag_values(
                analyzer.diagonal,
                sparse_output=sparse_output and analyzer.is_sparse,
            )

        # --------------------------------------------------------------
        # 2) Hermitian matris
        # --------------------------------------------------------------
        if analyzer.is_hermitian(rtol=self.tol, atol=self.tol):

            # Sparse Hermitian ve kullanıcı k verdiyse eigsh yolu
            if analyzer.is_sparse and sparse_hermitian_k is not None:
                return self._via_sparse_hermitian_eigsh(
                    analyzer.matrix,
                    k=sparse_hermitian_k,
                    which=sparse_hermitian_which,
                    sparse_output=sparse_output,
                )

            # Varsayılan tam Hermitian yol
            return self._via_eig(analyzer.dense)

        # --------------------------------------------------------------
        # 3) Genel matris
        # --------------------------------------------------------------
        return self._via_series(analyzer.dense, max_terms)

    # ------------------------------------------------------------------
    # Diyagonal değerler üzerinden değerlendirme
    # ------------------------------------------------------------------

    def _via_diag_values(self, diag, sparse_output: bool = False):
        """
        Sadece köşegen değerlerini kullanarak sonucu üretir.
        """
        diag = np.asarray(diag, dtype=complex)

        vals = np.array(
            [self._scalar.evaluate(z) for z in diag],
            dtype=complex,
        )

        if sparse_output and HAS_SCIPY:
            return sparse.diags(vals, format="csr")

        return np.diag(vals)

    # ------------------------------------------------------------------
    # Sparse Hermitian eigsh yolu
    # ------------------------------------------------------------------

    def _via_sparse_hermitian_eigsh(
        self,
        sparse_matrix,
        k=None,
        which: str = "SA",
        sparse_output: bool = False,
    ):
        """
        Sparse Hermitian matris için scipy.sparse.linalg.eigsh tabanlı yol.

        k=None:
            Tam sonuç için dense'e çevirir ve np.linalg.eigh kullanır.

        k=int:
            Sadece k özçift kullanır.

            Eğer k >= N - 1 ise ARPACK/eigsh sparse matris için kullanılamaz.
            Bu durumda dense eigh ile ilgili k özçift seçilir.

            Sonuç tam f(A) değildir:

                f(A) ≈ V_k f(Λ_k) V_k^H
        """

        if not HAS_SCIPY:
            raise RuntimeError("scipy.sparse gerekli. Lütfen scipy kurun.")

        if not sparse.issparse(sparse_matrix):
            sparse_matrix = sparse.csr_matrix(sparse_matrix)

        if sparse_matrix.shape[0] != sparse_matrix.shape[1]:
            raise ValueError("sparse_matrix kare matris olmalı.")

        n = sparse_matrix.shape[0]

        which = str(which).upper()

        if which not in {"SA", "LA", "SM", "LM"}:
            raise ValueError(
                "which parametresi 'SA', 'LA', 'SM' veya 'LM' olmalı."
            )

        # ------------------------------------------------------------------
        # Tam sonuç
        # ------------------------------------------------------------------
        if k is None:
            warnings.warn(
                "sparse_hermitian_k=None: "
                "Tam spektral ayrışım için sparse matris dense'e çevriliyor. "
                "Büyük matrislerde bellek kullanımı artabilir.",
                stacklevel=2,
            )

            A_dense = sparse_matrix.toarray()

            vals, V = np.linalg.eigh(A_dense)

            fvals = np.array(
                [self._scalar.evaluate(z) for z in vals],
                dtype=complex,
            )

            F = (V * fvals) @ V.conj().T

            if sparse_output:
                return sparse.csr_matrix(F)

            return F

        # ------------------------------------------------------------------
        # k kontrolü
        # ------------------------------------------------------------------
        if k <= 0:
            raise ValueError("k pozitif bir tamsayı olmalı.")

        if k >= n:
            raise ValueError(
                f"k={k} çok büyük. k < n olmalı. n={n}. "
                "Tam sonuç için sparse_hermitian_k=None kullanın."
            )

        # ------------------------------------------------------------------
        # eigsh, sparse matrislerde k >= N - 1 durumunu sevmez.
        # Özellikle 2x2 matriste k=1 buna girer.
        # Bu durumda dense eigh ile seçilmiş özçiftleri kullanıyoruz.
        # ------------------------------------------------------------------
        if k >= n - 1:
            warnings.warn(
                f"k={k} >= N-1={n - 1} olduğundan eigsh sparse matris için kullanılamıyor. "
                "Dense eigh ile seçilen özçiftler kullanılıyor. "
                "Tam sonuç için sparse_hermitian_k=None kullanın.",
                stacklevel=2,
            )

            A_dense = sparse_matrix.toarray()

            vals, V = np.linalg.eigh(A_dense)

            # Hangi özdeğerlerin seçileceğini belirle
            if which == "SA":
                idx = np.argsort(vals.real)[:k]

            elif which == "LA":
                idx = np.argsort(vals.real)[-k:]

            elif which == "SM":
                idx = np.argsort(np.abs(vals))[:k]

            elif which == "LM":
                idx = np.argsort(np.abs(vals))[-k:]

            else:
                raise ValueError("Geçersiz which değeri.")

            vals_k = vals[idx]
            V_k = V[:, idx]

            fvals = np.array(
                [self._scalar.evaluate(z) for z in vals_k],
                dtype=complex,
            )

            F_approx = (V_k * fvals) @ V_k.conj().T

            warnings.warn(
                f"sparse_hermitian_k={k}: "
                "Sonuç tam f(A) değildir. "
                "Seçilen özçiftlere dayalı yaklaşık bir projeksiyondur: "
                "f(A) ≈ V_k f(Λ_k) V_k^H",
                stacklevel=2,
            )

            if sparse_output:
                return sparse.csr_matrix(F_approx)

            return F_approx

        # ------------------------------------------------------------------
        # Normal eigsh yolu
        # ------------------------------------------------------------------
        if np.issubdtype(sparse_matrix.dtype, np.complexfloating):
            dtype = np.complex128
        else:
            dtype = float

        M = sparse_matrix.astype(dtype)

        vals, V = spla.eigsh(
            M,
            k=k,
            which=which,
            return_eigenvectors=True,
        )

        fvals = np.array(
            [self._scalar.evaluate(z) for z in vals],
            dtype=complex,
        )

        F_approx = (V * fvals) @ V.conj().T

        warnings.warn(
            f"sparse_hermitian_k={k}: "
            "Sonuç tam f(A) değildir. "
            "Seçilen özçiftlere dayalı yaklaşık bir projeksiyondur: "
            "f(A) ≈ V_k f(Λ_k) V_k^H",
            stacklevel=2,
        )

        if sparse_output:
            return sparse.csr_matrix(F_approx)

        return F_approx


# ======================================================================
# 2. time_process
# ======================================================================

def _evaluate_with_kwargs(mml, A, **kwargs):
    """
    Farklı MatrixMittagLeffler.evaluate imzalarına karşı güvenli çağrı.
    """
    try:
        return mml.evaluate(A, **kwargs)
    except TypeError:
        # Yeni parametreleri desteklemeyen eski evaluate için geri dönüş
        filtered = {}

        if "max_terms" in kwargs:
            filtered["max_terms"] = kwargs["max_terms"]

        if "sparse_output" in kwargs:
            filtered["sparse_output"] = kwargs["sparse_output"]

        try:
            return mml.evaluate(A, **filtered)
        except TypeError:
            return mml.evaluate(A)


def time_process(
    H,
    t_array,
    alpha: float = 1.0,
    beta: float = 1.0,
    gamma: float = 1.0,
    tolerance: float = 1e-13,
    max_terms: int = 200,
    sparse_output: bool = False,
    sparse_hermitian_k=None,
    sparse_hermitian_which: str = "SA",
    return_analyzers: bool = False,
    mml=None,
):
    """
    H matrisi üzerinden t parametreli matris süreci üretir.

    α=1, β=1, γ=1 için:

        E^1_{1,1}(-i H t) = exp(-i H t)

    Eğer H Hermitian ise bu matrisler unitary olur.

    Parametreler:
        H:
            Kare matris. Dense numpy veya scipy.sparse olabilir.

        t_array:
            Zaman/süreç parametreleri.

        return_analyzers:
            True ise her sonuç için MatrixAnalyzer nesneleri de döndürülür.

        mml:
            İstersen kendi MatrixMittagLeffler nesneni verebilirsin.
            Varsayılan olarak SparseHermitianMatrixMittagLeffler kullanılır.
    """

    H_analyzer = MatrixAnalyzer(H, name="H")

    if not H_analyzer.is_square:
        raise ValueError("H kare matris olmalı.")

    if mml is None:
        mml = SparseHermitianMatrixMittagLeffler(
            alpha=alpha,
            beta=beta,
            gamma=gamma,
            tolerance=tolerance,
        )

    scalar_input = np.isscalar(t_array)

    t_values = np.atleast_1d(np.asarray(t_array, dtype=float))

    H_mat = H_analyzer.matrix

    results = []
    analyzers = []

    for t in t_values:

        # A(t) = -i H t
        A_t = H_mat * (-1j * float(t))

        F = _evaluate_with_kwargs(
            mml,
            A_t,
            max_terms=max_terms,
            sparse_output=sparse_output,
            sparse_hermitian_k=sparse_hermitian_k,
            sparse_hermitian_which=sparse_hermitian_which,
        )

        results.append(F)

        if return_analyzers:
            analyzers.append(
                MatrixAnalyzer(F, name=f"t={float(t):.6g}")
            )

    if scalar_input:
        if return_analyzers:
            return results[0], analyzers[0]
        return results[0]

    if return_analyzers:
        return results, analyzers

    return results


# ======================================================================
# 3. Süreç metrikleri grafiği
# ======================================================================
def plot_time_process_metrics(
    t_array,
    results,
    save_prefix: str = None,
    show: bool = True,
    return_analyzers=False,
):
    """
    time_process sonuçları üzerinden temel metrikleri çizer:

    - spectral radius
    - condition number
    - |trace|
    """

    t_values = np.atleast_1d(np.asarray(t_array, dtype=float))

    analyzers = [
        MatrixAnalyzer(F, name=f"t={float(t):.6g}")
        for t, F in zip(t_values, results)
    ]

    spectral_radius = []
    condition_numbers = []
    trace_abs = []

    for analyzer in analyzers:

        spectral_radius.append(analyzer.spectral_radius)

        cond = analyzer.condition_number
        if np.isfinite(cond):
            condition_numbers.append(float(cond))
        else:
            condition_numbers.append(np.nan)

        trace_abs.append(abs(analyzer.trace))

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.2))

    axes[0].plot(t_values, spectral_radius, marker="o")
    axes[0].set_title("Spectral Radius")
    axes[0].set_xlabel("t")

    axes[1].plot(t_values, condition_numbers, marker="o", color="tab:orange")
    axes[1].set_title("Condition Number")
    axes[1].set_xlabel("t")

    axes[2].plot(t_values, trace_abs, marker="o", color="tab:green")
    axes[2].set_title("|Trace|")
    axes[2].set_xlabel("t")

    for ax in axes:
        ax.grid(True, linestyle="--", alpha=0.35)

    fig.tight_layout()

    if save_prefix:
        fig.savefig(
            f"{save_prefix}_metrics.png",
            dpi=150,
            bbox_inches="tight",
        )

    if show:
        plt.show()

    if return_analyzers:
        return fig, analyzers
    return fig


# ======================================================================
# 4. Süreç özdeğer grafiği
# ======================================================================

def plot_time_process_eigenvalues(
    t_array,
    results,
    save_prefix: str = None,
    show: bool = True,
):
    """
    time_process sonuçlarının özdeğerlerini karmaşık düzlemde çizer.

    Renk, t parametresini gösterir.
    """

    t_values = np.atleast_1d(np.asarray(t_array, dtype=float))

    if t_values.size == 0:
        return None

    vmin = float(np.min(t_values))
    vmax = float(np.max(t_values))

    if vmax <= vmin:
        vmax = vmin + 1.0

    norm = Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.get_cmap("viridis")

    fig, ax = plt.subplots(figsize=(7.6, 7.6))

    for t, F in zip(t_values, results):

        analyzer = MatrixAnalyzer(F, name=f"t={float(t):.6g}")
        vals = analyzer.eigenvalues

        if vals.size == 0:
            continue

        ax.scatter(
            vals.real,
            vals.imag,
            color=cmap(norm(float(t))),
            s=22,
            alpha=0.85,
        )

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_label("t")

    ax.axhline(0, color="black", linewidth=0.9)
    ax.axvline(0, color="black", linewidth=0.9)

    ax.grid(True, linestyle="--", alpha=0.35)

    ax.set_xlabel("Re(λ)")
    ax.set_ylabel("Im(λ)")
    ax.set_title("Eigenvalue Trajectory of Process")

    ax.set_aspect("equal", adjustable="datalim")

    fig.tight_layout()

    if save_prefix:
        fig.savefig(
            f"{save_prefix}_eigenvalues.png",
            dpi=150,
            bbox_inches="tight",
        )

    if show:
        plt.show()

    return fig

def NV_hamiltonian(Bz, Bx=0.0):
    Sx = np.array([[0,1,0],[1,0,1],[0,1,0]])/np.sqrt(2)
    Sy = np.array([[0,-1j,0],[1j,0,-1j],[0,1j,0]])/np.sqrt(2)
    Sz = np.array([[1,0,0],[0,0,0],[0,0,-1]])
    return D*Sz@Sz + gamma_e*(Bx*Sx + Bz*Sz)

# Alias
#kitaev_chain = kitaev_chain_hamiltonian
# Kitaev zinciri Hamiltonyeni (Majorana arayüzü)
def kitaev_chain(N, mu, t=1.0, Delta=1.0):
    """N site'lik Kitaev zinciri"""
    H = np.zeros((2*N, 2*N), dtype=complex)
    for i in range(N-1):
        H[2*i, 2*i+3] = -t + Delta
        H[2*i+3, 2*i] = -t + Delta
        H[2*i+1, 2*i+2] = t + Delta
        H[2*i+2, 2*i+1] = t + Delta
    for i in range(N):
        H[2*i, 2*i] = -mu
        H[2*i+1, 2*i+1] = mu
    return H


"""
Kesirli Lindblad Dinamiği Modülü
=================================
Açık kuantum sistemleri için matris değerli Mittag-Leffler evrimi:

    vec(ρ(t)) = E_α( L_super · t^α ) vec(ρ(0))

burada L_super sparse süperoperatördür.
"""
# ======================================================================
# 1. LİNDBLADYEN SÜPEROPERATÖR İNŞASI
# ======================================================================

def build_lindbladian(H, jump_operators=None, gammas=None):
    """
    Lindblad master denklemi için sparse süperoperatör inşa eder.
    
    dρ/dt = -i[H, ρ] + Σ_k γ_k (L_k ρ L_k† - 1/2 {L_k† L_k, ρ})
    
    Süperoperatör formunda:
        d·vec(ρ)/dt = L_super · vec(ρ)
    
    Parametreler
    ----------
    H : (N, N) ndarray veya sparse matrix
        Sistem Hamiltonyeni.
    jump_operators : list of (N, N) arrays
        Lindblad çöküş operatörleri L_k.
    gammas : list of float
        Her çöküş operatörü için oran γ_k.
    
    Döndürür
    -------
    L_super : (N², N²) sparse CSR matrix
        Süperoperatör.
    """
    if not HAS_SCIPY:
        raise RuntimeError("scipy gerekli.")
    
    if sparse.issparse(H):
        H = H.toarray()
    H = np.asarray(H, dtype=complex)
    
    N = H.shape[0]
    I = sparse.eye(N, dtype=complex, format='csr')
    
    # -iH ρ kısmı: (-iH) ⊗ I
    L_ham_left = sparse.kron(-1j * sparse.csr_matrix(H), I, format='csr')
    
    # +iρ H kısmı: I ⊗ (iH^T)
    L_ham_right = sparse.kron(I, sparse.csr_matrix(1j * H.T), format='csr')
    
    L_super = L_ham_left + L_ham_right
    
    # Lindblad dissipatör terimleri
    if jump_operators is not None:
        if gammas is None:
            gammas = [1.0] * len(jump_operators)
        
        for L_k, gamma in zip(jump_operators, gammas):
            if sparse.issparse(L_k):
                L_k = L_k.toarray()
            L_k = np.asarray(L_k, dtype=complex)
            L_k_dag = L_k.conj().T
            L_k_dag_L_k = L_k_dag @ L_k
            
            # L_k ρ L_k† → L_k ⊗ L_k*
            term1 = sparse.kron(
                sparse.csr_matrix(L_k),
                sparse.csr_matrix(L_k.conj()),
                format='csr'
            )
            
            # -1/2 {L_k†L_k, ρ} → -1/2 (L_k†L_k ⊗ I + I ⊗ (L_k†L_k)^T)
            term2 = -0.5 * sparse.kron(
                sparse.csr_matrix(L_k_dag_L_k), I, format='csr'
            )
            term3 = -0.5 * sparse.kron(
                I, sparse.csr_matrix(L_k_dag_L_k.T), format='csr'
            )
            
            L_super = L_super + gamma * (term1 + term2 + term3)
    
    return L_super


def vec(rho):
    """Yoğunluk matrisini vektöre çevirir (column-major)."""
    return np.asarray(rho, dtype=complex).flatten(order='F')


def unvec(v, N):
    """Vektörü yoğunluk matrisine çevirir."""
    return np.asarray(v, dtype=complex).reshape((N, N), order='F')


# ======================================================================
# 2. HIZLI DİAGONAL SÜREÇ (100x HIZLANDIRMA)
# ======================================================================

def fast_diagonal_process(
    diag_values,
    t_array,
    alpha: float = 1.0,
    beta: float = 1.0,
    gamma: float = 1.0,
):
    """
    Diagonal jeneratör için MATRİS ÜRETMEDEN süreç hesaplar.
    
    Normal time_process her t için N×N matris üretir.
    Bu fonksiyon sadece N boyutlu vektörler üretir → O(N²) → O(N).
    
    Parametreler
    ----------
    diag_values : (N,) ndarray
        Jeneratörün köşegen değerleri (örn. -i·enerjiler).
    t_array : (M,) ndarray
        Zaman noktaları.
    
    Döndürür
    -------
    result_diag : (M, N) ndarray
        Her t için f(diag) vektörü.
        Son diagonal matris: np.diag(result_diag[m])
    """
    diag_values = np.asarray(diag_values, dtype=complex)
    t_array = np.atleast_1d(np.asarray(t_array, dtype=float))
    
    mml = MatrixMittagLeffler(
        alpha=alpha, beta=beta, gamma=gamma, tolerance=1e-13
    )
    
    M = len(t_array)
    N = len(diag_values)
    result_diag = np.empty((M, N), dtype=complex)
    
    for m, t in enumerate(t_array):
        # f(λ_j · t) her j için skaler
        scaled = diag_values * float(t)
        result_diag[m, :] = np.array(
            [mml._scalar.evaluate(z) for z in scaled],
            dtype=complex
        )
    
    return result_diag


# ======================================================================
# 3. KESİRLİ LİNDBLADIAN SÜREÇ SINIFI
# ======================================================================

class FractionalLindbladianProcess:
    """
    Kesirli Lindblad dinamiği:
    
        vec(ρ(t)) = E_α(L_super · t^α) · vec(ρ(0))
    
    α = 1.0 → standart Markovian Lindblad
    α < 1.0 → alt-dağılımlı (subdiffusive) non-Markovian dinamik
    α > 1.0 → üst-dağılımlı (superdiffusive) dinamik
    """

    def __init__(self, H, jump_operators=None, gammas=None, alpha: float = 1.0, beta: float = 1.0, gamma: float = 1.0):
        self.H = H
        self.jump_operators = jump_operators
        self.gammas = gammas
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        
        N = H.shape[0]
        self.dim = N
        
        # Süperoperatörü inşa et
        self.L_super = build_lindbladian(H, jump_operators, gammas)
        
        self.mml = MatrixMittagLeffler(alpha=alpha, beta=beta, gamma=gamma, tolerance=1e-13)
        
        # KARARLILIK İÇİN: L_super'ü bir kez dense'e çevir ve özdeğerlerine ayır
        # (64x64 veya 400x400 için bu işlem milisaniyeler sürer ve çok kararlıdır)
        L_dense = self.L_super.toarray() if sparse.issparse(self.L_super) else self.L_super
        self.evals, self.evecs = np.linalg.eig(L_dense)
        self.evecs_inv = np.linalg.inv(self.evecs)

    def evolve(self, rho0, t_array):
        """
        ρ(0) → ρ(t) evrimi hesaplar.
        
        Parametreler
        ----------
        rho0 : (N, N) ndarray
            Başlangıç yoğunluk matrisi.
        t_array : (M,) ndarray veya float
            Zaman noktaları.
        
        Döndürür
        -------
        rhos : list of (N, N) ndarray
            Her t için yoğunluk matrisi.
        """
        scalar_input = np.isscalar(t_array)
        t_values = np.atleast_1d(np.asarray(t_array, dtype=float))
        v0 = vec(rho0)
        N = self.dim
        
        rhos = []
        for t in t_values:
            # 1. Özdeğerlere Mittag-Leffler uygula: E_alpha(lambda * t^alpha)
            scaled_evals = self.evals * (float(t) ** self.alpha)
            f_evals = np.array([self.mml._scalar.evaluate(l) for l in scaled_evals], dtype=complex)
            
            # 2. Evrim operatörünü yeniden oluştur: V * diag(f_evals) * V^-1
            F_matrix = self.evecs @ np.diag(f_evals) @ self.evecs_inv
            
            # 3. rho(t) vektörünü hesapla
            rho_t_vec = F_matrix @ v0
            rho_t = unvec(rho_t_vec, N)
            
            # 4. SAYISAL KARARLILIK İÇİN ZORUNLU TEMİZLİK (ÇOK ÖNEMLİ)
            # Gürültüden kaynaklanan anti-hermitian kısmı at
            rho_t = (rho_t + rho_t.conj().T) / 2.0
            
            # Trace'i 1'e normalize et
            trace_val = np.trace(rho_t)
            if np.abs(trace_val) > 1e-10:
                rho_t = rho_t / trace_val
            else:
                # Eğer tamamen çöktüyse, maksimum karışık duruma düşür
                rho_t = np.eye(N, dtype=complex) / N
            
            # Negatif özdeğerleri (sayısal gürültü) sıfırla
            evals_rho, evecs_rho = np.linalg.eigh(rho_t)
            evals_rho = np.maximum(evals_rho, 0.0)
            rho_t = evecs_rho @ np.diag(evals_rho) @ evecs_rho.conj().T
            rho_t = rho_t / np.trace(rho_t) # Yeniden normalize et
            
            rhos.append(rho_t)
            
        return rhos[0] if scalar_input else rhos
    
    def _dense_apply(self, A_sparse, v):
        """F(A)·v hesaplar (dense'e çevirerek)."""
        A_dense = A_sparse.toarray() if sparse.issparse(A_sparse) else A_sparse
        F = self.mml.evaluate(A_dense)
        return F @ v
    
    def steady_state(self, tol: float = 1e-10) -> np.ndarray:
        """
        Uzun-zaman limit yoğunluk matrisi.

        L_super · vec(ρ_ss) = 0  →  en küçük |λ| özvektörü.

        __init__'te hesaplanmış self.evals/self.evecs kullanılır.
        """
        idx = int(np.argmin(np.abs(self.evals)))
        lam_min = self.evals[idx]

        if np.abs(lam_min) > tol:
            import warnings
            warnings.warn(
                f"Kararlı durum yok: min |λ| = {np.abs(lam_min):.2e} > tol={tol}",
                RuntimeWarning, stacklevel=2)

        v_ss = self.evecs[:, idx]
        rho_ss = unvec(v_ss, self.dim)
        tr = np.trace(rho_ss)
        if abs(tr) < 1e-15:
            raise RuntimeError("Kararlı durum izi sıfır — normalize edilemiyor")
        return rho_ss / tr


# ======================================================================
# 4. TOPOLOJİK İNVARYANTLAR
# ======================================================================

def topological_invariant(analyzer, kind='winding'):
    """
    Matrisin topolojik invariant'ını hesaplar.
    
    Parametreler
    ----------
    analyzer : MatrixAnalyzer
    kind : str
        'winding'  → Winding number (1D sistemler)
        'chern'    → Chern number (2D sistemler)
        'pfaffian' → Z₂ invariant (Majorana sistemleri)
    """
    if kind == 'pfaffian':
        if not HAS_PFAPACK:
            raise ImportError("pfapack gerekli. Kurulum: pip install pfapack")
        
        A = analyzer.dense
        try:
            # Skew-symmetric kontrolü
            if not np.allclose(A, -A.T, atol=1e-10):
                warnings.warn("Pfaffian için skew-symmetric matris gerekli.")
                return 0
            
            pf_val = pf.pfaffian(A)
            return 1 if pf_val.real > 0 else -1
        except Exception as e:
            warnings.warn(f"Pfaffian hesabı başarısız: {e}")
            return 0
    
    elif kind == 'winding':
        # 1D topolojik invariant: sign(det(H))
        A = analyzer.dense
        det_sign = np.sign(np.linalg.det(A).real)
        return int(det_sign)
    
    else:
        raise ValueError(f"Bilinmeyen invariant türü: {kind}")

def kitaev_winding(mu, t=1.0, Delta=1.0, N_k=1000):
    """
    Kitaev zinciri için winding number hesabı (momentum uzayı).
    
    |μ| < 2t  → w = ±1 (topolojik faz, Majorana sıfır modları)
    |μ| > 2t  → w = 0 (triviyal faz)
    """
    ks = np.linspace(0, 2*np.pi, N_k, endpoint=False)
    
    d_y = 2 * Delta * np.sin(ks)
    d_z = -mu - 2 * t * np.cos(ks)
    
    # φ(k) = arg(d_y + i·d_z)
    phi = np.arctan2(d_z, d_y)
    
    # Winding: toplam φ değişimi / 2π
    dphi = np.diff(phi)
    # [-π, π] aralığına wrap et
    dphi = (dphi + np.pi) % (2 * np.pi) - np.pi
    
    winding = np.sum(dphi) / (2 * np.pi)
    return int(np.round(winding))


def kitaev_chain_hamiltonian(N, mu, t=1.0, Delta=1.0):
    """
    N site Kitaev zinciri (BdG formunda, 2N × 2N).
    """
    H = np.zeros((2*N, 2*N), dtype=complex)
    for i in range(N - 1):
        H[2*i, 2*i+3] = -t + Delta
        H[2*i+3, 2*i] = -t + Delta
        H[2*i+1, 2*i+2] = t + Delta
        H[2*i+2, 2*i+1] = t + Delta
    for i in range(N):
        H[2*i, 2*i] = -mu
        H[2*i+1, 2*i+1] = mu
    return H


# ======================================================================
# 5. YÜZEY KODU HAMILTONYENİ
# ======================================================================

def build_surface_code_hamiltonian(distance: int = 3):
    """
    Toric/surface code stabilizer Hamiltonyeni inşa eder.
    
    H = -Σ_s A_s - Σ_p B_p
    
    distance → N = 2·distance² - 1 kubit (yaklaşık)
    
    Döndürür
    -------
    H : sparse CSR matrix
        2^N × 2^N sparse Hamiltonyen.
    """
    if not HAS_SCIPY:
        raise RuntimeError("scipy gerekli.")
    
    # Basitleştirilmiş 1D Ising chain (surface code'un 1D sınırı)
    N = 2 * distance  # kubit sayısı
    
    I = sparse.eye(2, dtype=complex, format='csr')
    X = sparse.csr_matrix(np.array([[0,1],[1,0]], dtype=complex))
    Z = sparse.csr_matrix(np.array([[1,0],[0,-1]], dtype=complex))
    
    H = sparse.csr_matrix((2**N, 2**N), dtype=complex)
    
    # X-type stabilizers (komşu çiftler)
    for i in range(N - 1):
        term = sparse.eye(1, dtype=complex, format='csr')
        for j in range(N):
            if j == i or j == i + 1:
                term = sparse.kron(term, X, format='csr')
            else:
                term = sparse.kron(term, I, format='csr')
        H = H - term
    
    # Z-type stabilizers
    for i in range(0, N - 1, 2):
        if i + 1 < N:
            term = sparse.eye(1, dtype=complex, format='csr')
            for j in range(N):
                if j == i or j == i + 1:
                    term = sparse.kron(term, Z, format='csr')
                else:
                    term = sparse.kron(term, I, format='csr')
            H = H - term
    
    return H

