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
from typing import ClassVar, Final

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.special import (
    gamma as gamma_func,
    digamma,
    gammaln,
    rgamma,
    erfcx, poch
)


import mpmath as _mp
# Not: global mp.dps DEĞİŞTİRİLMEZ. _mpmath_fallback kendi
# workdps(dps) bağlamını kullanır; çağıran tarafın mpmath
# hassasiyetini etkilemez.


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
