import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Linarith

-- Algebra behind the second order inversion of G(y)=y(1-f(y)).
theorem peeling_taylor_identity (Z F Fp Fpp s : Real) :
    (Z-s)*(1-(F-Fp*s+Fpp/2*s^2)) =
    Z*(1-F) - (1-F-Z*Fp)*s - (Fp+Z*Fpp/2)*s^2 + Fpp/2*s^3 := by
  ring

theorem peeling_critical_taylor (Z F Fp Fpp s : Real)
    (h0 : Z*(1-F)=1) (h1 : 1-F-Z*Fp=0) :
    (Z-s)*(1-(F-Fp*s+Fpp/2*s^2)) =
    1-(Fp+Z*Fpp/2)*s^2+Fpp/2*s^3 := by
  rw [peeling_taylor_identity, h0, h1]
  ring

theorem peeling_quadratic_positive (Z Fpp : Real)
    (hZ : 0 ≤ Z) (hFpp : 0 ≤ Fpp) :
    1 ≤ 1 + Z^3*Fpp/2 := by
  have hp : 0 ≤ Z^3*Fpp := mul_nonneg (pow_nonneg hZ 3) hFpp
  linarith

-- Pointwise denominator comparison used before summing the kernel.
theorem kernel_denominator_comparison (m k p : Real)
    (hm : 0 ≤ m) (hkp : k ≤ p) :
    k * (m+p) ≤ p * (m+k) := by
  nlinarith [mul_nonneg hm (sub_nonneg.mpr hkp)]

theorem kernel_large_atom_comparison (m p : Real)
    (hp : p ≤ 2*m) : m+p ≤ 3*m := by
  linarith

-- Recurrence translating the partition-function integral into the kernel.
theorem kernel_factor_identity (k a b : Real)
    (hk : k ≠ 0) (hden : 2*k-1 ≠ 0)
    (hrec : (2*k)*b = (2*k-1)*a) :
    (2*k-1)*(a/(2*k)) = b := by
  rw [← mul_div_assoc]
  apply (div_eq_iff (mul_ne_zero (by norm_num) hk)).mpr
  nlinarith
