#include "gmp_math.h"

#include <cmath>
#include <algorithm>
#include <stdexcept>

mpz_class mpz_sqrt(const mpz_class& x) { return sqrt(x); }

mpz_pair mpz_sqrtrem(const mpz_class& x)
{
  mpz_pair result;
  mpz_sqrtrem(result.first.get_mpz_t(), result.second.get_mpz_t(),
              x.get_mpz_t());
  return result;
}

mpz_class mpz_root(const mpz_class& x, unsigned long y)
{
  mpz_class result;
  mpz_root(result.get_mpz_t(), x.get_mpz_t(), y);
  return result;
}

mpz_pair mpz_rootrem(const mpz_class& x, unsigned long y)
{
  mpz_pair result;
  mpz_rootrem(result.first.get_mpz_t(), result.second.get_mpz_t(),
              x.get_mpz_t(), y);
  return result;
}

bool mpz_perfect_power(const mpz_class& x)
{
  return mpz_perfect_power_p(x.get_mpz_t());
}

bool mpz_perfect_square(const mpz_class& x)
{
  return mpz_perfect_square_p(x.get_mpz_t());
}

mpz_class mpz_gcd(const mpz_class& x, const mpz_class& y) { return gcd(x, y); }

mpz_class mpz_gcdext(const mpz_class& x, const mpz_class& y,
                     mpz_class& a, mpz_class& b)
{
  mpz_class result;
  mpz_gcdext(result.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t(),
             x.get_mpz_t(), y.get_mpz_t());
  return result;
}

mpz_class mpz_lcm(const mpz_class& x, const mpz_class& y) { return lcm(x, y); }

mpz_class mpz_pow(const mpz_class& x, unsigned long y)
{
  mpz_class result;
  mpz_pow_ui(result.get_mpz_t(), x.get_mpz_t(), y);
  return result;
}

mpz_class mpz_powm(const mpz_class& x, const mpz_class& y,
                   const mpz_class& m)
{
  mpz_class result;
  mpz_powm(result.get_mpz_t(), x.get_mpz_t(), y.get_mpz_t(), m.get_mpz_t());
  return result;
}

mpf_class mpf_e(mp_bitcnt_t prec)
{
  mpz_class denom(1);
  mpf_class prv(0, prec), result(2, prec), inc(1, prec);
  do {
    prv = result;
    result += inc /= ++denom;
  } while (prv != result);
  return result;
}

mpf_class mpf_pi(mp_bitcnt_t prec)
{
  mpf_class a(1, prec), b(mpf_class(1, prec) / sqrt(mpf_class(2, prec))), t("0.25", prec), p(a);
  while (a != b) {
    mpf_class A((a + b) / 2), B(sqrt(a * b)), T(t - p * (a - A) * (a - A));
    a = A, b = B, t = T; p *= 2;
  }
  return (a + b) * (a + b) / 4 / t;
}

mpf_class mpf_sqrt(const mpf_class& x) { return sqrt(x); }

mpf_class mpf_cbrt(const mpf_class& x)
{
  long exp_part = mpf_log2floor(x) / 3;
  mpf_class prv(x), now(x), result(mpf_mulexp2(
      mpf_class(cbrt(mpf_mulexp2(x, -exp_part * 3).get_d()),
                x.get_prec()), exp_part));
  do {
    prv = now;
    result -= now = (result - x / (result * result)) / 3;
  } while (abs(now) < abs(prv));
  return result;
}

mpf_class mpf_root(const mpf_class& x, unsigned long y)
{
  long exp_part = mpf_log2floor(x) / (long)y;
  mpf_class prv(x), now(x), result(mpf_mulexp2(
      mpf_class(pow(mpf_mulexp2(x, -exp_part * y).get_d(), 1. / y),
                x.get_prec()), exp_part));
  do {
    prv = now;
    result -= now = (result - x / mpf_powint(result, y - 1)) / y;
  } while (abs(now) < abs(prv));
  return result;
}

mpf_class mpf_mulexp2(const mpf_class& x, long y)
{
  if (!y) return x;
  mpf_class result(0, x.get_prec());
  if (y < 0) {
    mpf_div_2exp(result.get_mpf_t(), x.get_mpf_t(), -y);
  } else {
    mpf_mul_2exp(result.get_mpf_t(), x.get_mpf_t(), y);
  }
  return result;
}

mpf_class _mpf_sin_taylor(mpf_class x)
{
  mpz_class counter(2);
  mpf_class prv(x), result(x), denom(1, x.get_prec()), x2(x * x);
  do {
    prv = result;
    denom *= -counter * (counter + 1);
    counter += 2;
    x *= x2;
    result += x / denom;
  } while (prv != result);
  return result;
}

mpf_class mpf_sin(const mpf_class& x)
{
  mpf_class pi2 = 2 * mpf_pi(x.get_prec());
  return _mpf_sin_taylor(x - floor(x / pi2 + .5) * pi2);
}

mpf_class _mpf_cos_taylor(mpf_class x)
{
  mpz_class counter(1);
  mpf_class prv(x), result(1, x.get_prec());
  mpf_class denom(1, x.get_prec()), x2(x * x);
  x = 1;
  do {
    prv = result;
    denom *= -counter * (counter + 1);
    counter += 2;
    x *= x2;
    result += x / denom;
  } while (prv != result);
  return result;
}

mpf_class mpf_cos(const mpf_class& x)
{
  mpf_class pi2 = 2 * mpf_pi(x.get_prec());
  return _mpf_cos_taylor(x - floor(x / pi2 + .5) * pi2);
}

mpf_class mpf_tan(const mpf_class& x)
{
  mpf_class pi2 = 2 * mpf_pi(x.get_prec());
  mpf_class xp = x - floor(x / pi2) * pi2;
  return _mpf_sin_taylor(xp) / _mpf_cos_taylor(xp);
}

mpf_class _mpf_arcsin_taylor_zero(mpf_class x)
{
  mpz_class denom(1);
  mpf_class prv(x), result(x), x2(x * x);
  do {
    prv = result;
    x *= x2 * denom / (denom + 1);
    denom += 2;
    result += x / denom;
  } while (prv != result);
  return result;
}

mpf_class _mpf_arcsin_half_pi_taylor_minus_one(mpf_class x)
{
  x += 1;
  mpz_class denom(1);
  mpf_class prv(x), result(1, x.get_prec()), num(1, x.get_prec());
  do {
    prv = result;
    num *= x * denom / (2 * (denom + 1));
    denom += 2;
    result += num / denom;
  } while (prv != result);
  return result * sqrt(2 * x);
}

mpf_class mpf_asin(const mpf_class& x)
{
  const double thresh = 0.5;
  if (x < thresh && x > -thresh)
    return _mpf_arcsin_taylor_zero(x);
  else if (x < 0)
    return _mpf_arcsin_half_pi_taylor_minus_one(x) - mpf_pi(x.get_prec()) / 2;
  else
    return mpf_pi(x.get_prec()) / 2 - _mpf_arcsin_half_pi_taylor_minus_one(-x);
}

mpf_class mpf_acos(const mpf_class& x)
{
  const double thresh = 0.5;
  if (x < thresh && x > -thresh)
    return mpf_pi(x.get_prec()) / 2 - _mpf_arcsin_taylor_zero(x);
  else if (x < 0)
    return mpf_pi(x.get_prec()) - _mpf_arcsin_half_pi_taylor_minus_one(x);
  else
    return _mpf_arcsin_half_pi_taylor_minus_one(-x);
}

mpf_class _mpf_arctan_taylor_zero(mpf_class x)
{
  mpz_class denom(1);
  mpf_class prv(x), result(x), nx2(-x * x);
  do {
    prv = result;
    x *= nx2;
    denom += 2;
    result += x / denom;
  } while (prv != result);
  return result;
}

mpf_class _mpf_arctan_taylor_one(mpf_class x)
{
  x -= 1;
  mpz_class denom_odd(1), denom_even(1);
  mpf_class prv(x), result(mpf_pi(x.get_prec()) / 4);
  mpf_class x2(x * x), num(x), pow2(2, x.get_prec());
  do {
    prv = result;
    result += num / (denom_odd * pow2);
    num *= x;
    result -= num / (denom_even * pow2 * 2);
    num *= x;
    result += num / ((denom_odd + 2) * pow2 * 2);
    num *= x2;
    denom_odd += 4;
    denom_even += 2;
    pow2 *= -4;
  } while (prv != result);
  return result;
}

mpf_class mpf_atan(const mpf_class& x)
{
  const double thresh = 0.6;
  if (x > 1 / thresh)
    return mpf_pi(x.get_prec()) / 2 - _mpf_arctan_taylor_zero(1 / x);
  else if (x < -1 / thresh)
    return -mpf_pi(x.get_prec()) / 2 - _mpf_arctan_taylor_zero(1 / x);
  else if (x > -thresh && x < thresh)
    return _mpf_arctan_taylor_zero(x);
  else if (x < 0)
    return -_mpf_arctan_taylor_one(-x);
  else
    return _mpf_arctan_taylor_one(x);
}

mpf_class mpf_atan2(const mpf_class& Y, const mpf_class& X)
{
  if (X > 0)
    return mpf_atan(Y / X);
  mpf_class pi = mpf_pi(std::max(X.get_prec(), Y.get_prec()));
  if (X == 0) {
    if (Y < 0) return -pi / 2;
    return pi / 2;
  }
  if (Y < 0) return mpf_atan(Y / X) - pi;
  return mpf_atan(Y / X) + pi;
}

long mpf_log2floor(const mpf_class& x)
{
  long result;
  mpf_get_d_2exp(&result, x.get_mpf_t());
  return result - 1;
}

mpf_class _mpf_ln1_taylor(const mpf_class& x)
{
  mpz_class denom(1);
  mpf_class prv(x), result(x), num(x);
  do {
    prv = result;
    num *= -x;
    result += num / ++denom;
  } while (prv != result);
  return result;
}

mpf_class mpf_ln(const mpf_class& x)
{
  if (x <= 0) throw std::domain_error("mpf_ln: non-positive logarithm argument");
  long exp_part = mpf_log2floor(x) + 1;
  mpf_class nhalf("-0.5", x.get_prec());
  return _mpf_ln1_taylor(mpf_mulexp2(x, -exp_part) - 1)
          - exp_part * _mpf_ln1_taylor(nhalf);
}

mpf_class mpf_log2(const mpf_class& x)
{
  if (x <= 0) throw std::domain_error("mpf_log2: non-positive logarithm argument");
  long exp_part = mpf_log2floor(x) + 1;
  mpf_class nhalf("-0.5", x.get_prec());
  return exp_part - _mpf_ln1_taylor(mpf_mulexp2(x, -exp_part) - 1)
          / _mpf_ln1_taylor(nhalf);
}

mpf_class _mpf_exp_taylor(mpf_class x)
{
  mpz_class denom(1);
  mpf_class prv(x), result(x + 1), num(x);
  do {
    prv = result;
    result += num *= x / ++denom;
  } while (prv != result);
  return result;
}

mpf_class _mpf_exp2(mpf_class x, const mpf_class& ln2)
{
  if (!x.fits_slong_p())
    throw std::out_of_range("_mpf_exp2: exponent too small or too large");
  mpf_class fl = floor(x);
  return mpf_mulexp2(_mpf_exp_taylor((x - fl) * ln2), fl.get_si());
}

mpf_class mpf_exp(const mpf_class& x)
{
  mpf_class ln2(-_mpf_ln1_taylor(mpf_class("-0.5", x.get_prec())));
  return _mpf_exp2(x / ln2, ln2);
}

mpf_class mpf_exp2(const mpf_class& x)
{
  mpf_class ln2(-_mpf_ln1_taylor(mpf_class("-0.5", x.get_prec())));
  return _mpf_exp2(x, ln2);
}

mpf_class mpf_powint(mpf_class x, long y)
{
  bool flag = false;
  if (y < 0) y = -y, flag = true;

  mpf_class result(1, x.get_prec());
  for (; y; (y >>= 1) && (x *= x))
    if (y & 1l) result *= x;
  return flag ? 1 / result : result;
}

mpf_class mpf_pow(const mpf_class& x, const mpf_class& y)
{
  if (floor(y) == y) {
    if (x == 0 && y == 0) throw std::domain_error("mpf_pow: zero base and exponent");
    if (!y.fits_slong_p())
      throw std::out_of_range("mpf_pow: exponent too small or too large");
    if (y.get_prec() > x.get_prec())
      return mpf_powint(mpf_class(x, y.get_prec()), y.get_si());
    return mpf_powint(x, y.get_si());
  }
  if (x < 0) throw std::domain_error("mpf_pow: negative base and non-integer exponent");
  if (x == 0) return mpf_class(0, std::max(x.get_prec(), y.get_prec()));

  mpf_class ln2(-_mpf_ln1_taylor(mpf_class("-0.5", x.get_prec())));
  return _mpf_exp2(y * mpf_log2(x), ln2);
}
