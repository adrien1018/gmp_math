#include "gmp_math.h"

#include <algorithm>
#include <stdexcept>

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
  long exp_part;
  mpf_get_d_2exp(&exp_part, x.get_mpf_t());
  mpf_t norm_x_t;
  mpf_init2(norm_x_t, x.get_prec());
  if (exp_part >= 0) {
    mpf_div_2exp(norm_x_t, x.get_mpf_t(), exp_part);
  } else {
    mpf_mul_2exp(norm_x_t, x.get_mpf_t(), -exp_part);
  }

  mpf_class norm_x(norm_x_t), nhalf("-0.5", x.get_prec());
  mpf_clear(norm_x_t);
  return _mpf_ln1_taylor(norm_x - 1) - exp_part * _mpf_ln1_taylor(nhalf);
}

mpf_class mpf_log2(const mpf_class& x)
{
  if (x <= 0) throw std::domain_error("mpf_log2: non-positive logarithm argument");
  long exp_part;
  mpf_get_d_2exp(&exp_part, x.get_mpf_t());
  mpf_t norm_x_t;
  mpf_init2(norm_x_t, x.get_prec());
  if (exp_part >= 0) {
    mpf_div_2exp(norm_x_t, x.get_mpf_t(), exp_part);
  } else {
    mpf_mul_2exp(norm_x_t, x.get_mpf_t(), -exp_part);
  }

  mpf_class norm_x(norm_x_t), nhalf("-0.5", x.get_prec());
  mpf_clear(norm_x_t);
  return exp_part - _mpf_ln1_taylor(norm_x - 1) / _mpf_ln1_taylor(nhalf);
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
  long int_part = mpf_class(floor(x)).get_si();
  x -= int_part;

  mpf_t result;
  mpf_init2(result, x.get_prec());
  if (int_part >= 0) {
    mpf_mul_2exp(result, _mpf_exp_taylor(x * ln2).get_mpf_t(), int_part);
  } else {
    mpf_div_2exp(result, _mpf_exp_taylor(x * ln2).get_mpf_t(), -int_part);
  }
  mpf_class ret(result);
  mpf_clear(result);
  return ret;
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
