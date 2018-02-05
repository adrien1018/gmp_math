#include "gmp_math.h"

mpf_class mpf_pi(mp_bitcnt_t prec)
{
  mpf_class a(1, prec), b(mpf_class(1, prec) / sqrt(mpf_class(2, prec))), t("0.25", prec), p(a);
  while (1) {
    mpf_class A((a + b) / 2), B(sqrt(a * b)), T(t - p * (a - A) * (a - A));
    if (A == B) break;
    a = A, b = B, t = T; p *= 2;
  }
  return (a + b) * (a + b) / 4 / t;
}

mpf_class mpf_sin_taylor(mpf_class x)
{
  mpz_class counter(2);
  mpf_class prv(9, x.get_prec()), result(x), denom(1, x.get_prec()), x2(x * x);
  while (prv != result) {
    prv = result;
    denom *= -counter * (counter + 1);
    counter += 2;
    x *= x2;
    result += x / denom;
  }
  return result;
}

mpf_class mpf_sin(const mpf_class& x)
{
  mpf_class pi2 = 2 * mpf_pi(x.get_prec());
  return mpf_sin_taylor(x - floor(x / pi2) * pi2);
}

mpf_class mpf_cos_taylor(mpf_class x)
{
  mpz_class counter(1);
  mpf_class prv(9, x.get_prec()), result(1, x.get_prec());
  mpf_class denom(1, x.get_prec()), x2(x * x);
  x = 1;
  while (prv != result) {
    prv = result;
    denom *= -counter * (counter + 1);
    counter += 2;
    x *= x2;
    result += x / denom;
  }
  return result;
}

mpf_class mpf_cos(const mpf_class& x)
{
  mpf_class pi2 = 2 * mpf_pi(x.get_prec());
  return mpf_cos_taylor(x - floor(x / pi2) * pi2);
}

mpf_class mpf_arctan_taylor_zero(mpf_class x)
{
  mpz_class denom(1);
  mpf_class prv(9, x.get_prec()), result(x), nx2(-x * x);
  while (prv != result) {
    prv = result;
    x *= nx2;
    denom += 2;
    result += x / denom;
  }
  return result;
}

mpf_class mpf_arctan_taylor_one(mpf_class x)
{
  mpz_class denom_odd(1), denom_even(1);
  mpf_class prv(9, x.get_prec()), result(mpf_pi(x.get_prec()) / 4);
  x -= 1;
  mpf_class x2(x * x), num(x), pow2(2, x.get_prec());
  while (prv != result) {
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
  }
  return result;
}

mpf_class mpf_atan(const mpf_class& x)
{
  const double thresh = 0.7;
  if (x > 1 / thresh)
    return mpf_pi(x.get_prec()) / 2 - mpf_arctan_taylor_zero(1 / x);
  else if (x < -1 / thresh)
    return -mpf_pi(x.get_prec()) / 2 - mpf_arctan_taylor_zero(1 / x);
  else if (x > -thresh && x < thresh)
    return mpf_arctan_taylor_zero(x);
  else if (x < 0)
    return -mpf_arctan_taylor_one(-x);
  else
    return mpf_arctan_taylor_one(x);
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
