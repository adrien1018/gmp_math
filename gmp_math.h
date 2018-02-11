#include <gmp.h>
#include <gmpxx.h>

typedef std::pair<mpz_class, mpz_class> mpz_pair;

mpz_class mpz_sqrt(const mpz_class& x);
mpz_pair mpz_sqrtrem(const mpz_class& x);
mpz_class mpz_root(const mpz_class& x, unsigned long y);
mpz_pair mpz_rootrem(const mpz_class& x, unsigned long y);
bool mpz_perfect_power(const mpz_class& x);
bool mpz_perfect_square(const mpz_class& x);
mpz_class mpz_gcd(const mpz_class& x, const mpz_class& y);
mpz_class mpz_gcdext(const mpz_class& x, const mpz_class& y,
                     mpz_class& a, mpz_class& b);
mpz_class mpz_lcm(const mpz_class& x, const mpz_class& y);
mpz_class mpz_pow(const mpz_class& x, unsigned long y);
mpz_class mpz_powm(const mpz_class& x, const mpz_class& y,
                   const mpz_class& m);

mpf_class mpf_e(mp_bitcnt_t prec);
mpf_class mpf_pi(mp_bitcnt_t prec);
mpf_class mpf_sqrt(const mpf_class& x);
mpf_class mpf_cbrt(const mpf_class& x);
mpf_class mpf_root(const mpf_class& x, unsigned long y);
mpf_class mpf_mulexp2(const mpf_class& x, long y);
mpf_class mpf_sin(const mpf_class& x);
mpf_class mpf_cos(const mpf_class& x);
mpf_class mpf_tan(const mpf_class& x);
mpf_class mpf_asin(const mpf_class& x);
mpf_class mpf_acos(const mpf_class& x);
mpf_class mpf_atan(const mpf_class& x);
mpf_class mpf_atan2(const mpf_class& Y, const mpf_class& X);
long mpf_log2floor(const mpf_class& x);
mpf_class mpf_ln(const mpf_class& x);
mpf_class mpf_log2(const mpf_class& x);
mpf_class mpf_exp(const mpf_class& x);
mpf_class mpf_exp2(const mpf_class& x);
mpf_class mpf_powint(mpf_class x, long y);
mpf_class mpf_pow(const mpf_class& x, const mpf_class& y);
