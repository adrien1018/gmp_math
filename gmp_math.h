#include <gmp.h>
#include <gmpxx.h>
#include <algorithm>

mpf_class mpf_pi(mp_bitcnt_t prec);
mpf_class mpf_sin(const mpf_class& x);
mpf_class mpf_cos(const mpf_class& x);
mpf_class mpf_atan(const mpf_class& x);
mpf_class mpf_atan2(const mpf_class& Y, const mpf_class& X);
