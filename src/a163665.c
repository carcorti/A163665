#define _POSIX_C_SOURCE 200809L

/*
 * OEIS A163665 official production code.
 *
 * Search variable: n = pi(k).  A hit satisfies phi(sigma(n)) == n, and
 * the reported OEIS term is k = prime(n).
 *
 * This official build intentionally implements only bounded segmented execution,
 * deterministic self-tests, and simple manifest summarization.  It is not a
 * campaign driver.
 */

#include <errno.h>
#include <fcntl.h>
#include <inttypes.h>
#include <limits.h>
#include <locale.h>
#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define A163665_VERSION "official"
#define A163665_SEQ_ID "A163665"
#define A163665_MAX_SEGMENT_WIDTH UINT64_C(10000000)
#define A163665_MAX_SIEVE_LIMIT UINT64_C(100000000)
#define A163665_LEHMER_SIEVE UINT32_C(5000000)
#define A163665_PHI_CACHE_X UINT32_C(100000)
#define A163665_PHI_CACHE_S UINT32_C(100)
#define A163665_FACTOR_CAP 96U
#define A163665_STR_CAP 512U
#define A163665_INTERNAL_NTH_PRIME_MAX UINT64_C(1000000000)

#ifndef A163665_SOURCE_PATH
#define A163665_SOURCE_PATH "src/a163665.c"
#endif

#ifndef A163665_BUILD_SOURCE_SHA256
#define A163665_BUILD_SOURCE_SHA256 ""
#endif

typedef struct {
    uint64_t p;
    uint32_t e;
} factor_t;

typedef struct {
    factor_t f[A163665_FACTOR_CAP];
    size_t len;
} factor_list_t;

typedef struct {
    uint64_t mr_count;
    uint64_t rho_count;
    bool used_mr;
    bool used_rho;
} factor_ctx_t;

typedef struct {
    uint64_t raw_width;
    uint64_t candidates_tested;
    uint64_t odd_skipped;
    uint64_t rejected_qminus1;
    uint64_t mr_count;
    uint64_t rho_count;
    uint64_t trial_only_count;
    uint64_t mr_path_count;
    uint64_t rho_path_count;
    uint64_t unresolved_count;
    uint64_t failure_count;
    uint64_t hit_count;
    uint64_t lpf_n_max;
    uint64_t lpf_sigma_max;
    uint64_t lpf_n_le_1e3;
    uint64_t lpf_n_le_1e6;
    uint64_t lpf_n_le_1e9;
    uint64_t lpf_n_gt_1e9;
    uint64_t lpf_sigma_le_1e3;
    uint64_t lpf_sigma_le_1e6;
    uint64_t lpf_sigma_le_1e9;
    uint64_t lpf_sigma_gt_1e9;
} telemetry_t;

typedef struct {
    uint64_t n;
    uint64_t k;
    uint64_t sigma_n;
    uint64_t phi_sigma_n;
    char factor_n[A163665_STR_CAP];
    char factor_sigma[A163665_STR_CAP];
} hit_t;

typedef struct {
    uint32_t *primes;
    uint32_t *pi;
    size_t count;
    uint32_t limit;
} prime_table_t;

typedef struct {
    uint64_t lo;
    uint64_t hi;
    uint64_t raw_width;
    uint64_t candidates_tested;
    uint64_t odd_skipped;
    uint64_t unresolved_count;
    uint64_t failure_count;
    uint64_t hit_count;
    uint64_t binary_fnv64;
    uint64_t source_fnv64;
    uint64_t hits_fnv64;
    char binary_sha256[65];
    char source_sha256[65];
    char hits_sha256[65];
    int status_class;
} manifest_summary_t;

typedef struct {
    uint32_t state[8];
    uint64_t bitlen;
    unsigned char data[64];
    size_t datalen;
} sha256_ctx_t;

static prime_table_t g_lehmer = {0};
static uint32_t g_phi_cache[A163665_PHI_CACHE_S][A163665_PHI_CACHE_X];
static bool g_phi_cache_ready = false;

static const uint64_t known_n[20] = {
    UINT64_C(1), UINT64_C(2), UINT64_C(8), UINT64_C(12),
    UINT64_C(128), UINT64_C(240), UINT64_C(720), UINT64_C(6912),
    UINT64_C(32768), UINT64_C(142560), UINT64_C(712800),
    UINT64_C(1140480), UINT64_C(1190400), UINT64_C(3345408),
    UINT64_C(3571200), UINT64_C(5702400), UINT64_C(14859936),
    UINT64_C(29719872), UINT64_C(50319360), UINT64_C(118879488)
};

static const uint64_t known_k[20] = {
    UINT64_C(2), UINT64_C(3), UINT64_C(19), UINT64_C(37),
    UINT64_C(719), UINT64_C(1511), UINT64_C(5443), UINT64_C(69709),
    UINT64_C(386093), UINT64_C(1907819), UINT64_C(10777931),
    UINT64_C(17819101), UINT64_C(18653749), UINT64_C(56125547),
    UINT64_C(60163267), UINT64_C(98911811), UINT64_C(272887613),
    UINT64_C(567611663), UINT64_C(989060309), UINT64_C(2444540149)
};

static const char manifest_header[] =
    "seq_id\tversion\tstatus\tlo\thi\tparity_policy\tthreads\tsegment_size"
    "\traw_range_width\tcandidates_tested\todd_skipped\trejected_qminus1"
    "\tmiller_rabin_count\tpollard_rho_count\ttrial_only_count\tmr_path_count"
    "\trho_path_count\tlpf_n_max\tlpf_sigma_max\tlpf_n_le_1e3\tlpf_n_le_1e6"
    "\tlpf_n_le_1e9\tlpf_n_gt_1e9\tlpf_sigma_le_1e3\tlpf_sigma_le_1e6"
    "\tlpf_sigma_le_1e9\tlpf_sigma_gt_1e9\tunresolved_count\tfailure_count"
    "\thit_count\twall_seconds\tcandidates_per_second\tbinary_fnv64"
    "\tsource_fnv64\thits_fnv64\tbinary_sha256\tsource_sha256\thits_sha256"
    "\tstarted_unix\tended_unix\n";

static bool line_has_cr(const char *s);

static void usage(FILE *out)
{
    fprintf(out,
            "Usage:\n"
            "  a163665 selftest\n"
            "  a163665 segment --lo N --hi M --threads T --out PREFIX\n"
            "  a163665 summarize MANIFEST.tsv [MANIFEST.tsv ...]\n");
}

static bool parse_u64(const char *s, uint64_t *out)
{
    char *end = NULL;
    unsigned long long v = 0ULL;

    if (s == NULL || s[0] < '0' || s[0] > '9') {
        return false;
    }
    errno = 0;
    v = strtoull(s, &end, 10);
    if (errno != 0 || end == s || *end != '\0') {
        return false;
    }
#if ULLONG_MAX > UINT64_MAX
    if (v > (unsigned long long)UINT64_MAX) {
        return false;
    }
#endif
    *out = (uint64_t)v;
    return true;
}

static bool parse_nonnegative_decimal(const char *s, double *out)
{
    bool saw_digit = false;
    bool saw_dot = false;
    char *end = NULL;
    double v = 0.0;

    if (s == NULL || s[0] == '\0') {
        return false;
    }
    for (const char *p = s; *p != '\0'; p++) {
        if (*p >= '0' && *p <= '9') {
            saw_digit = true;
        } else if (*p == '.' && !saw_dot) {
            saw_dot = true;
        } else {
            return false;
        }
    }
    if (!saw_digit) {
        return false;
    }
    errno = 0;
    v = strtod(s, &end);
    if (errno != 0 || end == s || *end != '\0' || !isfinite(v) || v < 0.0) {
        return false;
    }
    *out = v;
    return true;
}

static bool first_even_in_range(uint64_t lo, uint64_t hi, uint64_t *first_even)
{
    if (lo <= 2U) {
        *first_even = 2U;
    } else if (lo == UINT64_MAX) {
        *first_even = 0U;
        return false;
    } else if ((lo & 1U) == 0U) {
        *first_even = lo;
    } else {
        *first_even = lo + 1U;
    }
    return *first_even <= hi;
}

static uint64_t expected_candidate_count(uint64_t lo, uint64_t hi)
{
    uint64_t first_even = 0;
    uint64_t count = 0;

    if (lo == 1U) {
        count++;
    }
    if (first_even_in_range(lo, hi, &first_even)) {
        count += ((hi - first_even) / 2U) + 1U;
    }
    return count;
}

static double monotonic_seconds(void)
{
    struct timespec ts;

    if (clock_gettime(CLOCK_MONOTONIC, &ts) != 0) {
        return 0.0;
    }
    return (double)ts.tv_sec + ((double)ts.tv_nsec / 1000000000.0);
}

static bool u64_add_checked(uint64_t a, uint64_t b, uint64_t *out)
{
    if (UINT64_MAX - a < b) {
        return false;
    }
    *out = a + b;
    return true;
}

static bool u64_mul_checked(uint64_t a, uint64_t b, uint64_t *out)
{
    __uint128_t z = ((__uint128_t)a) * ((__uint128_t)b);
    if (z > (__uint128_t)UINT64_MAX) {
        return false;
    }
    *out = (uint64_t)z;
    return true;
}

static uint64_t gcd_u64(uint64_t a, uint64_t b)
{
    while (b != 0U) {
        uint64_t t = a % b;
        a = b;
        b = t;
    }
    return a;
}

static uint64_t isqrt_u64(uint64_t n)
{
    uint64_t lo = 0;
    uint64_t hi = UINT64_C(4294967296);

    while (lo + 1U < hi) {
        uint64_t mid = lo + ((hi - lo) / 2U);
        __uint128_t sq = ((__uint128_t)mid) * ((__uint128_t)mid);
        if (sq <= (__uint128_t)n) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    return lo;
}

static uint64_t icbrt_u64(uint64_t n)
{
    uint64_t lo = 0;
    uint64_t hi = UINT64_C(2642246);

    while (lo + 1U < hi) {
        uint64_t mid = lo + ((hi - lo) / 2U);
        __uint128_t cube = ((__uint128_t)mid) * ((__uint128_t)mid) * ((__uint128_t)mid);
        if (cube <= (__uint128_t)n) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    return lo;
}

static uint64_t mod_mul_u64(uint64_t a, uint64_t b, uint64_t m)
{
    return (uint64_t)((((__uint128_t)a) * ((__uint128_t)b)) % ((__uint128_t)m));
}

static uint64_t mod_add_u64(uint64_t a, uint64_t b, uint64_t m)
{
    return (uint64_t)((((__uint128_t)a) + ((__uint128_t)b)) % ((__uint128_t)m));
}

static uint64_t mod_pow_u64(uint64_t a, uint64_t e, uint64_t m)
{
    uint64_t r = 1U % m;
    uint64_t x = a % m;

    while (e != 0U) {
        if ((e & 1U) != 0U) {
            r = mod_mul_u64(r, x, m);
        }
        e >>= 1U;
        if (e != 0U) {
            x = mod_mul_u64(x, x, m);
        }
    }
    return r;
}

static bool is_prime64(uint64_t n, factor_ctx_t *ctx)
{
    static const uint64_t bases[] = {
        UINT64_C(2), UINT64_C(325), UINT64_C(9375), UINT64_C(28178),
        UINT64_C(450775), UINT64_C(9780504), UINT64_C(1795265022)
    };
    uint64_t d = 0;
    unsigned int s = 0;

    if (ctx != NULL) {
        ctx->mr_count++;
        ctx->used_mr = true;
    }
    if (n < 2U) {
        return false;
    }
    if ((n % 2U) == 0U) {
        return n == 2U;
    }
    if ((n % 3U) == 0U) {
        return n == 3U;
    }

    d = n - 1U;
    while ((d & 1U) == 0U) {
        d >>= 1U;
        s++;
    }

    for (size_t i = 0; i < (sizeof(bases) / sizeof(bases[0])); i++) {
        uint64_t a = bases[i] % n;
        uint64_t x = 0;
        bool witness_ok = false;

        if (a == 0U) {
            continue;
        }
        x = mod_pow_u64(a, d, n);
        if (x == 1U || x == n - 1U) {
            continue;
        }
        for (unsigned int r = 1; r < s; r++) {
            x = mod_mul_u64(x, x, n);
            if (x == n - 1U) {
                witness_ok = true;
                break;
            }
        }
        if (!witness_ok) {
            return false;
        }
    }
    return true;
}

static uint64_t pollard_rho_one(uint64_t n, factor_ctx_t *ctx)
{
    if ((n % 2U) == 0U) {
        return 2U;
    }
    if ((n % 3U) == 0U) {
        return 3U;
    }
    if (ctx != NULL) {
        ctx->rho_count++;
        ctx->used_rho = true;
    }
    for (uint64_t c = 1; c <= 32U; c++) {
        uint64_t x = 2U + c;
        uint64_t y = x;

        for (uint32_t iter = 0; iter < UINT32_C(200000); iter++) {
            uint64_t d = 0;
            x = mod_add_u64(mod_mul_u64(x, x, n), c, n);
            y = mod_add_u64(mod_mul_u64(y, y, n), c, n);
            y = mod_add_u64(mod_mul_u64(y, y, n), c, n);
            d = gcd_u64((x > y) ? (x - y) : (y - x), n);
            if (d == n) {
                break;
            }
            if (d > 1U) {
                return d;
            }
        }
    }
    return 0U;
}

static bool factor_add(factor_list_t *list, uint64_t p, uint32_t e)
{
    for (size_t i = 0; i < list->len; i++) {
        if (list->f[i].p == p) {
            if (UINT32_MAX - list->f[i].e < e) {
                return false;
            }
            list->f[i].e += e;
            return true;
        }
    }
    if (list->len >= A163665_FACTOR_CAP) {
        return false;
    }
    list->f[list->len].p = p;
    list->f[list->len].e = e;
    list->len++;
    return true;
}

static bool factor_rec(uint64_t n, factor_list_t *out, factor_ctx_t *ctx)
{
    uint64_t d = 0;

    if (n == 1U) {
        return true;
    }
    if (is_prime64(n, ctx)) {
        return factor_add(out, n, 1U);
    }
    d = pollard_rho_one(n, ctx);
    if (d == 0U || d == n) {
        return false;
    }
    return factor_rec(d, out, ctx) && factor_rec(n / d, out, ctx);
}

static bool factor_u64(uint64_t n, factor_list_t *out, factor_ctx_t *ctx)
{
    static const uint32_t small_primes[] = {
        2U, 3U, 5U, 7U, 11U, 13U, 17U, 19U, 23U, 29U, 31U, 37U,
        41U, 43U, 47U, 53U, 59U, 61U, 67U, 71U, 73U, 79U, 83U, 89U, 97U
    };

    out->len = 0;
    if (n == 1U) {
        return true;
    }
    for (size_t i = 0; i < (sizeof(small_primes) / sizeof(small_primes[0])); i++) {
        uint32_t p = small_primes[i];
        if ((n % (uint64_t)p) == 0U) {
            uint32_t e = 0;
            do {
                n /= (uint64_t)p;
                e++;
            } while ((n % (uint64_t)p) == 0U);
            if (!factor_add(out, (uint64_t)p, e)) {
                return false;
            }
        }
    }
    return factor_rec(n, out, ctx);
}

static int factor_cmp(const void *a, const void *b)
{
    const factor_t *fa = (const factor_t *)a;
    const factor_t *fb = (const factor_t *)b;

    if (fa->p < fb->p) {
        return -1;
    }
    if (fa->p > fb->p) {
        return 1;
    }
    return 0;
}

static void factor_sort(factor_list_t *list)
{
    qsort(list->f, list->len, sizeof(list->f[0]), factor_cmp);
}

static uint64_t factor_largest_prime(const factor_list_t *list)
{
    uint64_t best = 1U;

    for (size_t i = 0; i < list->len; i++) {
        if (list->f[i].p > best) {
            best = list->f[i].p;
        }
    }
    return best;
}

static bool factor_to_string(const factor_list_t *list, char *buf, size_t cap)
{
    size_t used = 0;

    if (cap == 0U) {
        return false;
    }
    if (list->len == 0U) {
        int w = snprintf(buf, cap, "1");
        return w > 0 && (size_t)w < cap;
    }
    buf[0] = '\0';
    for (size_t i = 0; i < list->len; i++) {
        int w = 0;
        if (i == 0U) {
            if (list->f[i].e == 1U) {
                w = snprintf(buf + used, cap - used, "%" PRIu64, list->f[i].p);
            } else {
                w = snprintf(buf + used, cap - used, "%" PRIu64 "^%" PRIu32,
                             list->f[i].p, list->f[i].e);
            }
        } else {
            if (list->f[i].e == 1U) {
                w = snprintf(buf + used, cap - used, "*%" PRIu64, list->f[i].p);
            } else {
                w = snprintf(buf + used, cap - used, "*%" PRIu64 "^%" PRIu32,
                             list->f[i].p, list->f[i].e);
            }
        }
        if (w <= 0 || (size_t)w >= cap - used) {
            return false;
        }
        used += (size_t)w;
    }
    return true;
}

static bool phi_from_factorization(uint64_t n, const factor_list_t *fac, uint64_t *phi)
{
    uint64_t result = n;

    for (size_t i = 0; i < fac->len; i++) {
        uint64_t divided = result / fac->f[i].p;
        if (!u64_mul_checked(divided, fac->f[i].p - 1U, &result)) {
            return false;
        }
    }
    *phi = result;
    return true;
}

static bool sigma_from_factorization(const factor_list_t *fac, uint64_t *sigma)
{
    uint64_t result = 1U;

    for (size_t i = 0; i < fac->len; i++) {
        uint64_t p = fac->f[i].p;
        uint64_t term = 1U;
        uint64_t pow = 1U;
        for (uint32_t e = 0; e < fac->f[i].e; e++) {
            uint64_t new_pow = 0;
            if (!u64_mul_checked(pow, p, &new_pow)) {
                return false;
            }
            pow = new_pow;
            if (!u64_add_checked(term, pow, &term)) {
                return false;
            }
        }
        if (!u64_mul_checked(result, term, &result)) {
            return false;
        }
    }
    *sigma = result;
    return true;
}

static void profile_lpf(uint64_t lpf, bool sigma_profile, telemetry_t *tel)
{
    if (sigma_profile) {
        if (lpf > tel->lpf_sigma_max) {
            tel->lpf_sigma_max = lpf;
        }
        if (lpf <= UINT64_C(1000)) {
            tel->lpf_sigma_le_1e3++;
        } else if (lpf <= UINT64_C(1000000)) {
            tel->lpf_sigma_le_1e6++;
        } else if (lpf <= UINT64_C(1000000000)) {
            tel->lpf_sigma_le_1e9++;
        } else {
            tel->lpf_sigma_gt_1e9++;
        }
    } else {
        if (lpf > tel->lpf_n_max) {
            tel->lpf_n_max = lpf;
        }
        if (lpf <= UINT64_C(1000)) {
            tel->lpf_n_le_1e3++;
        } else if (lpf <= UINT64_C(1000000)) {
            tel->lpf_n_le_1e6++;
        } else if (lpf <= UINT64_C(1000000000)) {
            tel->lpf_n_le_1e9++;
        } else {
            tel->lpf_n_gt_1e9++;
        }
    }
}

static uint32_t rotr32(uint32_t x, unsigned int n)
{
    if (n == 0U) {
        return x;
    }
    return (x >> n) | (x << (32U - n));
}

static void sha256_transform(sha256_ctx_t *ctx, const unsigned char data[64])
{
    static const uint32_t k[64] = {
        UINT32_C(0x428a2f98), UINT32_C(0x71374491), UINT32_C(0xb5c0fbcf), UINT32_C(0xe9b5dba5),
        UINT32_C(0x3956c25b), UINT32_C(0x59f111f1), UINT32_C(0x923f82a4), UINT32_C(0xab1c5ed5),
        UINT32_C(0xd807aa98), UINT32_C(0x12835b01), UINT32_C(0x243185be), UINT32_C(0x550c7dc3),
        UINT32_C(0x72be5d74), UINT32_C(0x80deb1fe), UINT32_C(0x9bdc06a7), UINT32_C(0xc19bf174),
        UINT32_C(0xe49b69c1), UINT32_C(0xefbe4786), UINT32_C(0x0fc19dc6), UINT32_C(0x240ca1cc),
        UINT32_C(0x2de92c6f), UINT32_C(0x4a7484aa), UINT32_C(0x5cb0a9dc), UINT32_C(0x76f988da),
        UINT32_C(0x983e5152), UINT32_C(0xa831c66d), UINT32_C(0xb00327c8), UINT32_C(0xbf597fc7),
        UINT32_C(0xc6e00bf3), UINT32_C(0xd5a79147), UINT32_C(0x06ca6351), UINT32_C(0x14292967),
        UINT32_C(0x27b70a85), UINT32_C(0x2e1b2138), UINT32_C(0x4d2c6dfc), UINT32_C(0x53380d13),
        UINT32_C(0x650a7354), UINT32_C(0x766a0abb), UINT32_C(0x81c2c92e), UINT32_C(0x92722c85),
        UINT32_C(0xa2bfe8a1), UINT32_C(0xa81a664b), UINT32_C(0xc24b8b70), UINT32_C(0xc76c51a3),
        UINT32_C(0xd192e819), UINT32_C(0xd6990624), UINT32_C(0xf40e3585), UINT32_C(0x106aa070),
        UINT32_C(0x19a4c116), UINT32_C(0x1e376c08), UINT32_C(0x2748774c), UINT32_C(0x34b0bcb5),
        UINT32_C(0x391c0cb3), UINT32_C(0x4ed8aa4a), UINT32_C(0x5b9cca4f), UINT32_C(0x682e6ff3),
        UINT32_C(0x748f82ee), UINT32_C(0x78a5636f), UINT32_C(0x84c87814), UINT32_C(0x8cc70208),
        UINT32_C(0x90befffa), UINT32_C(0xa4506ceb), UINT32_C(0xbef9a3f7), UINT32_C(0xc67178f2)
    };
    uint32_t m[64];
    uint32_t a = 0;
    uint32_t b = 0;
    uint32_t c = 0;
    uint32_t d = 0;
    uint32_t e = 0;
    uint32_t f = 0;
    uint32_t g = 0;
    uint32_t h = 0;

    for (size_t i = 0; i < 16U; i++) {
        size_t j = i * 4U;
        m[i] = ((uint32_t)data[j] << 24U) | ((uint32_t)data[j + 1U] << 16U) |
               ((uint32_t)data[j + 2U] << 8U) | (uint32_t)data[j + 3U];
    }
    for (size_t i = 16U; i < 64U; i++) {
        uint32_t s0 = rotr32(m[i - 15U], 7U) ^ rotr32(m[i - 15U], 18U) ^ (m[i - 15U] >> 3U);
        uint32_t s1 = rotr32(m[i - 2U], 17U) ^ rotr32(m[i - 2U], 19U) ^ (m[i - 2U] >> 10U);
        m[i] = m[i - 16U] + s0 + m[i - 7U] + s1;
    }

    a = ctx->state[0];
    b = ctx->state[1];
    c = ctx->state[2];
    d = ctx->state[3];
    e = ctx->state[4];
    f = ctx->state[5];
    g = ctx->state[6];
    h = ctx->state[7];

    for (size_t i = 0; i < 64U; i++) {
        uint32_t s1 = rotr32(e, 6U) ^ rotr32(e, 11U) ^ rotr32(e, 25U);
        uint32_t ch = (e & f) ^ ((~e) & g);
        uint32_t temp1 = h + s1 + ch + k[i] + m[i];
        uint32_t s0 = rotr32(a, 2U) ^ rotr32(a, 13U) ^ rotr32(a, 22U);
        uint32_t maj = (a & b) ^ (a & c) ^ (b & c);
        uint32_t temp2 = s0 + maj;
        h = g;
        g = f;
        f = e;
        e = d + temp1;
        d = c;
        c = b;
        b = a;
        a = temp1 + temp2;
    }

    ctx->state[0] += a;
    ctx->state[1] += b;
    ctx->state[2] += c;
    ctx->state[3] += d;
    ctx->state[4] += e;
    ctx->state[5] += f;
    ctx->state[6] += g;
    ctx->state[7] += h;
}

static void sha256_init(sha256_ctx_t *ctx)
{
    ctx->datalen = 0;
    ctx->bitlen = 0;
    ctx->state[0] = UINT32_C(0x6a09e667);
    ctx->state[1] = UINT32_C(0xbb67ae85);
    ctx->state[2] = UINT32_C(0x3c6ef372);
    ctx->state[3] = UINT32_C(0xa54ff53a);
    ctx->state[4] = UINT32_C(0x510e527f);
    ctx->state[5] = UINT32_C(0x9b05688c);
    ctx->state[6] = UINT32_C(0x1f83d9ab);
    ctx->state[7] = UINT32_C(0x5be0cd19);
}

static void sha256_update(sha256_ctx_t *ctx, const unsigned char *data, size_t len)
{
    for (size_t i = 0; i < len; i++) {
        ctx->data[ctx->datalen++] = data[i];
        if (ctx->datalen == 64U) {
            sha256_transform(ctx, ctx->data);
            ctx->bitlen += 512U;
            ctx->datalen = 0;
        }
    }
}

static void sha256_final(sha256_ctx_t *ctx, unsigned char hash[32])
{
    size_t i = ctx->datalen;

    ctx->data[i++] = 0x80U;
    if (i > 56U) {
        while (i < 64U) {
            ctx->data[i++] = 0U;
        }
        sha256_transform(ctx, ctx->data);
        i = 0;
    }
    while (i < 56U) {
        ctx->data[i++] = 0U;
    }

    ctx->bitlen += (uint64_t)ctx->datalen * 8U;
    for (unsigned int j = 0; j < 8U; j++) {
        ctx->data[63U - j] = (unsigned char)(ctx->bitlen >> (8U * j));
    }
    sha256_transform(ctx, ctx->data);

    for (size_t j = 0; j < 8U; j++) {
        hash[j * 4U] = (unsigned char)(ctx->state[j] >> 24U);
        hash[j * 4U + 1U] = (unsigned char)(ctx->state[j] >> 16U);
        hash[j * 4U + 2U] = (unsigned char)(ctx->state[j] >> 8U);
        hash[j * 4U + 3U] = (unsigned char)ctx->state[j];
    }
}

static void bytes_to_hex(const unsigned char *bytes, size_t len, char *hex, size_t cap)
{
    static const char digits[] = "0123456789abcdef";

    if (cap < (len * 2U) + 1U) {
        return;
    }
    for (size_t i = 0; i < len; i++) {
        hex[i * 2U] = digits[bytes[i] >> 4U];
        hex[i * 2U + 1U] = digits[bytes[i] & 0x0fU];
    }
    hex[len * 2U] = '\0';
}

static bool sha256_file_hex(const char *path, char hex[65])
{
    FILE *f = fopen(path, "rb");
    sha256_ctx_t ctx;
    unsigned char buf[8192];
    unsigned char hash[32];

    if (f == NULL) {
        hex[0] = '\0';
        return false;
    }
    sha256_init(&ctx);
    for (;;) {
        size_t n = fread(buf, 1U, sizeof(buf), f);
        if (n > 0U) {
            sha256_update(&ctx, buf, n);
        }
        if (n < sizeof(buf)) {
            if (ferror(f) != 0) {
                fclose(f);
                hex[0] = '\0';
                return false;
            }
            break;
        }
    }
    fclose(f);
    sha256_final(&ctx, hash);
    bytes_to_hex(hash, sizeof(hash), hex, 65U);
    return true;
}

static uint64_t fnv1a_file64(const char *path, bool *ok);

static bool resolve_source_checksums(uint64_t *source_sum, char source_sha256[65])
{
    const char *source_paths[] = {
        A163665_SOURCE_PATH,
        "src/a163665.c",
        "a163665.c"
    };

    for (size_t i = 0; i < (sizeof(source_paths) / sizeof(source_paths[0])); i++) {
        bool source_ok = false;
        uint64_t sum = fnv1a_file64(source_paths[i], &source_ok);
        if (source_ok && sha256_file_hex(source_paths[i], source_sha256)) {
            *source_sum = sum;
            return true;
        }
    }
    source_sha256[0] = '\0';
    *source_sum = 0U;
    return false;
}

static bool is_sha256_hex(const char *s)
{
    bool any_nonzero = false;

    if (s == NULL || strlen(s) != 64U) {
        return false;
    }
    for (size_t i = 0; i < 64U; i++) {
        char c = s[i];
        bool hex = (c >= '0' && c <= '9') || (c >= 'a' && c <= 'f');
        if (!hex) {
            return false;
        }
        if (c != '0') {
            any_nonzero = true;
        }
    }
    return any_nonzero;
}

static uint64_t fnv1a_file64(const char *path, bool *ok)
{
    FILE *f = fopen(path, "rb");
    uint64_t h = UINT64_C(14695981039346656037);
    unsigned char buf[8192];

    if (f == NULL) {
        *ok = false;
        return 0U;
    }
    for (;;) {
        size_t n = fread(buf, 1U, sizeof(buf), f);
        for (size_t i = 0; i < n; i++) {
            h ^= (uint64_t)buf[i];
            h *= UINT64_C(1099511628211);
        }
        if (n < sizeof(buf)) {
            if (ferror(f) != 0) {
                fclose(f);
                *ok = false;
                return 0U;
            }
            break;
        }
    }
    fclose(f);
    *ok = true;
    return h;
}

static bool fsync_parent_dir(const char *path)
{
    char dir_path[4096];
    const char *slash = strrchr(path, '/');
    size_t len;
    int flags = O_RDONLY;
    int fd = -1;
    bool ok = false;

    if (slash == NULL) {
        if (snprintf(dir_path, sizeof(dir_path), ".") >= (int)sizeof(dir_path)) {
            return false;
        }
    } else if (slash == path) {
        if (snprintf(dir_path, sizeof(dir_path), "/") >= (int)sizeof(dir_path)) {
            return false;
        }
    } else {
        len = (size_t)(slash - path);
        if (len + 1U > sizeof(dir_path)) {
            return false;
        }
        memcpy(dir_path, path, len);
        dir_path[len] = '\0';
    }
#ifdef O_DIRECTORY
    flags |= O_DIRECTORY;
#endif
    fd = open(dir_path, flags);
    if (fd < 0) {
        return false;
    }
    ok = (fsync(fd) == 0);
    if (close(fd) != 0) {
        ok = false;
    }
    return ok;
}

static bool publish_file_noreplace(const char *tmp_path, const char *final_path)
{
    if (link(tmp_path, final_path) != 0) {
        fprintf(stderr, "publish(%s,%s) failed: %s\n", tmp_path, final_path, strerror(errno));
        return false;
    }
    if (!fsync_parent_dir(final_path)) {
        fprintf(stderr, "fsync parent for %s failed: %s\n", final_path, strerror(errno));
        return false;
    }
    if (unlink(tmp_path) != 0) {
        fprintf(stderr, "warning: could not remove temporary file %s: %s\n",
                tmp_path, strerror(errno));
    } else if (!fsync_parent_dir(tmp_path)) {
        fprintf(stderr, "warning: fsync parent after removing %s failed: %s\n",
                tmp_path, strerror(errno));
    }
    return true;
}

static FILE *open_exclusive_temp(const char *tmp_path)
{
    int flags = O_WRONLY | O_CREAT | O_EXCL;
    int fd = -1;
    FILE *out = NULL;

#ifdef O_NOFOLLOW
    flags |= O_NOFOLLOW;
#endif
    fd = open(tmp_path, flags, S_IRUSR | S_IWUSR);
    if (fd < 0) {
        return NULL;
    }
    out = fdopen(fd, "wb");
    if (out == NULL) {
        close(fd);
        unlink(tmp_path);
        return NULL;
    }
    return out;
}

static bool output_paths_available(const char *prefix)
{
    char path[4096];

    if (snprintf(path, sizeof(path), "%s.hits.tsv", prefix) >= (int)sizeof(path) ||
        access(path, F_OK) == 0) {
        return false;
    }
    if (errno != ENOENT) {
        return false;
    }
    if (snprintf(path, sizeof(path), "%s.manifest.tsv", prefix) >= (int)sizeof(path) ||
        access(path, F_OK) == 0) {
        return false;
    }
    if (errno != ENOENT) {
        return false;
    }
    if (snprintf(path, sizeof(path), "%s.hits.tsv.tmp", prefix) >= (int)sizeof(path) ||
        access(path, F_OK) == 0) {
        return false;
    }
    if (errno != ENOENT) {
        return false;
    }
    if (snprintf(path, sizeof(path), "%s.manifest.tsv.tmp", prefix) >= (int)sizeof(path) ||
        access(path, F_OK) == 0) {
        return false;
    }
    return errno == ENOENT;
}

static bool segment_domain_supported(uint64_t lo, uint64_t hi, char *reason, size_t reason_cap)
{
    uint64_t raw_width = 0;
    uint64_t first_even = 0;
    bool has_candidate = false;

    if (lo == 0U) {
        snprintf(reason, reason_cap, "lo must be >= 1");
        return false;
    }
    if (hi < lo) {
        snprintf(reason, reason_cap, "hi must be >= lo");
        return false;
    }
    raw_width = hi - lo + 1U;
    if (raw_width > A163665_MAX_SEGMENT_WIDTH) {
        snprintf(reason, reason_cap, "segment width %" PRIu64 " exceeds production limit %" PRIu64,
                 raw_width, (uint64_t)A163665_MAX_SEGMENT_WIDTH);
        return false;
    }
    has_candidate = (lo == 1U) || first_even_in_range(lo, hi, &first_even);
    if (has_candidate && isqrt_u64(hi) > A163665_MAX_SIEVE_LIMIT) {
        snprintf(reason, reason_cap,
                 "segment has retained candidates but sqrt(hi) exceeds production sieve limit %" PRIu64,
                 (uint64_t)A163665_MAX_SIEVE_LIMIT);
        return false;
    }
    if (reason_cap > 0U) {
        reason[0] = '\0';
    }
    return true;
}

static bool lpf_bucket_max_consistent(uint64_t max_lpf, uint64_t le_1e3,
                                      uint64_t le_1e6, uint64_t le_1e9,
                                      uint64_t gt_1e9, uint64_t item_count)
{
    if (item_count == 0U) {
        return max_lpf == 0U && le_1e3 == 0U && le_1e6 == 0U &&
               le_1e9 == 0U && gt_1e9 == 0U;
    }
    if (max_lpf == 0U) {
        return false;
    }
    if (max_lpf <= UINT64_C(1000)) {
        return le_1e3 > 0U && le_1e6 == 0U && le_1e9 == 0U && gt_1e9 == 0U;
    }
    if (max_lpf <= UINT64_C(1000000)) {
        return le_1e6 > 0U && le_1e9 == 0U && gt_1e9 == 0U;
    }
    if (max_lpf <= UINT64_C(1000000000)) {
        return le_1e9 > 0U && gt_1e9 == 0U;
    }
    return gt_1e9 > 0U;
}

static bool build_prime_table(uint64_t limit, uint32_t **primes_out, size_t *count_out)
{
    bool *is_composite = NULL;
    uint32_t *primes = NULL;
    size_t count = 0;
    size_t cap = 0;
    size_t n = 0;

    if (limit > A163665_MAX_SIEVE_LIMIT || limit > (uint64_t)SIZE_MAX - 1U) {
        return false;
    }
    n = (size_t)limit + 1U;
    is_composite = (bool *)calloc(n, sizeof(*is_composite));
    if (is_composite == NULL) {
        return false;
    }
    cap = (n / 8U) + 16U;
    primes = (uint32_t *)malloc(cap * sizeof(*primes));
    if (primes == NULL) {
        free(is_composite);
        return false;
    }
    for (size_t i = 2; i < n; i++) {
        if (!is_composite[i]) {
            if (count == cap) {
                size_t new_cap = cap * 2U;
                uint32_t *new_primes = (uint32_t *)realloc(primes, new_cap * sizeof(*primes));
                if (new_primes == NULL) {
                    free(primes);
                    free(is_composite);
                    return false;
                }
                primes = new_primes;
                cap = new_cap;
            }
            primes[count++] = (uint32_t)i;
            if (i <= (n - 1U) / i) {
                for (size_t j = i * i; j < n; j += i) {
                    is_composite[j] = true;
                }
            }
        }
    }
    free(is_composite);
    *primes_out = primes;
    *count_out = count;
    return true;
}

static bool init_lehmer_table(void)
{
    bool *is_composite = NULL;
    uint32_t *primes = NULL;
    uint32_t *pi = NULL;
    size_t count = 0;
    uint32_t limit = A163665_LEHMER_SIEVE;

    if (g_lehmer.primes != NULL) {
        return true;
    }
    is_composite = (bool *)calloc((size_t)limit + 1U, sizeof(*is_composite));
    primes = (uint32_t *)malloc(((size_t)limit / 8U + 16U) * sizeof(*primes));
    pi = (uint32_t *)calloc((size_t)limit + 1U, sizeof(*pi));
    if (is_composite == NULL || primes == NULL || pi == NULL) {
        free(is_composite);
        free(primes);
        free(pi);
        return false;
    }
    for (uint32_t i = 2; i <= limit; i++) {
        if (!is_composite[i]) {
            primes[count++] = i;
            if (i <= limit / i) {
                for (uint32_t j = i * i; j <= limit; j += i) {
                    is_composite[j] = true;
                }
            }
        }
        pi[i] = (uint32_t)count;
    }
    free(is_composite);
    g_lehmer.primes = primes;
    g_lehmer.pi = pi;
    g_lehmer.count = count;
    g_lehmer.limit = limit;
    for (uint32_t x = 0; x < A163665_PHI_CACHE_X; x++) {
        g_phi_cache[0][x] = x;
    }
    for (uint32_t s = 1; s < A163665_PHI_CACHE_S; s++) {
        uint32_t p = g_lehmer.primes[s - 1U];
        for (uint32_t x = 0; x < A163665_PHI_CACHE_X; x++) {
            g_phi_cache[s][x] = g_phi_cache[s - 1U][x] - g_phi_cache[s - 1U][x / p];
        }
    }
    g_phi_cache_ready = true;
    return true;
}

static void free_lehmer_table(void)
{
    free(g_lehmer.primes);
    free(g_lehmer.pi);
    g_lehmer.primes = NULL;
    g_lehmer.pi = NULL;
    g_lehmer.count = 0U;
    g_lehmer.limit = 0U;
    g_phi_cache_ready = false;
}

static bool lehmer_pi_domain_ok(uint64_t x)
{
    return x <= (uint64_t)g_lehmer.limit || isqrt_u64(x) <= (uint64_t)g_lehmer.limit;
}

static uint64_t phi_lehmer(uint64_t x, uint64_t s)
{
    if (g_phi_cache_ready && s < A163665_PHI_CACHE_S && x < A163665_PHI_CACHE_X) {
        return (uint64_t)g_phi_cache[s][x];
    }
    if (s == 0U) {
        return x;
    }
    if (s == 1U) {
        return x - (x / 2U);
    }
    return phi_lehmer(x, s - 1U) - phi_lehmer(x / g_lehmer.primes[s - 1U], s - 1U);
}

static uint64_t lehmer_pi(uint64_t x)
{
    uint64_t a = 0;
    uint64_t b = 0;
    uint64_t c = 0;
    uint64_t sum = 0;

    if (g_lehmer.primes == NULL || g_lehmer.pi == NULL || g_lehmer.count == 0U) {
        return 0;
    }
    if (x <= (uint64_t)g_lehmer.limit) {
        return (uint64_t)g_lehmer.pi[x];
    }
    a = lehmer_pi(isqrt_u64(isqrt_u64(x)));
    b = lehmer_pi(isqrt_u64(x));
    c = lehmer_pi(icbrt_u64(x));
    sum = phi_lehmer(x, a) + (((b + a - 2U) * (b - a + 1U)) / 2U);

    for (uint64_t i = a; i < b; i++) {
        uint64_t w = x / g_lehmer.primes[i];
        sum -= lehmer_pi(w);
        if (i < c) {
            uint64_t lim = lehmer_pi(isqrt_u64(w));
            for (uint64_t j = i; j < lim; j++) {
                sum -= lehmer_pi(w / g_lehmer.primes[j]) - j;
            }
        }
    }
    return sum;
}

static bool nth_prime_primecount(uint64_t n, uint64_t *p)
{
    char command[128];
    char line[128];
    FILE *pipe = NULL;
    int status = 0;
    size_t len = 0;

    if (snprintf(command, sizeof(command),
                 "primecount %" PRIu64 " --nth-prime --threads=1", n) >=
        (int)sizeof(command)) {
        return false;
    }
    pipe = popen(command, "r");
    if (pipe == NULL) {
        return false;
    }
    if (fgets(line, sizeof(line), pipe) == NULL) {
        status = pclose(pipe);
        (void)status;
        return false;
    }
    status = pclose(pipe);
    if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
        return false;
    }
    len = strlen(line);
    if (len == 0U || line_has_cr(line)) {
        return false;
    }
    if (line[len - 1U] == '\n') {
        line[len - 1U] = '\0';
    }
    return parse_u64(line, p);
}

static bool nth_prime(uint64_t n, uint64_t *p)
{
    uint64_t lo = 0;
    uint64_t hi = 0;
    long double nd = 0.0L;
    long double estimate = 0.0L;

    if (!init_lehmer_table()) {
        return false;
    }
    if (n == 0U) {
        return false;
    }
    if (n == 1U) {
        *p = 2U;
        return true;
    }
    if (n == 2U) {
        *p = 3U;
        return true;
    }
    if (n > A163665_INTERNAL_NTH_PRIME_MAX) {
        return nth_prime_primecount(n, p);
    }
    nd = (long double)n;
    estimate = nd * (logl(nd) + logl(logl(nd))) + 16.0L;
    if (estimate > (long double)(UINT64_MAX / 2U)) {
        return false;
    }
    hi = (uint64_t)estimate;
    if (hi < 16U) {
        hi = 16U;
    }
    if (!lehmer_pi_domain_ok(hi)) {
        return false;
    }
    while (lehmer_pi(hi) < n) {
        if (hi > UINT64_MAX / 2U) {
            return false;
        }
        hi *= 2U;
        if (!lehmer_pi_domain_ok(hi)) {
            return false;
        }
    }
    while (lo + 1U < hi) {
        uint64_t mid = lo + ((hi - lo) / 2U);
        if (lehmer_pi(mid) >= n) {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    *p = hi;
    return true;
}

static bool direct_phi_sigma(uint64_t n, uint64_t *sigma_n, uint64_t *phi_sigma_n)
{
    factor_list_t nf;
    factor_list_t sf;
    factor_ctx_t ctx = {0};

    if (!factor_u64(n, &nf, &ctx)) {
        return false;
    }
    if (!sigma_from_factorization(&nf, sigma_n)) {
        return false;
    }
    ctx.used_mr = false;
    ctx.used_rho = false;
    if (!factor_u64(*sigma_n, &sf, &ctx)) {
        return false;
    }
    if (!phi_from_factorization(*sigma_n, &sf, phi_sigma_n)) {
        return false;
    }
    return true;
}

static bool segment_build_sigma(uint64_t lo, uint64_t hi, uint64_t *nums,
                                uint64_t *sigma, uint64_t *lpf_n,
                                size_t count, telemetry_t *tel)
{
    uint64_t first_even = 0;
    bool has_even = false;
    uint32_t *primes = NULL;
    size_t prime_count = 0;
    uint64_t limit = isqrt_u64(hi);

    (void)count;
    if (!build_prime_table(limit, &primes, &prime_count)) {
        return false;
    }
    has_even = first_even_in_range(lo, hi, &first_even);

    for (size_t i = 0; nums[i] != 0U; i++) {
        sigma[i] = 1U;
        lpf_n[i] = (nums[i] == 1U) ? 1U : 0U;
    }

    for (size_t pi = 0; pi < prime_count; pi++) {
        uint64_t p = (uint64_t)primes[pi];
        uint64_t r = lo % p;
        uint64_t gap = (r == 0U) ? 0U : p - r;
        uint64_t m = 0;

        if (gap != 0U && lo > UINT64_MAX - gap) {
            continue;
        }
        m = lo + gap;

        while (m <= hi) {
            size_t idx = SIZE_MAX;
            if (m == 1U && lo == 1U) {
                idx = 0U;
            } else if (has_even && m >= first_even && (m & 1U) == 0U) {
                idx = (size_t)((m - first_even) / 2U) + ((lo == 1U) ? 1U : 0U);
            }
            if (idx != SIZE_MAX && nums[idx] != 0U && (nums[idx] % p) == 0U) {
                uint64_t ppow = 1U;
                uint64_t term = 1U;
                while ((nums[idx] % p) == 0U) {
                    uint64_t new_pow = 0;
                    nums[idx] /= p;
                    if (!u64_mul_checked(ppow, p, &new_pow)) {
                        free(primes);
                        return false;
                    }
                    ppow = new_pow;
                    if (!u64_add_checked(term, ppow, &term)) {
                        free(primes);
                        return false;
                    }
                }
                if (!u64_mul_checked(sigma[idx], term, &sigma[idx])) {
                    free(primes);
                    return false;
                }
                if (p > lpf_n[idx]) {
                    lpf_n[idx] = p;
                }
            }
            if (UINT64_MAX - m < p) {
                break;
            }
            m += p;
        }
    }

    for (size_t i = 0; nums[i] != 0U; i++) {
        if (nums[i] > 1U) {
            uint64_t term = nums[i] + 1U;
            if (!u64_mul_checked(sigma[i], term, &sigma[i])) {
                free(primes);
                return false;
            }
            if (nums[i] > lpf_n[i]) {
                lpf_n[i] = nums[i];
            }
        }
        profile_lpf(lpf_n[i], false, tel);
    }

    free(primes);
    return true;
}

static bool prepare_candidates(uint64_t lo, uint64_t hi, uint64_t **nums_out,
                               uint64_t **sigma_out, uint64_t **lpf_out,
                               unsigned char **hit_flags_out, size_t *count_out,
                               telemetry_t *tel)
{
    uint64_t first_even = 0;
    uint64_t even_count = 0;
    uint64_t count64 = 0;
    bool has_even = false;
    size_t count = 0;
    uint64_t *nums = NULL;
    uint64_t *sigma = NULL;
    uint64_t *lpf = NULL;
    unsigned char *hit_flags = NULL;

    if (lo == 0U || hi < lo) {
        return false;
    }
    tel->raw_width = hi - lo + 1U;
    if (tel->raw_width > A163665_MAX_SEGMENT_WIDTH) {
        fprintf(stderr, "segment width %" PRIu64 " exceeds production limit %" PRIu64 "\n",
                tel->raw_width, (uint64_t)A163665_MAX_SEGMENT_WIDTH);
        return false;
    }
    has_even = first_even_in_range(lo, hi, &first_even);
    if (has_even) {
        even_count = ((hi - first_even) / 2U) + 1U;
    }
    count64 = even_count + ((lo == 1U) ? 1U : 0U);
    if (count64 > (uint64_t)(SIZE_MAX - 1U)) {
        return false;
    }
    count = (size_t)count64;
    nums = (uint64_t *)calloc(count + 1U, sizeof(*nums));
    sigma = (uint64_t *)calloc(count + 1U, sizeof(*sigma));
    lpf = (uint64_t *)calloc(count + 1U, sizeof(*lpf));
    hit_flags = (unsigned char *)calloc(count + 1U, sizeof(*hit_flags));
    if (nums == NULL || sigma == NULL || lpf == NULL || hit_flags == NULL) {
        free(nums);
        free(sigma);
        free(lpf);
        free(hit_flags);
        return false;
    }

    count = 0U;
    if (lo == 1U) {
        nums[count++] = 1U;
    }
    if (has_even) {
        for (uint64_t n = first_even; n <= hi; n += 2U) {
            nums[count++] = n;
            if (UINT64_MAX - n < 2U) {
                break;
            }
        }
    }
    nums[count] = 0U;
    tel->candidates_tested = (uint64_t)count;
    tel->odd_skipped = tel->raw_width - tel->candidates_tested;

    *nums_out = nums;
    *sigma_out = sigma;
    *lpf_out = lpf;
    *hit_flags_out = hit_flags;
    *count_out = count;
    return true;
}

typedef struct {
    uint64_t rejected_qminus1;
    uint64_t mr_count;
    uint64_t rho_count;
    uint64_t trial_only_count;
    uint64_t mr_path_count;
    uint64_t rho_path_count;
    uint64_t unresolved_count;
    uint64_t failure_count;
    uint64_t hit_count;
    uint64_t lpf_sigma_max;
    uint64_t lpf_sigma_le_1e3;
    uint64_t lpf_sigma_le_1e6;
    uint64_t lpf_sigma_le_1e9;
    uint64_t lpf_sigma_gt_1e9;
} sigma_worker_stats_t;

static void sigma_stats_merge(sigma_worker_stats_t *dst, const sigma_worker_stats_t *src)
{
    dst->rejected_qminus1 += src->rejected_qminus1;
    dst->mr_count += src->mr_count;
    dst->rho_count += src->rho_count;
    dst->trial_only_count += src->trial_only_count;
    dst->mr_path_count += src->mr_path_count;
    dst->rho_path_count += src->rho_path_count;
    dst->unresolved_count += src->unresolved_count;
    dst->failure_count += src->failure_count;
    dst->hit_count += src->hit_count;
    if (src->lpf_sigma_max > dst->lpf_sigma_max) {
        dst->lpf_sigma_max = src->lpf_sigma_max;
    }
    dst->lpf_sigma_le_1e3 += src->lpf_sigma_le_1e3;
    dst->lpf_sigma_le_1e6 += src->lpf_sigma_le_1e6;
    dst->lpf_sigma_le_1e9 += src->lpf_sigma_le_1e9;
    dst->lpf_sigma_gt_1e9 += src->lpf_sigma_gt_1e9;
}

static void process_sigma_one(const uint64_t *orig_n, const uint64_t *sigma,
                              unsigned char *hit_flags, size_t i,
                              sigma_worker_stats_t *stats)
{
    factor_list_t sf;
    factor_ctx_t ctx = {0};
    uint64_t phi = 0;
    uint64_t lpf_sigma = 1U;
    bool rejected = false;

    if (!factor_u64(sigma[i], &sf, &ctx)) {
        stats->unresolved_count++;
        return;
    }
    factor_sort(&sf);
    lpf_sigma = factor_largest_prime(&sf);
    if (lpf_sigma > stats->lpf_sigma_max) {
        stats->lpf_sigma_max = lpf_sigma;
    }
    if (lpf_sigma <= UINT64_C(1000)) {
        stats->lpf_sigma_le_1e3++;
    } else if (lpf_sigma <= UINT64_C(1000000)) {
        stats->lpf_sigma_le_1e6++;
    } else if (lpf_sigma <= UINT64_C(1000000000)) {
        stats->lpf_sigma_le_1e9++;
    } else {
        stats->lpf_sigma_gt_1e9++;
    }

    stats->mr_count += ctx.mr_count;
    stats->rho_count += ctx.rho_count;
    if (ctx.used_rho) {
        stats->rho_path_count++;
    } else if (ctx.used_mr) {
        stats->mr_path_count++;
    } else {
        stats->trial_only_count++;
    }

    for (size_t j = 0; j < sf.len; j++) {
        uint64_t qminus1 = sf.f[j].p - 1U;
        /*
         * Necessary condition: if phi(sigma(n)) == n and q divides
         * sigma(n), then (q - 1) divides phi(sigma(n)) and therefore n.
         * This is an exact reject, not a heuristic filter.
         */
        if (qminus1 != 0U && (orig_n[i] % qminus1) != 0U) {
            rejected = true;
            break;
        }
    }
    if (rejected) {
        stats->rejected_qminus1++;
        return;
    }
    if (!phi_from_factorization(sigma[i], &sf, &phi)) {
        stats->failure_count++;
        return;
    }
    if (phi == orig_n[i]) {
        hit_flags[i] = 1U;
        stats->hit_count++;
    }
}

static bool process_sigma_candidates(const uint64_t *orig_n, const uint64_t *sigma,
                                     uint64_t *lpf_n, unsigned char *hit_flags,
                                     size_t count, int threads, telemetry_t *tel)
{
    sigma_worker_stats_t totals = {0};

    (void)lpf_n;

#ifdef _OPENMP
    if (threads > 0) {
        omp_set_num_threads(threads);
    }
#else
    (void)threads;
#endif

#ifdef _OPENMP
#pragma omp parallel if(count > 512U)
    {
        sigma_worker_stats_t local = {0};

#pragma omp for schedule(dynamic, 64)
        for (size_t i = 0; i < count; i++) {
            process_sigma_one(orig_n, sigma, hit_flags, i, &local);
        }

#pragma omp critical(a163665_sigma_stats_merge)
        {
            sigma_stats_merge(&totals, &local);
        }
    }
#else
    {
        sigma_worker_stats_t local = {0};

        for (size_t i = 0; i < count; i++) {
            process_sigma_one(orig_n, sigma, hit_flags, i, &local);
        }
        sigma_stats_merge(&totals, &local);
    }
#endif

    tel->rejected_qminus1 += totals.rejected_qminus1;
    tel->mr_count += totals.mr_count;
    tel->rho_count += totals.rho_count;
    tel->trial_only_count += totals.trial_only_count;
    tel->mr_path_count += totals.mr_path_count;
    tel->rho_path_count += totals.rho_path_count;
    tel->unresolved_count += totals.unresolved_count;
    tel->failure_count += totals.failure_count;
    tel->hit_count += totals.hit_count;
    if (totals.lpf_sigma_max > tel->lpf_sigma_max) {
        tel->lpf_sigma_max = totals.lpf_sigma_max;
    }
    tel->lpf_sigma_le_1e3 += totals.lpf_sigma_le_1e3;
    tel->lpf_sigma_le_1e6 += totals.lpf_sigma_le_1e6;
    tel->lpf_sigma_le_1e9 += totals.lpf_sigma_le_1e9;
    tel->lpf_sigma_gt_1e9 += totals.lpf_sigma_gt_1e9;
    return true;
}

static bool collect_hits(const uint64_t *orig_n, const uint64_t *sigma,
                         const unsigned char *hit_flags, size_t count,
                         hit_t **hits_out, size_t *hit_count_out)
{
    size_t hit_count = 0;
    size_t pos = 0;
    hit_t *hits = NULL;

    for (size_t i = 0; i < count; i++) {
        if (hit_flags[i] != 0U) {
            hit_count++;
        }
    }
    hits = (hit_t *)calloc(hit_count + 1U, sizeof(*hits));
    if (hits == NULL) {
        return false;
    }
    for (size_t i = 0; i < count; i++) {
        factor_list_t nf;
        factor_list_t sf;
        factor_ctx_t ctx = {0};

        if (hit_flags[i] == 0U) {
            continue;
        }
        hits[pos].n = orig_n[i];
        hits[pos].sigma_n = sigma[i];
        hits[pos].phi_sigma_n = orig_n[i];
        if (!nth_prime(orig_n[i], &hits[pos].k)) {
            free(hits);
            return false;
        }
        if (!factor_u64(orig_n[i], &nf, &ctx) || !factor_u64(sigma[i], &sf, &ctx)) {
            free(hits);
            return false;
        }
        factor_sort(&nf);
        factor_sort(&sf);
        if (!factor_to_string(&nf, hits[pos].factor_n, sizeof(hits[pos].factor_n)) ||
            !factor_to_string(&sf, hits[pos].factor_sigma, sizeof(hits[pos].factor_sigma))) {
            free(hits);
            return false;
        }
        pos++;
    }
    *hits_out = hits;
    *hit_count_out = hit_count;
    return true;
}

static bool write_hits_file(const char *prefix, const hit_t *hits, size_t hit_count,
                            char *hits_path, size_t path_cap)
{
    char tmp_path[4096];
    FILE *out = NULL;

    if (snprintf(hits_path, path_cap, "%s.hits.tsv", prefix) >= (int)path_cap ||
        snprintf(tmp_path, sizeof(tmp_path), "%s.hits.tsv.tmp", prefix) >= (int)sizeof(tmp_path)) {
        return false;
    }
    out = open_exclusive_temp(tmp_path);
    if (out == NULL) {
        return false;
    }
    fprintf(out, "n\tk\tsigma_n\tphi_sigma_n\tfactor_n\tfactor_sigma_n\n");
    for (size_t i = 0; i < hit_count; i++) {
        fprintf(out, "%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%s\t%s\n",
                hits[i].n, hits[i].k, hits[i].sigma_n, hits[i].phi_sigma_n,
                hits[i].factor_n, hits[i].factor_sigma);
    }
    bool write_ok = true;
    if (fflush(out) != 0) {
        write_ok = false;
    }
    if (fsync(fileno(out)) != 0) {
        write_ok = false;
    }
    if (fclose(out) != 0) {
        write_ok = false;
    }
    if (!write_ok) {
        unlink(tmp_path);
        return false;
    }
    return publish_file_noreplace(tmp_path, hits_path);
}

static bool write_manifest(const char *prefix, const char *status, uint64_t lo,
                           uint64_t hi, int threads, const telemetry_t *tel,
                           double wall_seconds, const char *argv0,
                           const char *hits_path, time_t started_unix)
{
    char manifest_path[4096];
    char tmp_path[4096];
    bool binary_ok = false;
    bool source_ok = false;
    bool hits_ok = false;
    uint64_t binary_sum = fnv1a_file64(argv0, &binary_ok);
    uint64_t source_sum = 0;
    uint64_t hits_sum = fnv1a_file64(hits_path, &hits_ok);
    char binary_sha256[65];
    char source_sha256[65];
    char hits_sha256[65];
    FILE *out = NULL;
    time_t ended = time(NULL);
    double cps = (wall_seconds > 0.0) ? ((double)tel->candidates_tested / wall_seconds) : 0.0;

    source_ok = resolve_source_checksums(&source_sum, source_sha256);
    if (!source_ok) {
        fprintf(stderr, "cannot checksum source path %s\n", A163665_SOURCE_PATH);
        return false;
    }
    if (A163665_BUILD_SOURCE_SHA256[0] != '\0' &&
        strcmp(source_sha256, A163665_BUILD_SOURCE_SHA256) != 0) {
        fprintf(stderr, "runtime source checksum does not match build-time source checksum\n");
        return false;
    }
    if (!binary_ok || !hits_ok ||
        !sha256_file_hex(argv0, binary_sha256) ||
        !sha256_file_hex(hits_path, hits_sha256)) {
        fprintf(stderr, "cannot compute artifact checksums\n");
        return false;
    }
    if (snprintf(manifest_path, sizeof(manifest_path), "%s.manifest.tsv", prefix) >=
            (int)sizeof(manifest_path) ||
        snprintf(tmp_path, sizeof(tmp_path), "%s.manifest.tsv.tmp", prefix) >=
            (int)sizeof(tmp_path)) {
        return false;
    }
    out = open_exclusive_temp(tmp_path);
    if (out == NULL) {
        return false;
    }
    fputs(manifest_header, out);
    fprintf(out,
            "%s\t%s\t%s\t%" PRIu64 "\t%" PRIu64 "\tn_eq_1_or_even_gt_1\t%d\t%" PRIu64
            "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64
            "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64
            "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64
            "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64
            "\t%" PRIu64 "\t%" PRIu64 "\t%.9f\t%.6f\t%016" PRIx64 "\t%016" PRIx64
            "\t%016" PRIx64 "\t%s\t%s\t%s\t%" PRIu64 "\t%" PRIu64 "\n",
            A163665_SEQ_ID, A163665_VERSION, status, lo, hi, threads,
            tel->raw_width, tel->raw_width, tel->candidates_tested, tel->odd_skipped,
            tel->rejected_qminus1, tel->mr_count, tel->rho_count, tel->trial_only_count,
            tel->mr_path_count, tel->rho_path_count, tel->lpf_n_max, tel->lpf_sigma_max,
            tel->lpf_n_le_1e3, tel->lpf_n_le_1e6, tel->lpf_n_le_1e9, tel->lpf_n_gt_1e9,
            tel->lpf_sigma_le_1e3, tel->lpf_sigma_le_1e6, tel->lpf_sigma_le_1e9,
            tel->lpf_sigma_gt_1e9, tel->unresolved_count, tel->failure_count,
            tel->hit_count, wall_seconds, cps, binary_ok ? binary_sum : 0U,
            source_ok ? source_sum : 0U, hits_ok ? hits_sum : 0U,
            binary_sha256, source_sha256, hits_sha256, (uint64_t)started_unix,
            (uint64_t)ended);
    bool write_ok = true;
    if (fflush(out) != 0) {
        write_ok = false;
    }
    if (fsync(fileno(out)) != 0) {
        write_ok = false;
    }
    if (fclose(out) != 0) {
        write_ok = false;
    }
    if (!write_ok) {
        unlink(tmp_path);
        return false;
    }
    return publish_file_noreplace(tmp_path, manifest_path);
}

static bool run_segment(uint64_t lo, uint64_t hi, int threads, const char *prefix,
                        const char *argv0, bool quiet)
{
    telemetry_t tel = {0};
    uint64_t *work_n = NULL;
    uint64_t *orig_n = NULL;
    uint64_t *sigma = NULL;
    uint64_t *lpf_n = NULL;
    unsigned char *hit_flags = NULL;
    hit_t *hits = NULL;
    size_t count = 0;
    size_t hit_count = 0;
    char hits_path[4096];
    const char *status = "FAILED";
    bool ok = false;
    double t0 = monotonic_seconds();
    double t1 = 0.0;
    time_t started_unix = time(NULL);
    double wall_seconds = 0.0;

    if (threads < 1 || threads > 256) {
        fprintf(stderr, "threads must be in [1,256]\n");
        return false;
    }
    if (!output_paths_available(prefix)) {
        fprintf(stderr, "output prefix already exists or is not writable: %s\n", prefix);
        tel.failure_count++;
        goto done;
    }
    if (!prepare_candidates(lo, hi, &work_n, &sigma, &lpf_n, &hit_flags, &count, &tel)) {
        tel.failure_count++;
        goto done;
    }
    orig_n = (uint64_t *)calloc(count + 1U, sizeof(*orig_n));
    if (orig_n == NULL) {
        tel.failure_count++;
        goto done;
    }
    memcpy(orig_n, work_n, (count + 1U) * sizeof(*orig_n));
    if (count > 0U) {
        if (!segment_build_sigma(lo, hi, work_n, sigma, lpf_n, count, &tel)) {
            tel.failure_count++;
            goto done;
        }
        if (!process_sigma_candidates(orig_n, sigma, lpf_n, hit_flags, count, threads, &tel)) {
            tel.failure_count++;
            goto done;
        }
    }
    if (!collect_hits(orig_n, sigma, hit_flags, count, &hits, &hit_count)) {
        tel.failure_count++;
        goto done;
    }
    tel.hit_count = (uint64_t)hit_count;
    status = (tel.failure_count > 0U) ? "FAILED" :
             ((tel.unresolved_count > 0U) ? "UNRESOLVED" :
              ((tel.hit_count > 0U) ? "OK_HITS" : "OK_EMPTY"));
    if (!write_hits_file(prefix, hits, hit_count, hits_path, sizeof(hits_path))) {
        tel.failure_count++;
        status = "FAILED";
        goto done;
    }
    t1 = monotonic_seconds();
    wall_seconds = (t1 >= t0) ? (t1 - t0) : 0.0;
    if (!write_manifest(prefix, status, lo, hi, threads, &tel, wall_seconds, argv0,
                        hits_path, started_unix)) {
        tel.failure_count++;
        status = "FAILED";
        goto done;
    }
    ok = (strcmp(status, "OK_EMPTY") == 0 || strcmp(status, "OK_HITS") == 0);

done:
    if (!quiet) {
        fprintf(stderr, "%s [%" PRIu64 ",%" PRIu64 "] tested=%" PRIu64
                " hits=%" PRIu64 " unresolved=%" PRIu64 " failures=%" PRIu64 "\n",
                status, lo, hi, tel.candidates_tested, tel.hit_count,
                tel.unresolved_count, tel.failure_count);
    }
    free(work_n);
    free(orig_n);
    free(sigma);
    free(lpf_n);
    free(hit_flags);
    free(hits);
    return ok;
}

static bool read_bfile_terms(void)
{
    static const uint64_t new_k[8] = {
        UINT64_C(50685770167), UINT64_C(94206075193),
        UINT64_C(106882008719), UINT64_C(571264143487),
        UINT64_C(1272242879459), UINT64_C(2536961181761),
        UINT64_C(3715374607207), UINT64_C(9576704442589)
    };
    FILE *f = fopen("data/b163665.txt", "rb");
    char line[256];
    char expected[256];
    unsigned int idx = 0;
    bool saw_final_blank = false;

    if (f == NULL) {
        f = fopen("../data/b163665.txt", "rb");
        if (f == NULL) {
            return false;
        }
    }
    while (fgets(line, sizeof(line), f) != NULL) {
        if (line_has_cr(line) || strchr(line, '\n') == NULL) {
            fclose(f);
            return false;
        }
        if (saw_final_blank) {
            fclose(f);
            return false;
        }
        if (idx < 20U) {
            if (snprintf(expected, sizeof(expected), "%u %" PRIu64 "\n",
                         idx + 1U, known_k[idx]) >= (int)sizeof(expected) ||
                strcmp(line, expected) != 0) {
                fclose(f);
                return false;
            }
            idx++;
            continue;
        }
        if (idx < 28U) {
            const uint64_t value = new_k[idx - 20U];
            if (snprintf(expected, sizeof(expected), "%u %" PRIu64 "\n",
                         idx + 1U, value) >= (int)sizeof(expected) ||
                strcmp(line, expected) != 0) {
                fclose(f);
                return false;
            }
            idx++;
            continue;
        }
        if (strcmp(line, "\n") != 0) {
            fclose(f);
            return false;
        }
        saw_final_blank = true;
    }
    if (ferror(f) != 0) {
        fclose(f);
        return false;
    }
    fclose(f);
    return idx == 28U && saw_final_blank;
}

static bool selftest_known_terms(void)
{
    for (size_t i = 0; i < 20U; i++) {
        uint64_t sigma_n = 0;
        uint64_t phi_sigma_n = 0;
        uint64_t k = 0;

        if (!direct_phi_sigma(known_n[i], &sigma_n, &phi_sigma_n)) {
            fprintf(stderr, "selftest: failed predicate computation at index %zu\n", i + 1U);
            return false;
        }
        if (phi_sigma_n != known_n[i]) {
            fprintf(stderr, "selftest: known n mismatch at index %zu\n", i + 1U);
            return false;
        }
        if (!nth_prime(known_n[i], &k) || k != known_k[i]) {
            fprintf(stderr, "selftest: prime(n) mismatch at index %zu\n", i + 1U);
            return false;
        }
    }
    return true;
}

static bool selftest_low_range(void)
{
    size_t hit_pos = 0;
    for (uint64_t n = 1; n <= UINT64_C(100000); n++) {
        uint64_t sigma_n = 0;
        uint64_t phi_sigma_n = 0;
        bool is_hit = false;

        if (n > 1U && (n & 1U) != 0U) {
            continue;
        }
        if (!direct_phi_sigma(n, &sigma_n, &phi_sigma_n)) {
            return false;
        }
        is_hit = (phi_sigma_n == n);
        if (is_hit) {
            if (hit_pos >= 9U || n != known_n[hit_pos]) {
                fprintf(stderr, "selftest: unexpected low-range hit n=%" PRIu64 "\n", n);
                return false;
            }
            hit_pos++;
        }
    }
    if (hit_pos != 9U) {
        fprintf(stderr, "selftest: low-range hit count mismatch\n");
        return false;
    }
    return true;
}

static bool selftest_segment_vs_direct(const char *argv0)
{
    char prefix[256];
    char tmp_dir[] = "/tmp/a163665_selftest_XXXXXX";
    char hits_path[320] = {0};
    char manifest_path[320] = {0};
    FILE *f = NULL;
    char line[1024];
    uint64_t expected[32];
    size_t expected_count = 0;
    size_t seen_count = 0;
    bool ok = false;

    for (uint64_t n = 1; n <= UINT64_C(100000); n++) {
        uint64_t sigma_n = 0;
        uint64_t phi_sigma_n = 0;
        if (n > 1U && (n & 1U) != 0U) {
            continue;
        }
        if (!direct_phi_sigma(n, &sigma_n, &phi_sigma_n)) {
            return false;
        }
        if (phi_sigma_n == n) {
            if (expected_count >= (sizeof(expected) / sizeof(expected[0]))) {
                return false;
            }
            expected[expected_count++] = n;
        }
    }
    if (mkdtemp(tmp_dir) == NULL) {
        return false;
    }
    if (snprintf(prefix, sizeof(prefix), "%s/segment", tmp_dir) >= (int)sizeof(prefix)) {
        goto cleanup;
    }
    if (!run_segment(1U, UINT64_C(100000), 1, prefix, argv0, true)) {
        goto cleanup;
    }
    if (snprintf(hits_path, sizeof(hits_path), "%s.hits.tsv", prefix) >=
            (int)sizeof(hits_path) ||
        snprintf(manifest_path, sizeof(manifest_path), "%s.manifest.tsv", prefix) >=
            (int)sizeof(manifest_path)) {
        goto cleanup;
    }
    f = fopen(hits_path, "rb");
    if (f == NULL) {
        goto cleanup;
    }
    if (fgets(line, sizeof(line), f) == NULL) {
        goto cleanup;
    }
    while (fgets(line, sizeof(line), f) != NULL) {
        uint64_t n = 0;
        if (sscanf(line, "%" SCNu64 "\t", &n) != 1) {
            goto cleanup;
        }
        if (seen_count >= expected_count || n != expected[seen_count]) {
            goto cleanup;
        }
        seen_count++;
    }
    ok = (seen_count == expected_count);

cleanup:
    if (f != NULL) {
        fclose(f);
    }
    if (hits_path[0] != '\0') {
        unlink(hits_path);
    }
    if (manifest_path[0] != '\0') {
        unlink(manifest_path);
    }
    rmdir(tmp_dir);
    return ok;
}

static bool selftest_numeric_edges(void)
{
    uint64_t lim_sq = 0;
    uint64_t next_lim = 0;
    uint64_t next_lim_sq = 0;

    if (mod_add_u64(UINT64_MAX - 1U, 31U, UINT64_MAX) != 30U) {
        fprintf(stderr, "selftest: modular addition edge failed\n");
        return false;
    }
    if (!init_lehmer_table()) {
        return false;
    }
    if (!u64_mul_checked((uint64_t)g_lehmer.limit, (uint64_t)g_lehmer.limit, &lim_sq)) {
        return false;
    }
    next_lim = (uint64_t)g_lehmer.limit + 1U;
    if (!u64_mul_checked(next_lim, next_lim, &next_lim_sq)) {
        return false;
    }
    if (!lehmer_pi_domain_ok(lim_sq - 1U) || lehmer_pi_domain_ok(next_lim_sq)) {
        fprintf(stderr, "selftest: Lehmer domain guard edge failed\n");
        return false;
    }
    return true;
}

static int cmd_selftest(const char *argv0)
{
    int rc = 0;

    if (!read_bfile_terms()) {
        fprintf(stderr, "selftest: b-file replay failed\n");
        rc = 1;
        goto done;
    }
    if (!selftest_known_terms()) {
        rc = 1;
        goto done;
    }
    if (!selftest_low_range()) {
        rc = 1;
        goto done;
    }
    if (!selftest_segment_vs_direct(argv0)) {
        fprintf(stderr, "selftest: segmented/direct comparison failed\n");
        rc = 1;
        goto done;
    }
    if (!selftest_numeric_edges()) {
        rc = 1;
        goto done;
    }
    printf("selftest OK: 20 b-file terms, known indices, low-range cross-check, segmented=direct\n");
done:
    free_lehmer_table();
    return rc;
}

static int cmd_segment(int argc, char **argv)
{
    uint64_t lo = 0;
    uint64_t hi = 0;
    uint64_t threads64 = 0;
    const char *out = NULL;
    char reason[256];
    bool have_lo = false;
    bool have_hi = false;
    bool have_threads = false;
    bool have_out = false;

    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "--lo") == 0 && i + 1 < argc) {
            if (have_lo) {
                usage(stderr);
                return 2;
            }
            have_lo = true;
            if (!parse_u64(argv[++i], &lo)) {
                usage(stderr);
                return 2;
            }
        } else if (strcmp(argv[i], "--hi") == 0 && i + 1 < argc) {
            if (have_hi) {
                usage(stderr);
                return 2;
            }
            have_hi = true;
            if (!parse_u64(argv[++i], &hi)) {
                usage(stderr);
                return 2;
            }
        } else if (strcmp(argv[i], "--threads") == 0 && i + 1 < argc) {
            if (have_threads) {
                usage(stderr);
                return 2;
            }
            have_threads = true;
            if (!parse_u64(argv[++i], &threads64)) {
                usage(stderr);
                return 2;
            }
        } else if (strcmp(argv[i], "--out") == 0 && i + 1 < argc) {
            if (have_out) {
                usage(stderr);
                return 2;
            }
            out = argv[++i];
            have_out = true;
        } else {
            usage(stderr);
            return 2;
        }
    }
    if (!have_lo || !have_hi || !have_threads || out == NULL || lo == 0U || threads64 < 1U ||
        threads64 > 256U || hi < lo) {
        usage(stderr);
        return 2;
    }
    if (!segment_domain_supported(lo, hi, reason, sizeof(reason))) {
        fprintf(stderr, "%s\n", reason);
        return 2;
    }
    int rc = run_segment(lo, hi, (int)threads64, out, argv[0], false) ? 0 : 1;
    free_lehmer_table();
    return rc;
}

static bool parse_hex64(const char *s, uint64_t *out)
{
    uint64_t value = 0;

    if (s == NULL || strlen(s) != 16U) {
        return false;
    }
    for (const char *p = s; *p != '\0'; p++) {
        unsigned int digit = 0;
        if (*p >= '0' && *p <= '9') {
            digit = (unsigned int)(*p - '0');
        } else if (*p >= 'a' && *p <= 'f') {
            digit = (unsigned int)(*p - 'a') + 10U;
        } else {
            return false;
        }
        value = (value << 4U) | (uint64_t)digit;
    }
    *out = value;
    return true;
}

static bool replace_suffix(const char *path, const char *old_suffix, const char *new_suffix,
                           char *out, size_t cap)
{
    size_t path_len = strlen(path);
    size_t old_len = strlen(old_suffix);
    size_t new_len = strlen(new_suffix);
    size_t prefix_len = 0;

    if (path_len < old_len || strcmp(path + path_len - old_len, old_suffix) != 0) {
        return false;
    }
    prefix_len = path_len - old_len;
    if (prefix_len + new_len + 1U > cap) {
        return false;
    }
    memcpy(out, path, prefix_len);
    memcpy(out + prefix_len, new_suffix, new_len + 1U);
    return true;
}

static bool line_has_cr(const char *s)
{
    return strchr(s, '\r') != NULL;
}

static bool split_tsv_line_exact(char *line, size_t expected_fields, char **cols)
{
    size_t len = strlen(line);
    size_t col_count = 1U;

    if (expected_fields == 0U || len == 0U || line[len - 1U] != '\n' || line_has_cr(line)) {
        return false;
    }
    line[--len] = '\0';
    if (len == 0U || line[0] == '\t' || line[len - 1U] == '\t') {
        return false;
    }
    cols[0] = line;
    for (size_t i = 0; i < len; i++) {
        if (line[i] != '\t') {
            continue;
        }
        if (line[i + 1U] == '\t' || line[i + 1U] == '\0' || col_count >= expected_fields) {
            return false;
        }
        line[i] = '\0';
        cols[col_count++] = line + i + 1U;
    }
    return col_count == expected_fields;
}

static bool validate_hits_file(const char *path, uint64_t lo, uint64_t hi,
                               uint64_t expected_rows, uint64_t *rows)
{
    FILE *f = fopen(path, "rb");
    char line[4096];
    uint64_t count = 0;
    uint64_t prev_n = 0;
    bool have_prev = false;

    if (f == NULL) {
        return false;
    }
    if (fgets(line, sizeof(line), f) == NULL || line_has_cr(line) ||
        strcmp(line, "n\tk\tsigma_n\tphi_sigma_n\tfactor_n\tfactor_sigma_n\n") != 0) {
        fclose(f);
        return false;
    }
    while (fgets(line, sizeof(line), f) != NULL) {
        char *cols[6];
        uint64_t n = 0;
        uint64_t k = 0;
        uint64_t sigma_n = 0;
        uint64_t phi_sigma_n = 0;
        uint64_t calc_sigma = 0;
        uint64_t calc_phi = 0;
        uint64_t calc_k = 0;
        factor_list_t nf;
        factor_list_t sf;
        factor_ctx_t ctx = {0};
        char factor_n[A163665_STR_CAP];
        char factor_sigma[A163665_STR_CAP];

        if (!split_tsv_line_exact(line, 6U, cols)) {
            fclose(f);
            return false;
        }
        if (!parse_u64(cols[0], &n) || !parse_u64(cols[1], &k) ||
            !parse_u64(cols[2], &sigma_n) || !parse_u64(cols[3], &phi_sigma_n)) {
            fclose(f);
            return false;
        }
        if (n < lo || n > hi || (n > 1U && (n & 1U) != 0U) ||
            phi_sigma_n != n || (have_prev && n <= prev_n)) {
            fclose(f);
            return false;
        }
        if (!direct_phi_sigma(n, &calc_sigma, &calc_phi) ||
            calc_sigma != sigma_n || calc_phi != phi_sigma_n ||
            !nth_prime(n, &calc_k) || calc_k != k) {
            fclose(f);
            return false;
        }
        if (!factor_u64(n, &nf, &ctx) || !factor_u64(sigma_n, &sf, &ctx)) {
            fclose(f);
            return false;
        }
        factor_sort(&nf);
        factor_sort(&sf);
        if (!factor_to_string(&nf, factor_n, sizeof(factor_n)) ||
            !factor_to_string(&sf, factor_sigma, sizeof(factor_sigma)) ||
            strcmp(factor_n, cols[4]) != 0 || strcmp(factor_sigma, cols[5]) != 0) {
            fclose(f);
            return false;
        }
        count++;
        if (count > expected_rows) {
            fclose(f);
            return false;
        }
        prev_n = n;
        have_prev = true;
    }
    if (ferror(f) != 0) {
        fclose(f);
        return false;
    }
    fclose(f);
    if (count != expected_rows) {
        return false;
    }
    *rows = count;
    return true;
}

static int manifest_summary_cmp(const void *a, const void *b)
{
    const manifest_summary_t *ma = (const manifest_summary_t *)a;
    const manifest_summary_t *mb = (const manifest_summary_t *)b;

    if (ma->lo < mb->lo) {
        return -1;
    }
    if (ma->lo > mb->lo) {
        return 1;
    }
    if (ma->hi < mb->hi) {
        return -1;
    }
    if (ma->hi > mb->hi) {
        return 1;
    }
    return 0;
}

static bool parse_manifest_summary(const char *path, const char *argv0,
                                   manifest_summary_t *summary)
{
    FILE *f = fopen(path, "rb");
    char header[4096];
    char row[4096];
    char extra[4096];
    char *cols[40];
    uint64_t added = 0;
    char hits_path[4096];
    bool hits_ok = false;
    bool binary_ok = false;
    uint64_t hits_sum = 0;
    uint64_t binary_sum = 0;
    uint64_t source_sum = 0;
    uint64_t hit_rows = 0;
    uint64_t threads = 0;
    uint64_t segment_size = 0;
    uint64_t rejected_qminus1 = 0;
    uint64_t mr_count = 0;
    uint64_t rho_count = 0;
    uint64_t trial_only_count = 0;
    uint64_t mr_path_count = 0;
    uint64_t rho_path_count = 0;
    uint64_t lpf_n_max = 0;
    uint64_t lpf_sigma_max = 0;
    uint64_t lpf_n_le_1e3 = 0;
    uint64_t lpf_n_le_1e6 = 0;
    uint64_t lpf_n_le_1e9 = 0;
    uint64_t lpf_n_gt_1e9 = 0;
    uint64_t lpf_sigma_le_1e3 = 0;
    uint64_t lpf_sigma_le_1e6 = 0;
    uint64_t lpf_sigma_le_1e9 = 0;
    uint64_t lpf_sigma_gt_1e9 = 0;
    uint64_t started_unix = 0;
    uint64_t ended_unix = 0;
    uint64_t sigma_factored_count = 0;
    double wall_seconds = 0.0;
    double candidates_per_second = 0.0;
    double expected_cps = 0.0;
    double cps_tolerance = 0.0;
    char hits_sha256[65];
    char binary_sha256[65];
    char source_sha256[65];

    if (f == NULL) {
        fprintf(stderr, "cannot open %s: %s\n", path, strerror(errno));
        return false;
    }
    if (fgets(header, sizeof(header), f) == NULL || fgets(row, sizeof(row), f) == NULL ||
        line_has_cr(header) || line_has_cr(row)) {
        fclose(f);
        return false;
    }
    if (fgets(extra, sizeof(extra), f) != NULL) {
        fclose(f);
        return false;
    }
    fclose(f);
    if (strcmp(header, manifest_header) != 0) {
        return false;
    }
    if (!split_tsv_line_exact(row, 40U, cols)) {
        return false;
    }
    if (strcmp(cols[0], A163665_SEQ_ID) != 0 || strcmp(cols[1], A163665_VERSION) != 0 ||
        strcmp(cols[5], "n_eq_1_or_even_gt_1") != 0) {
        return false;
    }
    if (strcmp(cols[2], "OK_EMPTY") == 0) {
        summary->status_class = 0;
    } else if (strcmp(cols[2], "OK_HITS") == 0) {
        summary->status_class = 1;
    } else if (strcmp(cols[2], "UNRESOLVED") == 0) {
        summary->status_class = 2;
    } else if (strcmp(cols[2], "FAILED") == 0) {
        summary->status_class = 3;
    } else {
        return false;
    }
    if (!parse_u64(cols[3], &summary->lo) || !parse_u64(cols[4], &summary->hi) ||
        !parse_u64(cols[6], &threads) ||
        !parse_u64(cols[7], &segment_size) ||
        !parse_u64(cols[8], &summary->raw_width) ||
        !parse_u64(cols[9], &summary->candidates_tested) ||
        !parse_u64(cols[10], &summary->odd_skipped) ||
        !parse_u64(cols[11], &rejected_qminus1) ||
        !parse_u64(cols[12], &mr_count) ||
        !parse_u64(cols[13], &rho_count) ||
        !parse_u64(cols[14], &trial_only_count) ||
        !parse_u64(cols[15], &mr_path_count) ||
        !parse_u64(cols[16], &rho_path_count) ||
        !parse_u64(cols[17], &lpf_n_max) ||
        !parse_u64(cols[18], &lpf_sigma_max) ||
        !parse_u64(cols[19], &lpf_n_le_1e3) ||
        !parse_u64(cols[20], &lpf_n_le_1e6) ||
        !parse_u64(cols[21], &lpf_n_le_1e9) ||
        !parse_u64(cols[22], &lpf_n_gt_1e9) ||
        !parse_u64(cols[23], &lpf_sigma_le_1e3) ||
        !parse_u64(cols[24], &lpf_sigma_le_1e6) ||
        !parse_u64(cols[25], &lpf_sigma_le_1e9) ||
        !parse_u64(cols[26], &lpf_sigma_gt_1e9) ||
        !parse_u64(cols[27], &summary->unresolved_count) ||
        !parse_u64(cols[28], &summary->failure_count) ||
        !parse_u64(cols[29], &summary->hit_count) ||
        !parse_nonnegative_decimal(cols[30], &wall_seconds) ||
        !parse_nonnegative_decimal(cols[31], &candidates_per_second) ||
        !parse_hex64(cols[32], &summary->binary_fnv64) ||
        !parse_hex64(cols[33], &summary->source_fnv64) ||
        !parse_hex64(cols[34], &summary->hits_fnv64) ||
        !is_sha256_hex(cols[35]) || !is_sha256_hex(cols[36]) ||
        !is_sha256_hex(cols[37]) ||
        !parse_u64(cols[38], &started_unix) ||
        !parse_u64(cols[39], &ended_unix)) {
        return false;
    }
    memcpy(summary->binary_sha256, cols[35], sizeof(summary->binary_sha256));
    memcpy(summary->source_sha256, cols[36], sizeof(summary->source_sha256));
    memcpy(summary->hits_sha256, cols[37], sizeof(summary->hits_sha256));
    if (summary->lo == 0U || summary->hi < summary->lo ||
        threads < 1U || threads > 256U ||
        segment_size != summary->raw_width ||
        summary->raw_width != summary->hi - summary->lo + 1U ||
        summary->raw_width > A163665_MAX_SEGMENT_WIDTH ||
        summary->candidates_tested != expected_candidate_count(summary->lo, summary->hi) ||
        !u64_add_checked(summary->candidates_tested, summary->odd_skipped, &added) ||
        added != summary->raw_width ||
        rejected_qminus1 > summary->candidates_tested ||
        summary->hit_count > summary->candidates_tested ||
        summary->unresolved_count > summary->candidates_tested ||
        summary->failure_count > summary->candidates_tested ||
        summary->binary_fnv64 == 0U || summary->source_fnv64 == 0U ||
        summary->hits_fnv64 == 0U ||
        ended_unix < started_unix) {
        return false;
    }
    sigma_factored_count = summary->candidates_tested - summary->unresolved_count;
    expected_cps = (wall_seconds > 0.0) ?
                   ((double)summary->candidates_tested / wall_seconds) : 0.0;
    cps_tolerance = 0.000001 + (((expected_cps > 1.0) ? expected_cps : 1.0) * 1.0e-9);
    if (fabs(candidates_per_second - expected_cps) > cps_tolerance) {
        return false;
    }
    if (!u64_add_checked(lpf_n_le_1e3, lpf_n_le_1e6, &added) ||
        !u64_add_checked(added, lpf_n_le_1e9, &added) ||
        !u64_add_checked(added, lpf_n_gt_1e9, &added) ||
        added != summary->candidates_tested) {
        return false;
    }
    if (!u64_add_checked(lpf_sigma_le_1e3, lpf_sigma_le_1e6, &added) ||
        !u64_add_checked(added, lpf_sigma_le_1e9, &added) ||
        !u64_add_checked(added, lpf_sigma_gt_1e9, &added) ||
        !u64_add_checked(added, summary->unresolved_count, &added) ||
        added != summary->candidates_tested) {
        return false;
    }
    if (!lpf_bucket_max_consistent(lpf_n_max, lpf_n_le_1e3, lpf_n_le_1e6,
                                   lpf_n_le_1e9, lpf_n_gt_1e9,
                                   summary->candidates_tested) ||
        !lpf_bucket_max_consistent(lpf_sigma_max, lpf_sigma_le_1e3,
                                   lpf_sigma_le_1e6, lpf_sigma_le_1e9,
                                   lpf_sigma_gt_1e9, sigma_factored_count)) {
        return false;
    }
    if (!u64_add_checked(trial_only_count, mr_path_count, &added) ||
        !u64_add_checked(added, rho_path_count, &added) ||
        !u64_add_checked(added, summary->unresolved_count, &added) ||
        added != summary->candidates_tested) {
        return false;
    }
    if ((summary->candidates_tested == 0U && (lpf_n_max != 0U || lpf_sigma_max != 0U)) ||
        (summary->candidates_tested > 0U && lpf_n_max == 0U) ||
        (summary->candidates_tested > summary->unresolved_count && lpf_sigma_max == 0U) ||
        (summary->candidates_tested == summary->unresolved_count && lpf_sigma_max != 0U) ||
        (summary->candidates_tested == 0U && candidates_per_second != 0.0) ||
        (wall_seconds == 0.0 && candidates_per_second != 0.0) ||
        (mr_path_count > mr_count) ||
        (rho_path_count > rho_count)) {
        return false;
    }
    if ((summary->status_class == 0 || summary->status_class == 1) &&
        (summary->unresolved_count != 0U || summary->failure_count != 0U)) {
        return false;
    }
    if ((summary->status_class == 0 && summary->hit_count != 0U) ||
        (summary->status_class == 1 && summary->hit_count == 0U)) {
        return false;
    }
    if (summary->status_class == 2 && summary->unresolved_count == 0U) {
        return false;
    }
    if (summary->status_class == 3 && summary->failure_count == 0U) {
        return false;
    }
    if (!replace_suffix(path, ".manifest.tsv", ".hits.tsv", hits_path, sizeof(hits_path))) {
        return false;
    }
    binary_sum = fnv1a_file64(argv0, &binary_ok);
    if (!binary_ok || binary_sum != summary->binary_fnv64 ||
        !sha256_file_hex(argv0, binary_sha256) ||
        strcmp(binary_sha256, summary->binary_sha256) != 0 ||
        !resolve_source_checksums(&source_sum, source_sha256) ||
        source_sum != summary->source_fnv64 ||
        strcmp(source_sha256, summary->source_sha256) != 0) {
        return false;
    }
    hits_sum = fnv1a_file64(hits_path, &hits_ok);
    if (!hits_ok || hits_sum != summary->hits_fnv64 ||
        !sha256_file_hex(hits_path, hits_sha256) ||
        strcmp(hits_sha256, summary->hits_sha256) != 0 ||
        !validate_hits_file(hits_path, summary->lo, summary->hi,
                            summary->hit_count, &hit_rows)) {
        return false;
    }
    return true;
}

static int cmd_summarize(int argc, char **argv)
{
    uint64_t ok_segments = 0;
    uint64_t unresolved_segments = 0;
    uint64_t failed_segments = 0;
    uint64_t total_hits = 0;
    manifest_summary_t *summaries = NULL;
    size_t summary_count = 0;
    bool coverage_ok = true;
    int rc = 0;

    if (argc < 3) {
        usage(stderr);
        return 2;
    }
    summary_count = (size_t)(argc - 2);
    summaries = (manifest_summary_t *)calloc(summary_count, sizeof(*summaries));
    if (summaries == NULL) {
        return 1;
    }
    for (size_t i = 0; i < summary_count; i++) {
        if (!parse_manifest_summary(argv[i + 2], argv[0], &summaries[i])) {
            fprintf(stderr, "manifest validation failed: %s\n", argv[i + 2]);
            free(summaries);
            return 1;
        }
        if (summaries[i].status_class == 0 || summaries[i].status_class == 1) {
            ok_segments++;
        } else if (summaries[i].status_class == 2) {
            unresolved_segments++;
        } else {
            failed_segments++;
        }
        if (!u64_add_checked(total_hits, summaries[i].hit_count, &total_hits)) {
            free(summaries);
            return 1;
        }
    }
    qsort(summaries, summary_count, sizeof(*summaries), manifest_summary_cmp);
    for (size_t i = 1; i < summary_count; i++) {
        uint64_t expected_next = 0;
        if (summaries[i].lo <= summaries[i - 1U].hi ||
            !u64_add_checked(summaries[i - 1U].hi, 1U, &expected_next) ||
            summaries[i].lo != expected_next) {
            coverage_ok = false;
            break;
        }
    }
    printf("segments_ok=%" PRIu64 "\tsegments_unresolved=%" PRIu64
           "\tsegments_failed=%" PRIu64 "\ttotal_hits=%" PRIu64
           "\tcoverage=%s\tlo=%" PRIu64 "\thi=%" PRIu64 "\n",
           ok_segments, unresolved_segments, failed_segments, total_hits,
           coverage_ok ? "CONTIGUOUS" : "BROKEN", summaries[0].lo,
           summaries[summary_count - 1U].hi);
    rc = (failed_segments == 0U && unresolved_segments == 0U && coverage_ok) ? 0 : 1;
    free(summaries);
    free_lehmer_table();
    return rc;
}

int main(int argc, char **argv)
{
    (void)setlocale(LC_ALL, "C");

    if (argc < 2) {
        usage(stderr);
        return 2;
    }
    if (strcmp(argv[1], "selftest") == 0) {
        if (argc != 2) {
            usage(stderr);
            return 2;
        }
        return cmd_selftest(argv[0]);
    }
    if (strcmp(argv[1], "segment") == 0) {
        return cmd_segment(argc, argv);
    }
    if (strcmp(argv[1], "summarize") == 0) {
        return cmd_summarize(argc, argv);
    }
    usage(stderr);
    return 2;
}
