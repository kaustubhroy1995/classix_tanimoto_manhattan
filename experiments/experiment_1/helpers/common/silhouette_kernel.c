/* Serial packed-binary Jaccard distances and full-reference silhouette reduction.
   Build with MSVC x64. No OpenMP, numerical libraries or clustering code. */
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#if defined(_MSC_VER)
#include <intrin.h>
#define API __declspec(dllexport)
#define INLINE __forceinline
#else
#define API __attribute__((visibility("default")))
#define INLINE inline
#endif

static INLINE uint64_t software_popcount(uint64_t x) {
    x -= (x >> 1) & UINT64_C(0x5555555555555555);
    x = (x & UINT64_C(0x3333333333333333)) + ((x >> 2) & UINT64_C(0x3333333333333333));
    x = (x + (x >> 4)) & UINT64_C(0x0f0f0f0f0f0f0f0f);
    return (x * UINT64_C(0x0101010101010101)) >> 56;
}

API int silhouette_has_hardware_popcount(void) {
#if defined(_MSC_VER) && defined(_M_X64)
    int info[4];
    __cpuid(info, 1);
    return (info[2] & (1 << 23)) != 0;
#elif defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))
    return __builtin_cpu_supports("popcnt") != 0;
#else
    return 0;
#endif
}

static INLINE uint64_t popcount(uint64_t x, int hardware) {
#if defined(_MSC_VER) && defined(_M_X64)
    if (hardware) return (uint64_t)__popcnt64(x);
#else
    (void)hardware;
#endif
    return software_popcount(x);
}

static INLINE uint64_t load_word(const uint8_t *p) {
    uint64_t result;
    memcpy(&result, p, sizeof(result));
    return result;
}

API int packed_row_popcounts(const uint8_t *packed, int64_t n, int64_t bytes,
                             int64_t *counts) {
    int64_t row, k;
    int hardware = silhouette_has_hardware_popcount();
    if (!packed || !counts || n < 1 || bytes < 1) return 1;
    if ((uint64_t)n > SIZE_MAX / (uint64_t)bytes) return 2;
    for (row = 0; row < n; ++row) {
        const uint8_t *p = packed + (size_t)row * (size_t)bytes;
        uint64_t total = 0;
        for (k = 0; k + 8 <= bytes; k += 8) total += popcount(load_word(p + k), hardware);
        for (; k < bytes; ++k) total += software_popcount(p[k]);
        counts[row] = (int64_t)total;
    }
    return 0;
}

API int tanimoto_distance_block(const uint8_t *packed, const int64_t *popcounts,
                                int64_t n, int64_t bytes, const int64_t *queries,
                                int64_t m, double *output) {
    int64_t qi, row, k;
    int hardware = silhouette_has_hardware_popcount();
    if (!packed || !popcounts || !queries || !output || n < 1 || bytes < 1 || m < 0 || m > 128) return 1;
    if ((uint64_t)n > SIZE_MAX / (uint64_t)bytes || (uint64_t)n > SIZE_MAX / sizeof(double) / (uint64_t)(m ? m : 1)) return 2;
    for (qi = 0; qi < m; ++qi) {
        int64_t q = queries[qi];
        const uint8_t *a;
        double *out = output + (size_t)qi * (size_t)n;
        if (q < 0 || q >= n) return 3;
        a = packed + (size_t)q * (size_t)bytes;
        for (row = 0; row < n; ++row) {
            const uint8_t *b = packed + (size_t)row * (size_t)bytes;
            uint64_t intersection = 0;
            int64_t union_count;
            for (k = 0; k + 32 <= bytes; k += 32) {
                intersection += popcount(load_word(a+k) & load_word(b+k), hardware);
                intersection += popcount(load_word(a+k+8) & load_word(b+k+8), hardware);
                intersection += popcount(load_word(a+k+16) & load_word(b+k+16), hardware);
                intersection += popcount(load_word(a+k+24) & load_word(b+k+24), hardware);
            }
            for (; k + 8 <= bytes; k += 8) intersection += popcount(load_word(a+k) & load_word(b+k), hardware);
            for (; k < bytes; ++k) intersection += software_popcount((uint64_t)(a[k] & b[k]));
            union_count = popcounts[q] + popcounts[row] - (int64_t)intersection;
            out[row] = union_count ? (double)(union_count - (int64_t)intersection) / (double)union_count : 0.0;
        }
    }
    return 0;
}

API int reduce_silhouette_block(const double *distances, int64_t m, int64_t n,
                                const int32_t *labels, const int64_t *counts,
                                int64_t clusters, const int64_t *queries,
                                double *silhouettes, double *within, double *nearest) {
    int64_t qi, row, cluster;
    double *sums;
    if (!distances || !labels || !counts || !queries || !silhouettes || !within || !nearest ||
        m < 0 || m > 128 || n < 1 || clusters < 2) return 1;
    if ((uint64_t)clusters > SIZE_MAX / sizeof(double)) return 2;
    for (cluster = 0; cluster < clusters; ++cluster) if (counts[cluster] < 1) return 3;
    sums = (double *)malloc((size_t)clusters * sizeof(double));
    if (!sums) return 4;
    for (qi = 0; qi < m; ++qi) {
        int64_t q = queries[qi];
        int32_t own;
        const double *distance = distances + (size_t)qi * (size_t)n;
        double a, b = DBL_MAX;
        if (q < 0 || q >= n) { free(sums); return 5; }
        own = labels[q];
        if (own < 0 || own >= clusters) { free(sums); return 6; }
        memset(sums, 0, (size_t)clusters * sizeof(double));
        for (row = 0; row < n; ++row) {
            int32_t label = labels[row];
            if (label < -1 || label >= clusters) { free(sums); return 7; }
            if (label >= 0 && row != q) sums[label] += distance[row];
        }
        a = counts[own] > 1 ? sums[own] / (double)(counts[own] - 1) : 0.0;
        for (cluster = 0; cluster < clusters; ++cluster) {
            double average;
            if (cluster == own) continue;
            average = sums[cluster] / (double)counts[cluster];
            if (average < b) b = average;
        }
        within[qi] = a;
        nearest[qi] = b;
        silhouettes[qi] = counts[own] == 1 || (a == 0.0 && b == 0.0) ? 0.0 : (b - a) / (a > b ? a : b);
    }
    free(sums);
    return 0;
}
