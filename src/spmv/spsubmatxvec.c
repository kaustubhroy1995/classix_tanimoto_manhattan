/* CLASSIX_T packed binary submatrix-vector and row-popcount operations.
 * Input: row-major np.packbits(..., axis=1, bitorder="little") uint8 arrays.
 * Padding bits must be zero. memcpy avoids alignment/aliasing assumptions.
 */
#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <intrin.h>

#define EXPORT __declspec(dllexport)

EXPORT int hardware_popcount_available(void) {
    int regs[4];
    __cpuid(regs, 1);
    return (regs[2] & (1 << 23)) != 0;
}

static __forceinline unsigned int portable_popcount64(uint64_t value) {
    value -= (value >> 1) & UINT64_C(0x5555555555555555);
    value = (value & UINT64_C(0x3333333333333333)) +
            ((value >> 2) & UINT64_C(0x3333333333333333));
    value = (value + (value >> 4)) & UINT64_C(0x0f0f0f0f0f0f0f0f);
    return (unsigned int)((value * UINT64_C(0x0101010101010101)) >> 56);
}

static __forceinline unsigned int popcount64(uint64_t value, int hardware) {
    return hardware ? (unsigned int)__popcnt64(value) : portable_popcount64(value);
}

EXPORT void row_popcounts(const uint8_t *packed, int nrows, int nbytes,
                          int32_t *output) {
    const int hardware = hardware_popcount_available();
    for (int row = 0; row < nrows; ++row) {
        const uint8_t *current = packed + (size_t)row * (size_t)nbytes;
        int32_t count = 0;
        int column = 0;
        for (; column + 8 <= nbytes; column += 8) {
            uint64_t word;
            memcpy(&word, current + column, sizeof(word));
            count += (int32_t)popcount64(word, hardware);
        }
        for (; column < nbytes; ++column) {
            count += (int32_t)popcount64(current[column], hardware);
        }
        output[row] = count;
    }
}

EXPORT void packed_submatxvec(const uint8_t *packed, int nbytes,
                             int start_row, int end_row,
                             const uint8_t *vector, int32_t *output) {
    const int hardware = hardware_popcount_available();
    for (int row = start_row; row < end_row; ++row) {
        const uint8_t *current = packed + (size_t)row * (size_t)nbytes;
        int32_t intersection = 0;
        int column = 0;
        for (; column + 8 <= nbytes; column += 8) {
            uint64_t lhs, rhs;
            memcpy(&lhs, current + column, sizeof(lhs));
            memcpy(&rhs, vector + column, sizeof(rhs));
            intersection += (int32_t)popcount64(lhs & rhs, hardware);
        }
        for (; column < nbytes; ++column) {
            intersection += (int32_t)popcount64(current[column] & vector[column], hardware);
        }
        output[row - start_row] = intersection;
    }
}
