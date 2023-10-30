#ifndef SSE2RVV_H
#define SSE2RVV_H

// This header file provides a simple API translation layer
// between SSE intrinsics to their corresponding RVV versions

/*
 * sse2rvv is freely redistributable under the MIT License.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/* Tunable configurations */

/* Enable precise implementation of math operations
 * This would slow down the computation a bit, but gives consistent result with
 * x86 SSE. (e.g. would solve a hole or NaN pixel in the rendering result)
 */
/* _mm_min|max_ps|ss|pd|sd */
#ifndef SSE2RVV_PRECISE_MINMAX
#define SSE2RVV_PRECISE_MINMAX (0)
#endif
/* _mm_rcp_ps and _mm_div_ps */
#ifndef SSE2RVV_PRECISE_DIV
#define SSE2RVV_PRECISE_DIV (0)
#endif
/* _mm_sqrt_ps and _mm_rsqrt_ps */
#ifndef SSE2RVV_PRECISE_SQRT
#define SSE2RVV_PRECISE_SQRT (0)
#endif
/* _mm_dp_pd */
#ifndef SSE2RVV_PRECISE_DP
#define SSE2RVV_PRECISE_DP (0)
#endif

/* Enable inclusion of windows.h on MSVC platforms
 * This makes _mm_clflush functional on windows, as there is no builtin.
 */
#ifndef SSE2RVV_INCLUDE_WINDOWS_H
#define SSE2RVV_INCLUDE_WINDOWS_H (0)
#endif

/* compiler specific definitions */
#if defined(__GNUC__) || defined(__clang__)
#pragma push_macro("FORCE_INLINE")
#pragma push_macro("ALIGN_STRUCT")
#define FORCE_INLINE static inline __attribute__((always_inline))
#define ALIGN_STRUCT(x) __attribute__((aligned(x)))
#define _sse2rvv_likely(x) __builtin_expect(!!(x), 1)
#define _sse2rvv_unlikely(x) __builtin_expect(!!(x), 0)
#elif defined(_MSC_VER)
#if _MSVC_TRADITIONAL
#error Using the traditional MSVC preprocessor is not supported! Use /Zc:preprocessor instead.
#endif
#ifndef FORCE_INLINE
#define FORCE_INLINE static inline
#endif
#ifndef ALIGN_STRUCT
#define ALIGN_STRUCT(x) __declspec(align(x))
#endif
#define _sse2rvv_likely(x) (x)
#define _sse2rvv_unlikely(x) (x)
#else
#pragma message("Macro name collisions may happen with unsupported compilers.")
#endif

/* C language does not allow initializing a variable with a function call. */
#ifdef __cplusplus
#define _sse2rvv_const static const
#else
#define _sse2rvv_const const
#endif

#include <riscv_vector.h>
#include <stdint.h>
#include <stdlib.h>

#if defined(__GNUC__) || defined(__clang__)
#define _sse2rvv_define0(type, s, body) \
    __extension__({                     \
        type _a = (s);                  \
        body                            \
    })
#define _sse2rvv_define1(type, s, body) \
    __extension__({                     \
        type _a = (s);                  \
        body                            \
    })
#define _sse2rvv_define2(type, a, b, body) \
    __extension__({                        \
        type _a = (a), _b = (b);           \
        body                               \
    })
#define _sse2rvv_return(ret) (ret)
#else
#define _sse2rvv_define0(type, a, body) [=](type _a) { body }(a)
#define _sse2rvv_define1(type, a, body) [](type _a) { body }(a)
#define _sse2rvv_define2(type, a, b, body) \
    [](type _a, type _b) { body }((a), (b))
#define _sse2rvv_return(ret) return ret
#endif

#define _sse2rvv_init(...) \
    {                      \
        __VA_ARGS__        \
    }

/* Compiler barrier */
#if defined(_MSC_VER)
#define SSE2RVV_BARRIER() _ReadWriteBarrier()
#else
#define SSE2RVV_BARRIER()                      \
    do {                                       \
        __asm__ __volatile__("" ::: "memory"); \
        (void) 0;                              \
    } while (0)
#endif

/* Memory barriers
 * __atomic_thread_fence does not include a compiler barrier; instead,
 * the barrier is part of __atomic_load/__atomic_store's "volatile-like"
 * semantics.
 */
#if defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 201112L)
#include <stdatomic.h>
#endif

FORCE_INLINE void _sse2rvv_smp_mb(void)
{
    SSE2RVV_BARRIER();
#if defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 201112L) && \
    !defined(__STDC_NO_ATOMICS__)
    atomic_thread_fence(memory_order_seq_cst);
#elif defined(__GNUC__) || defined(__clang__)
    __atomic_thread_fence(__ATOMIC_SEQ_CST);
#endif
}

/* "__has_builtin" can be used to query support for built-in functions
 * provided by gcc/clang and other compilers that support it.
 */
#ifndef __has_builtin /* GCC prior to 10 or non-clang compilers */
/* Compatibility with gcc <= 9 */
#if defined(__GNUC__) && (__GNUC__ <= 9)
#define __has_builtin(x) HAS##x
#define HAS__builtin_popcount 1
#define HAS__builtin_popcountll 1

// __builtin_shuffle introduced in GCC 4.7.0
#if (__GNUC__ >= 5) || ((__GNUC__ == 4) && (__GNUC_MINOR__ >= 7))
#define HAS__builtin_shuffle 1
#else
#define HAS__builtin_shuffle 0
#endif

#define HAS__builtin_shufflevector 0
#define HAS__builtin_nontemporal_store 0
#else
#define __has_builtin(x) 0
#endif
#endif

/**
 * MACRO for shuffle parameter for _mm_shuffle_ps().
 * Argument fp3 is a digit[0123] that represents the fp from argument "b"
 * of mm_shuffle_ps that will be placed in fp3 of result. fp2 is the same
 * for fp2 in result. fp1 is a digit[0123] that represents the fp from
 * argument "a" of mm_shuffle_ps that will be places in fp1 of result.
 * fp0 is the same for fp0 of result.
 */
#define _MM_SHUFFLE(fp3, fp2, fp1, fp0) \
    (((fp3) << 6) | ((fp2) << 4) | ((fp1) << 2) | ((fp0)))

#if __has_builtin(__builtin_shufflevector)
#define _sse2rvv_shuffle(type, a, b, ...) \
    __builtin_shufflevector(a, b, __VA_ARGS__)
#elif __has_builtin(__builtin_shuffle)
#define _sse2rvv_shuffle(type, a, b, ...) \
    __extension__({                       \
        type tmp = {__VA_ARGS__};         \
        __builtin_shuffle(a, b, tmp);     \
    })
#endif

#ifdef _sse2rvv_shuffle
#define vshuffle_s16(a, b, ...) _sse2rvv_shuffle(int16x4_t, a, b, __VA_ARGS__)
#define vshuffleq_s16(a, b, ...) _sse2rvv_shuffle(int16x8_t, a, b, __VA_ARGS__)
#define vshuffle_s32(a, b, ...) _sse2rvv_shuffle(int32x2_t, a, b, __VA_ARGS__)
#define vshuffleq_s32(a, b, ...) _sse2rvv_shuffle(int32x4_t, a, b, __VA_ARGS__)
#define vshuffle_s64(a, b, ...) _sse2rvv_shuffle(int64x1_t, a, b, __VA_ARGS__)
#define vshuffleq_s64(a, b, ...) _sse2rvv_shuffle(int64x2_t, a, b, __VA_ARGS__)
#endif

/* Rounding mode macros. */
#define _MM_FROUND_TO_NEAREST_INT 0x00
#define _MM_FROUND_TO_NEG_INF 0x01
#define _MM_FROUND_TO_POS_INF 0x02
#define _MM_FROUND_TO_ZERO 0x03
#define _MM_FROUND_CUR_DIRECTION 0x04
#define _MM_FROUND_NO_EXC 0x08
#define _MM_FROUND_RAISE_EXC 0x00
#define _MM_FROUND_NINT (_MM_FROUND_TO_NEAREST_INT | _MM_FROUND_RAISE_EXC)
#define _MM_FROUND_FLOOR (_MM_FROUND_TO_NEG_INF | _MM_FROUND_RAISE_EXC)
#define _MM_FROUND_CEIL (_MM_FROUND_TO_POS_INF | _MM_FROUND_RAISE_EXC)
#define _MM_FROUND_TRUNC (_MM_FROUND_TO_ZERO | _MM_FROUND_RAISE_EXC)
#define _MM_FROUND_RINT (_MM_FROUND_CUR_DIRECTION | _MM_FROUND_RAISE_EXC)
#define _MM_FROUND_NEARBYINT (_MM_FROUND_CUR_DIRECTION | _MM_FROUND_NO_EXC)
#define _MM_ROUND_NEAREST 0x0000
#define _MM_ROUND_DOWN 0x2000
#define _MM_ROUND_UP 0x4000
#define _MM_ROUND_TOWARD_ZERO 0x6000
/* Flush zero mode macros. */
#define _MM_FLUSH_ZERO_MASK 0x8000
#define _MM_FLUSH_ZERO_ON 0x8000
#define _MM_FLUSH_ZERO_OFF 0x0000
/* Denormals are zeros mode macros. */
#define _MM_DENORMALS_ZERO_MASK 0x0040
#define _MM_DENORMALS_ZERO_ON 0x0040
#define _MM_DENORMALS_ZERO_OFF 0x0000

/* indicate immediate constant argument in a given range */
#define __constrange(a, b) const

/* A few intrinsics accept traditional data types like ints or floats, but
 * most operate on data types that are specific to SSE.
 * If a vector type ends in d, it contains doubles, and if it does not have
 * a suffix, it contains floats. An integer vector type can contain any type
 * of integer, from chars to shorts to unsigned long longs.
 */
typedef vint64m1_t __m64;
typedef vfloat32m1_t __m128;  /* 128-bit vector containing 4 floats */
typedef vfloat64m1_t __m128d; /* 128-bit vector containing 2 doubles */
typedef vint64m1_t __m128i;   /* 128-bit vector containing integers */

// __int64 is defined in the Intrinsics Guide which maps to different datatype
// in different data model
#if !(defined(_WIN32) || defined(_WIN64) || defined(__int64))
#if (defined(__x86_64__) || defined(__i386__))
#define __int64 long long
#else
#define __int64 int64_t
#endif
#endif

// A struct is defined in this header file called 'SIMDVec' which can be used
// by applications which attempt to access the contents of an __m128 struct
// directly.  It is important to note that accessing the __m128 struct directly
// is bad coding practice by Microsoft: @see:
// https://learn.microsoft.com/en-us/cpp/cpp/m128
//
// However, some legacy source code may try to access the contents of an __m128
// struct directly so the developer can use the SIMDVec as an alias for it.  Any
// casting must be done manually by the developer, as you cannot cast or
// otherwise alias the base NEON data type for intrinsic operations.
//
// union intended to allow direct access to an __m128 variable using the names
// that the MSVC compiler provides.  This union should really only be used when
// trying to access the members of the vector as integer values.  GCC/clang
// allow native access to the float members through a simple array access
// operator (in C since 4.6, in C++ since 4.8).
//
// Ideally direct accesses to SIMD vectors should not be used since it can cause
// a performance hit.  If it really is needed however, the original __m128
// variable can be aliased with a pointer to this union and used to access
// individual components.  The use of this union should be hidden behind a macro
// that is used throughout the codebase to access the members instead of always
// declaring this type of variable.
typedef union ALIGN_STRUCT(16) SIMDVec {
    float m128_f32[4];     // as floats - DON'T USE. Added for convenience.
    int8_t m128_i8[16];    // as signed 8-bit integers.
    int16_t m128_i16[8];   // as signed 16-bit integers.
    int32_t m128_i32[4];   // as signed 32-bit integers.
    int64_t m128_i64[2];   // as signed 64-bit integers.
    uint8_t m128_u8[16];   // as unsigned 8-bit integers.
    uint16_t m128_u16[8];  // as unsigned 16-bit integers.
    uint32_t m128_u32[4];  // as unsigned 32-bit integers.
    uint64_t m128_u64[2];  // as unsigned 64-bit integers.
} SIMDVec;

/* SSE macros */
// #define _MM_GET_FLUSH_ZERO_MODE _sse2rvv_mm_get_flush_zero_mode
// #define _MM_SET_FLUSH_ZERO_MODE _sse2rvv_mm_set_flush_zero_mode
// #define _MM_GET_DENORMALS_ZERO_MODE _sse2rvv_mm_get_denormals_zero_mode
// #define _MM_SET_DENORMALS_ZERO_MODE _sse2rvv_mm_set_denormals_zero_mode

// Function declaration
// SSE
FORCE_INLINE unsigned int _MM_GET_ROUNDING_MODE(void);
// FORCE_INLINE __m128 _mm_move_ss(__m128, __m128);
// FORCE_INLINE __m128 _mm_or_ps(__m128, __m128);
// FORCE_INLINE __m128 _mm_set_ps1(float);
// FORCE_INLINE __m128 _mm_setzero_ps(void);
// SSE2
// FORCE_INLINE __m128i _mm_and_si128(__m128i, __m128i);
// FORCE_INLINE __m128i _mm_castps_si128(__m128);
// FORCE_INLINE __m128i _mm_cmpeq_epi32(__m128i, __m128i);
// FORCE_INLINE __m128i _mm_cvtps_epi32(__m128);
// FORCE_INLINE __m128d _mm_move_sd(__m128d, __m128d);
// FORCE_INLINE __m128i _mm_or_si128(__m128i, __m128i);
// FORCE_INLINE __m128i _mm_set_epi32(int, int, int, int);
// FORCE_INLINE __m128i _mm_set_epi64x(int64_t, int64_t);
// FORCE_INLINE __m128d _mm_set_pd(double, double);
// FORCE_INLINE __m128i _mm_set1_epi32(int);
// FORCE_INLINE __m128i _mm_setzero_si128(void);
// SSE4.1
// FORCE_INLINE __m128d _mm_ceil_pd(__m128d);
// FORCE_INLINE __m128 _mm_ceil_ps(__m128);
// FORCE_INLINE __m128d _mm_floor_pd(__m128d);
// FORCE_INLINE __m128 _mm_floor_ps(__m128);
// FORCE_INLINE __m128d _mm_round_pd(__m128d, int);
// FORCE_INLINE __m128 _mm_round_ps(__m128, int);
// SSE4.2
// FORCE_INLINE uint32_t _mm_crc32_u8(uint32_t, uint8_t);

/* Backwards compatibility for compilers with lack of specific type support */

// Older gcc does not define vld1q_u8_x4 type
#if defined(__GNUC__) && !defined(__clang__) &&                        \
    ((__GNUC__ <= 13 && defined(__arm__)) ||                           \
     (__GNUC__ == 10 && __GNUC_MINOR__ < 3 && defined(__aarch64__)) || \
     (__GNUC__ <= 9 && defined(__aarch64__)))
FORCE_INLINE uint8x16x4_t _sse2rvv_vld1q_u8_x4(const uint8_t *p)
{
    uint8x16x4_t ret;
    ret.val[0] = vld1q_u8(p + 0);
    ret.val[1] = vld1q_u8(p + 16);
    ret.val[2] = vld1q_u8(p + 32);
    ret.val[3] = vld1q_u8(p + 48);
    return ret;
}
#else
// Wraps vld1q_u8_x4
FORCE_INLINE uint8x16x4_t _sse2rvv_vld1q_u8_x4(const uint8_t *p)
{
    return vld1q_u8_x4(p);
}
#endif

#if !defined(__aarch64__) && !defined(_M_ARM64)
/* emulate vaddv u8 variant */
FORCE_INLINE uint8_t _sse2rvv_vaddv_u8(uint8x8_t v8)
{
    const uint64x1_t v1 = vpaddl_u32(vpaddl_u16(vpaddl_u8(v8)));
    return vget_lane_u8(vreinterpret_u8_u64(v1), 0);
}
#else
// Wraps vaddv_u8
FORCE_INLINE uint8_t _sse2rvv_vaddv_u8(uint8x8_t v8)
{
    return vaddv_u8(v8);
}
#endif

#if !defined(__aarch64__) && !defined(_M_ARM64)
/* emulate vaddvq u8 variant */
FORCE_INLINE uint8_t _sse2rvv_vaddvq_u8(uint8x16_t a)
{
    uint8x8_t tmp = vpadd_u8(vget_low_u8(a), vget_high_u8(a));
    uint8_t res = 0;
    for (int i = 0; i < 8; ++i)
        res += tmp[i];
    return res;
}
#else
// Wraps vaddvq_u8
FORCE_INLINE uint8_t _sse2rvv_vaddvq_u8(uint8x16_t a)
{
    return vaddvq_u8(a);
}
#endif

#if !defined(__aarch64__) && !defined(_M_ARM64)
/* emulate vaddvq u16 variant */
FORCE_INLINE uint16_t _sse2rvv_vaddvq_u16(uint16x8_t a)
{
    uint32x4_t m = vpaddlq_u16(a);
    uint64x2_t n = vpaddlq_u32(m);
    uint64x1_t o = vget_low_u64(n) + vget_high_u64(n);

    return vget_lane_u32((uint32x2_t) o, 0);
}
#else
// Wraps vaddvq_u16
FORCE_INLINE uint16_t _sse2rvv_vaddvq_u16(uint16x8_t a)
{
    return vaddvq_u16(a);
}
#endif

/* Function Naming Conventions
 * The naming convention of SSE intrinsics is straightforward. A generic SSE
 * intrinsic function is given as follows:
 *   _mm_<name>_<data_type>
 *
 * The parts of this format are given as follows:
 * 1. <name> describes the operation performed by the intrinsic
 * 2. <data_type> identifies the data type of the function's primary arguments
 *
 * This last part, <data_type>, is a little complicated. It identifies the
 * content of the input values, and can be set to any of the following values:
 * + ps - vectors contain floats (ps stands for packed single-precision)
 * + pd - vectors contain doubles (pd stands for packed double-precision)
 * + epi8/epi16/epi32/epi64 - vectors contain 8-bit/16-bit/32-bit/64-bit
 *                            signed integers
 * + epu8/epu16/epu32/epu64 - vectors contain 8-bit/16-bit/32-bit/64-bit
 *                            unsigned integers
 * + si128 - unspecified 128-bit vector or 256-bit vector
 * + m128/m128i/m128d - identifies input vector types when they are different
 *                      than the type of the returned vector
 *
 * For example, _mm_setzero_ps. The _mm implies that the function returns
 * a 128-bit vector. The _ps at the end implies that the argument vectors
 * contain floats.
 *
 * A complete example: Byte Shuffle - pshufb (_mm_shuffle_epi8)
 *   // Set packed 16-bit integers. 128 bits, 8 short, per 16 bits
 *   __m128i v_in = _mm_setr_epi16(1, 2, 3, 4, 5, 6, 7, 8);
 *   // Set packed 8-bit integers
 *   // 128 bits, 16 chars, per 8 bits
 *   __m128i v_perm = _mm_setr_epi8(1, 0,  2,  3, 8, 9, 10, 11,
 *                                  4, 5, 12, 13, 6, 7, 14, 15);
 *   // Shuffle packed 8-bit integers
 *   __m128i v_out = _mm_shuffle_epi8(v_in, v_perm); // pshufb
 */

/* Constants for use with _mm_prefetch. */
enum _mm_hint {
    _MM_HINT_NTA = 0, /* load data to L1 and L2 cache, mark it as NTA */
    _MM_HINT_T0 = 1,  /* load data to L1 and L2 cache */
    _MM_HINT_T1 = 2,  /* load data to L2 cache only */
    _MM_HINT_T2 = 3,  /* load data to L2 cache only, mark it as NTA */
};

// The bit field mapping to the FPCR(floating-point control register)
typedef struct {
    uint16_t res0;
    uint8_t res1 : 6;
    uint8_t bit22 : 1;
    uint8_t bit23 : 1;
    uint8_t bit24 : 1;
    uint8_t res2 : 7;
#if defined(__aarch64__) || defined(_M_ARM64)
    uint32_t res3;
#endif
} fpcr_bitfield;

// Takes the upper 64 bits of a and places it in the low end of the result
// Takes the lower 64 bits of b and places it into the high end of the result.
// FORCE_INLINE __m128 _mm_shuffle_ps_1032(__m128 a, __m128 b) {}

// takes the lower two 32-bit values from a and swaps them and places in high
// end of result takes the higher two 32 bit values from b and swaps them and
// places in low end of result.
// FORCE_INLINE __m128 _mm_shuffle_ps_2301(__m128 a, __m128 b) {}

// FORCE_INLINE __m128 _mm_shuffle_ps_0321(__m128 a, __m128 b) {}

// FORCE_INLINE __m128 _mm_shuffle_ps_2103(__m128 a, __m128 b) {}

// FORCE_INLINE __m128 _mm_shuffle_ps_1010(__m128 a, __m128 b) {}

// FORCE_INLINE __m128 _mm_shuffle_ps_1001(__m128 a, __m128 b) {}

// FORCE_INLINE __m128 _mm_shuffle_ps_0101(__m128 a, __m128 b) {}

// keeps the low 64 bits of b in the low and puts the high 64 bits of a in the
// high
// FORCE_INLINE __m128 _mm_shuffle_ps_3210(__m128 a, __m128 b) {}

// FORCE_INLINE __m128 _mm_shuffle_ps_0011(__m128 a, __m128 b) {}

// FORCE_INLINE __m128 _mm_shuffle_ps_0022(__m128 a, __m128 b) {}

// FORCE_INLINE __m128 _mm_shuffle_ps_2200(__m128 a, __m128 b) {}

// FORCE_INLINE __m128 _mm_shuffle_ps_3202(__m128 a, __m128 b) {}

// FORCE_INLINE __m128 _mm_shuffle_ps_1133(__m128 a, __m128 b) {}

// FORCE_INLINE __m128 _mm_shuffle_ps_2010(__m128 a, __m128 b) {}

// FORCE_INLINE __m128 _mm_shuffle_ps_2001(__m128 a, __m128 b) {}

// FORCE_INLINE __m128 _mm_shuffle_ps_2032(__m128 a, __m128 b) {}

// C equivalent:
//   __m128i _mm_shuffle_epi32_default(__m128i a,
//                                     __constrange(0, 255) int imm) {
//       __m128i ret;
//       ret[0] = a[imm        & 0x3];   ret[1] = a[(imm >> 2) & 0x3];
//       ret[2] = a[(imm >> 4) & 0x03];  ret[3] = a[(imm >> 6) & 0x03];
//       return ret;
//   }
// #define _mm_shuffle_epi32_default(a, imm)

// Takes the upper 64 bits of a and places it in the low end of the result
// Takes the lower 64 bits of a and places it into the high end of the result.
// FORCE_INLINE __m128i _mm_shuffle_epi_1032(__m128i a) {}

// takes the lower two 32-bit values from a and swaps them and places in low end
// of result takes the higher two 32 bit values from a and swaps them and places
// in high end of result.
// FORCE_INLINE __m128i _mm_shuffle_epi_2301(__m128i a) {}

// rotates the least significant 32 bits into the most significant 32 bits, and
// shifts the rest down
// FORCE_INLINE __m128i _mm_shuffle_epi_0321(__m128i a) {}

// rotates the most significant 32 bits into the least significant 32 bits, and
// shifts the rest up
// FORCE_INLINE __m128i _mm_shuffle_epi_2103(__m128i a) {}

// gets the lower 64 bits of a, and places it in the upper 64 bits
// gets the lower 64 bits of a and places it in the lower 64 bits
// FORCE_INLINE __m128i _mm_shuffle_epi_1010(__m128i a) {}

// gets the lower 64 bits of a, swaps the 0 and 1 elements, and places it in the
// lower 64 bits gets the lower 64 bits of a, and places it in the upper 64 bits
// FORCE_INLINE __m128i _mm_shuffle_epi_1001(__m128i a) {}

// gets the lower 64 bits of a, swaps the 0 and 1 elements and places it in the
// upper 64 bits gets the lower 64 bits of a, swaps the 0 and 1 elements, and
// places it in the lower 64 bits
// FORCE_INLINE __m128i _mm_shuffle_epi_0101(__m128i a) {}

// FORCE_INLINE __m128i _mm_shuffle_epi_2211(__m128i a) {}

// FORCE_INLINE __m128i _mm_shuffle_epi_0122(__m128i a) {}

// FORCE_INLINE __m128i _mm_shuffle_epi_3332(__m128i a) {}

#if defined(__aarch64__) || defined(_M_ARM64)
// #define _mm_shuffle_epi32_splat(a, imm)
#else
// #define _mm_shuffle_epi32_splat(a, imm)
#endif

// NEON does not support a general purpose permute intrinsic.
// Shuffle single-precision (32-bit) floating-point elements in a using the
// control in imm8, and store the results in dst.
//
// C equivalent:
//   __m128 _mm_shuffle_ps_default(__m128 a, __m128 b,
//                                 __constrange(0, 255) int imm) {
//       __m128 ret;
//       ret[0] = a[imm        & 0x3];   ret[1] = a[(imm >> 2) & 0x3];
//       ret[2] = b[(imm >> 4) & 0x03];  ret[3] = b[(imm >> 6) & 0x03];
//       return ret;
//   }
//
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_shuffle_ps
// #define _mm_shuffle_ps_default(a, b, imm)

// Shuffle 16-bit integers in the low 64 bits of a using the control in imm8.
// Store the results in the low 64 bits of dst, with the high 64 bits being
// copied from a to dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_shufflelo_epi16
// #define _mm_shufflelo_epi16_function(a, imm)

// Shuffle 16-bit integers in the high 64 bits of a using the control in imm8.
// Store the results in the high 64 bits of dst, with the low 64 bits being
// copied from a to dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_shufflehi_epi16
// #define _mm_shufflehi_epi16_function(a, imm)

/* MMX */

//_mm_empty is a no-op on arm
// FORCE_INLINE void _mm_empty(void) {} 
/* SSE */

// Add packed single-precision (32-bit) floating-point elements in a and b, and
// store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_add_ps
// FORCE_INLINE __m128 _mm_add_ps(__m128 a, __m128 b) {}

// Add the lower single-precision (32-bit) floating-point element in a and b,
// store the result in the lower element of dst, and copy the upper 3 packed
// elements from a to the upper elements of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_add_ss
// FORCE_INLINE __m128 _mm_add_ss(__m128 a, __m128 b) {}

// Compute the bitwise AND of packed single-precision (32-bit) floating-point
// elements in a and b, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_and_ps
// FORCE_INLINE __m128 _mm_and_ps(__m128 a, __m128 b) {}

// Compute the bitwise NOT of packed single-precision (32-bit) floating-point
// elements in a and then AND with b, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_andnot_ps
// FORCE_INLINE __m128 _mm_andnot_ps(__m128 a, __m128 b) {}

// Average packed unsigned 16-bit integers in a and b, and store the results in
// dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_avg_pu16
// FORCE_INLINE __m64 _mm_avg_pu16(__m64 a, __m64 b) {}

// Average packed unsigned 8-bit integers in a and b, and store the results in
// dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_avg_pu8
// FORCE_INLINE __m64 _mm_avg_pu8(__m64 a, __m64 b) {}

// Compare packed single-precision (32-bit) floating-point elements in a and b
// for equality, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpeq_ps
// FORCE_INLINE __m128 _mm_cmpeq_ps(__m128 a, __m128 b) {}

// Compare the lower single-precision (32-bit) floating-point elements in a and
// b for equality, store the result in the lower element of dst, and copy the
// upper 3 packed elements from a to the upper elements of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpeq_ss
// FORCE_INLINE __m128 _mm_cmpeq_ss(__m128 a, __m128 b) {}

// Compare packed single-precision (32-bit) floating-point elements in a and b
// for greater-than-or-equal, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpge_ps
// FORCE_INLINE __m128 _mm_cmpge_ps(__m128 a, __m128 b) {}

// Compare the lower single-precision (32-bit) floating-point elements in a and
// b for greater-than-or-equal, store the result in the lower element of dst,
// and copy the upper 3 packed elements from a to the upper elements of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpge_ss
// FORCE_INLINE __m128 _mm_cmpge_ss(__m128 a, __m128 b) {}

// Compare packed single-precision (32-bit) floating-point elements in a and b
// for greater-than, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpgt_ps
// FORCE_INLINE __m128 _mm_cmpgt_ps(__m128 a, __m128 b) {}

// Compare the lower single-precision (32-bit) floating-point elements in a and
// b for greater-than, store the result in the lower element of dst, and copy
// the upper 3 packed elements from a to the upper elements of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpgt_ss
// FORCE_INLINE __m128 _mm_cmpgt_ss(__m128 a, __m128 b) {}

// Compare packed single-precision (32-bit) floating-point elements in a and b
// for less-than-or-equal, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmple_ps
// FORCE_INLINE __m128 _mm_cmple_ps(__m128 a, __m128 b) {}

// Compare the lower single-precision (32-bit) floating-point elements in a and
// b for less-than-or-equal, store the result in the lower element of dst, and
// copy the upper 3 packed elements from a to the upper elements of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmple_ss
// FORCE_INLINE __m128 _mm_cmple_ss(__m128 a, __m128 b) {}

// Compare packed single-precision (32-bit) floating-point elements in a and b
// for less-than, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmplt_ps
// FORCE_INLINE __m128 _mm_cmplt_ps(__m128 a, __m128 b) {}

// Compare the lower single-precision (32-bit) floating-point elements in a and
// b for less-than, store the result in the lower element of dst, and copy the
// upper 3 packed elements from a to the upper elements of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmplt_ss
// FORCE_INLINE __m128 _mm_cmplt_ss(__m128 a, __m128 b) {}

// Compare packed single-precision (32-bit) floating-point elements in a and b
// for not-equal, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpneq_ps
// FORCE_INLINE __m128 _mm_cmpneq_ps(__m128 a, __m128 b) {}

// Compare the lower single-precision (32-bit) floating-point elements in a and
// b for not-equal, store the result in the lower element of dst, and copy the
// upper 3 packed elements from a to the upper elements of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpneq_ss
// FORCE_INLINE __m128 _mm_cmpneq_ss(__m128 a, __m128 b) {}

// Compare packed single-precision (32-bit) floating-point elements in a and b
// for not-greater-than-or-equal, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpnge_ps
// FORCE_INLINE __m128 _mm_cmpnge_ps(__m128 a, __m128 b) {}

// Compare the lower single-precision (32-bit) floating-point elements in a and
// b for not-greater-than-or-equal, store the result in the lower element of
// dst, and copy the upper 3 packed elements from a to the upper elements of
// dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpnge_ss
// FORCE_INLINE __m128 _mm_cmpnge_ss(__m128 a, __m128 b) {}

// Compare packed single-precision (32-bit) floating-point elements in a and b
// for not-greater-than, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpngt_ps
// FORCE_INLINE __m128 _mm_cmpngt_ps(__m128 a, __m128 b) {}

// Compare the lower single-precision (32-bit) floating-point elements in a and
// b for not-greater-than, store the result in the lower element of dst, and
// copy the upper 3 packed elements from a to the upper elements of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpngt_ss
// FORCE_INLINE __m128 _mm_cmpngt_ss(__m128 a, __m128 b) {}

// Compare packed single-precision (32-bit) floating-point elements in a and b
// for not-less-than-or-equal, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpnle_ps
// FORCE_INLINE __m128 _mm_cmpnle_ps(__m128 a, __m128 b) {}

// Compare the lower single-precision (32-bit) floating-point elements in a and
// b for not-less-than-or-equal, store the result in the lower element of dst,
// and copy the upper 3 packed elements from a to the upper elements of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpnle_ss
// FORCE_INLINE __m128 _mm_cmpnle_ss(__m128 a, __m128 b) {}

// Compare packed single-precision (32-bit) floating-point elements in a and b
// for not-less-than, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpnlt_ps
// FORCE_INLINE __m128 _mm_cmpnlt_ps(__m128 a, __m128 b) {}

// Compare the lower single-precision (32-bit) floating-point elements in a and
// b for not-less-than, store the result in the lower element of dst, and copy
// the upper 3 packed elements from a to the upper elements of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpnlt_ss
// FORCE_INLINE __m128 _mm_cmpnlt_ss(__m128 a, __m128 b) {}

// Compare packed single-precision (32-bit) floating-point elements in a and b
// to see if neither is NaN, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpord_ps
//
// See also:
// http://stackoverflow.com/questions/8627331/what-does-ordered-unordered-comparison-mean
// http://stackoverflow.com/questions/29349621/neon-isnanval-intrinsics
// FORCE_INLINE __m128 _mm_cmpord_ps(__m128 a, __m128 b) {}

// Compare the lower single-precision (32-bit) floating-point elements in a and
// b to see if neither is NaN, store the result in the lower element of dst, and
// copy the upper 3 packed elements from a to the upper elements of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpord_ss
// FORCE_INLINE __m128 _mm_cmpord_ss(__m128 a, __m128 b) {}

// Compare packed single-precision (32-bit) floating-point elements in a and b
// to see if either is NaN, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpunord_ps
// FORCE_INLINE __m128 _mm_cmpunord_ps(__m128 a, __m128 b) {}

// Compare the lower single-precision (32-bit) floating-point elements in a and
// b to see if either is NaN, store the result in the lower element of dst, and
// copy the upper 3 packed elements from a to the upper elements of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpunord_ss
// FORCE_INLINE __m128 _mm_cmpunord_ss(__m128 a, __m128 b) {}

// Compare the lower single-precision (32-bit) floating-point element in a and b
// for equality, and return the boolean result (0 or 1).
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_comieq_ss
// FORCE_INLINE int _mm_comieq_ss(__m128 a, __m128 b) {}

// Compare the lower single-precision (32-bit) floating-point element in a and b
// for greater-than-or-equal, and return the boolean result (0 or 1).
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_comige_ss
// FORCE_INLINE int _mm_comige_ss(__m128 a, __m128 b) {}

// Compare the lower single-precision (32-bit) floating-point element in a and b
// for greater-than, and return the boolean result (0 or 1).
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_comigt_ss
// FORCE_INLINE int _mm_comigt_ss(__m128 a, __m128 b) {}

// Compare the lower single-precision (32-bit) floating-point element in a and b
// for less-than-or-equal, and return the boolean result (0 or 1).
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_comile_ss
// FORCE_INLINE int _mm_comile_ss(__m128 a, __m128 b) {}

// Compare the lower single-precision (32-bit) floating-point element in a and b
// for less-than, and return the boolean result (0 or 1).
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_comilt_ss
// FORCE_INLINE int _mm_comilt_ss(__m128 a, __m128 b) {}

// Compare the lower single-precision (32-bit) floating-point element in a and b
// for not-equal, and return the boolean result (0 or 1).
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_comineq_ss
// FORCE_INLINE int _mm_comineq_ss(__m128 a, __m128 b) {}

// Convert packed signed 32-bit integers in b to packed single-precision
// (32-bit) floating-point elements, store the results in the lower 2 elements
// of dst, and copy the upper 2 packed elements from a to the upper elements of
// dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvt_pi2ps
// FORCE_INLINE __m128 _mm_cvt_pi2ps(__m128 a, __m64 b) {}

// Convert packed single-precision (32-bit) floating-point elements in a to
// packed 32-bit integers, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvt_ps2pi
// FORCE_INLINE __m64 _mm_cvt_ps2pi(__m128 a) {}

// Convert the signed 32-bit integer b to a single-precision (32-bit)
// floating-point element, store the result in the lower element of dst, and
// copy the upper 3 packed elements from a to the upper elements of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvt_si2ss
// FORCE_INLINE __m128 _mm_cvt_si2ss(__m128 a, int b) {}

// Convert the lower single-precision (32-bit) floating-point element in a to a
// 32-bit integer, and store the result in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvt_ss2si
// FORCE_INLINE int _mm_cvt_ss2si(__m128 a) {}

// Convert packed 16-bit integers in a to packed single-precision (32-bit)
// floating-point elements, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtpi16_ps
// FORCE_INLINE __m128 _mm_cvtpi16_ps(__m64 a) {}

// Convert packed 32-bit integers in b to packed single-precision (32-bit)
// floating-point elements, store the results in the lower 2 elements of dst,
// and copy the upper 2 packed elements from a to the upper elements of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtpi32_ps
// FORCE_INLINE __m128 _mm_cvtpi32_ps(__m128 a, __m64 b) {}

// Convert packed signed 32-bit integers in a to packed single-precision
// (32-bit) floating-point elements, store the results in the lower 2 elements
// of dst, then convert the packed signed 32-bit integers in b to
// single-precision (32-bit) floating-point element, and store the results in
// the upper 2 elements of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtpi32x2_ps
// FORCE_INLINE __m128 _mm_cvtpi32x2_ps(__m64 a, __m64 b) {}

// Convert the lower packed 8-bit integers in a to packed single-precision
// (32-bit) floating-point elements, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtpi8_ps
// FORCE_INLINE __m128 _mm_cvtpi8_ps(__m64 a) {}

// Convert packed single-precision (32-bit) floating-point elements in a to
// packed 16-bit integers, and store the results in dst. Note: this intrinsic
// will generate 0x7FFF, rather than 0x8000, for input values between 0x7FFF and
// 0x7FFFFFFF.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtps_pi16
// FORCE_INLINE __m64 _mm_cvtps_pi16(__m128 a) {}

// Convert packed single-precision (32-bit) floating-point elements in a to
// packed 32-bit integers, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtps_pi32
// #define _mm_cvtps_pi32(a) _mm_cvt_ps2pi(a)

// Convert packed single-precision (32-bit) floating-point elements in a to
// packed 8-bit integers, and store the results in lower 4 elements of dst.
// Note: this intrinsic will generate 0x7F, rather than 0x80, for input values
// between 0x7F and 0x7FFFFFFF.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtps_pi8
// FORCE_INLINE __m64 _mm_cvtps_pi8(__m128 a) {}

// Convert packed unsigned 16-bit integers in a to packed single-precision
// (32-bit) floating-point elements, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtpu16_ps
// FORCE_INLINE __m128 _mm_cvtpu16_ps(__m64 a) {}

// Convert the lower packed unsigned 8-bit integers in a to packed
// single-precision (32-bit) floating-point elements, and store the results in
// dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtpu8_ps
// FORCE_INLINE __m128 _mm_cvtpu8_ps(__m64 a) {}

// Convert the signed 32-bit integer b to a single-precision (32-bit)
// floating-point element, store the result in the lower element of dst, and
// copy the upper 3 packed elements from a to the upper elements of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtsi32_ss
// #define _mm_cvtsi32_ss(a, b) _mm_cvt_si2ss(a, b)

// Convert the signed 64-bit integer b to a single-precision (32-bit)
// floating-point element, store the result in the lower element of dst, and
// copy the upper 3 packed elements from a to the upper elements of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtsi64_ss
// FORCE_INLINE __m128 _mm_cvtsi64_ss(__m128 a, int64_t b) {}

// Copy the lower single-precision (32-bit) floating-point element of a to dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtss_f32
// FORCE_INLINE float _mm_cvtss_f32(__m128 a) {}

// Convert the lower single-precision (32-bit) floating-point element in a to a
// 32-bit integer, and store the result in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtss_si32
// #define _mm_cvtss_si32(a) _mm_cvt_ss2si(a)

// Convert the lower single-precision (32-bit) floating-point element in a to a
// 64-bit integer, and store the result in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtss_si64
// FORCE_INLINE int64_t _mm_cvtss_si64(__m128 a) {}

// Convert packed single-precision (32-bit) floating-point elements in a to
// packed 32-bit integers with truncation, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtt_ps2pi
// FORCE_INLINE __m64 _mm_cvtt_ps2pi(__m128 a) {}

// Convert the lower single-precision (32-bit) floating-point element in a to a
// 32-bit integer with truncation, and store the result in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtt_ss2si
// FORCE_INLINE int _mm_cvtt_ss2si(__m128 a) {}

// Convert packed single-precision (32-bit) floating-point elements in a to
// packed 32-bit integers with truncation, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvttps_pi32
// #define _mm_cvttps_pi32(a) _mm_cvtt_ps2pi(a)

// Convert the lower single-precision (32-bit) floating-point element in a to a
// 32-bit integer with truncation, and store the result in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvttss_si32
// #define _mm_cvttss_si32(a) _mm_cvtt_ss2si(a)

// Convert the lower single-precision (32-bit) floating-point element in a to a
// 64-bit integer with truncation, and store the result in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvttss_si64
// FORCE_INLINE int64_t _mm_cvttss_si64(__m128 a) {}

// Divide packed single-precision (32-bit) floating-point elements in a by
// packed elements in b, and store the results in dst.
// Due to ARMv7-A NEON's lack of a precise division intrinsic, we implement
// division by multiplying a by b's reciprocal before using the Newton-Raphson
// method to approximate the results.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_div_ps
// FORCE_INLINE __m128 _mm_div_ps(__m128 a, __m128 b) {}

// Divide the lower single-precision (32-bit) floating-point element in a by the
// lower single-precision (32-bit) floating-point element in b, store the result
// in the lower element of dst, and copy the upper 3 packed elements from a to
// the upper elements of dst.
// Warning: ARMv7-A does not produce the same result compared to Intel and not
// IEEE-compliant.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_div_ss
// FORCE_INLINE __m128 _mm_div_ss(__m128 a, __m128 b) {}

// Extract a 16-bit integer from a, selected with imm8, and store the result in
// the lower element of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_extract_pi16
// #define _mm_extract_pi16(a, imm)

// Free aligned memory that was allocated with _mm_malloc.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_free
#if !defined(SSE2RVV_ALLOC_DEFINED)
// FORCE_INLINE void _mm_free(void *addr) {}
#endif

FORCE_INLINE uint64_t _sse2rvv_get_fpcr(void)
{
    uint64_t value;
#if defined(_MSC_VER)
    value = _ReadStatusReg(ARM64_FPCR);
#else
    __asm__ __volatile__("mrs %0, FPCR" : "=r"(value)); /* read */
#endif
    return value;
}

FORCE_INLINE void _sse2rvv_set_fpcr(uint64_t value)
{
#if defined(_MSC_VER)
    _WriteStatusReg(ARM64_FPCR, value);
#else
    __asm__ __volatile__("msr FPCR, %0" ::"r"(value)); /* write */
#endif
}

// Macro: Get the flush zero bits from the MXCSR control and status register.
// The flush zero may contain any of the following flags: _MM_FLUSH_ZERO_ON or
// _MM_FLUSH_ZERO_OFF
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_MM_GET_FLUSH_ZERO_MODE
// FORCE_INLINE unsigned int _sse2rvv_mm_get_flush_zero_mode(void) {}

// Macro: Get the rounding mode bits from the MXCSR control and status register.
// The rounding mode may contain any of the following flags: _MM_ROUND_NEAREST,
// _MM_ROUND_DOWN, _MM_ROUND_UP, _MM_ROUND_TOWARD_ZERO
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_MM_GET_ROUNDING_MODE
FORCE_INLINE unsigned int _MM_GET_ROUNDING_MODE(void)
{
    union {
        fpcr_bitfield field;
#if defined(__aarch64__) || defined(_M_ARM64)
        uint64_t value;
#else
        uint32_t value;
#endif
    } r;

#if defined(__aarch64__) || defined(_M_ARM64)
    r.value = _sse2rvv_get_fpcr();
#else
    __asm__ __volatile__("vmrs %0, FPSCR" : "=r"(r.value)); /* read */
#endif

    if (r.field.bit22) {
        return r.field.bit23 ? _MM_ROUND_TOWARD_ZERO : _MM_ROUND_UP;
    } else {
        return r.field.bit23 ? _MM_ROUND_DOWN : _MM_ROUND_NEAREST;
    }
}

// Copy a to dst, and insert the 16-bit integer i into dst at the location
// specified by imm8.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_insert_pi16
// #define _mm_insert_pi16(a, b, imm)

// Load 128-bits (composed of 4 packed single-precision (32-bit) floating-point
// elements) from memory into dst. mem_addr must be aligned on a 16-byte
// boundary or a general-protection exception may be generated.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_load_ps
// FORCE_INLINE __m128 _mm_load_ps(const float *p) {}

// Load a single-precision (32-bit) floating-point element from memory into all
// elements of dst.
//
//   dst[31:0] := MEM[mem_addr+31:mem_addr]
//   dst[63:32] := MEM[mem_addr+31:mem_addr]
//   dst[95:64] := MEM[mem_addr+31:mem_addr]
//   dst[127:96] := MEM[mem_addr+31:mem_addr]
//
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_load_ps1
// #define _mm_load_ps1 _mm_load1_ps

// Load a single-precision (32-bit) floating-point element from memory into the
// lower of dst, and zero the upper 3 elements. mem_addr does not need to be
// aligned on any particular boundary.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_load_ss
// FORCE_INLINE __m128 _mm_load_ss(const float *p) {}

// Load a single-precision (32-bit) floating-point element from memory into all
// elements of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_load1_ps
// FORCE_INLINE __m128 _mm_load1_ps(const float *p) {}

// Load 2 single-precision (32-bit) floating-point elements from memory into the
// upper 2 elements of dst, and copy the lower 2 elements from a to dst.
// mem_addr does not need to be aligned on any particular boundary.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_loadh_pi
// FORCE_INLINE __m128 _mm_loadh_pi(__m128 a, __m64 const *p) {}

// Load 2 single-precision (32-bit) floating-point elements from memory into the
// lower 2 elements of dst, and copy the upper 2 elements from a to dst.
// mem_addr does not need to be aligned on any particular boundary.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_loadl_pi
// FORCE_INLINE __m128 _mm_loadl_pi(__m128 a, __m64 const *p) {}

// Load 4 single-precision (32-bit) floating-point elements from memory into dst
// in reverse order. mem_addr must be aligned on a 16-byte boundary or a
// general-protection exception may be generated.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_loadr_ps
// FORCE_INLINE __m128 _mm_loadr_ps(const float *p) {}

// Load 128-bits (composed of 4 packed single-precision (32-bit) floating-point
// elements) from memory into dst. mem_addr does not need to be aligned on any
// particular boundary.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_loadu_ps
// FORCE_INLINE __m128 _mm_loadu_ps(const float *p) {}

// Load unaligned 16-bit integer from memory into the first element of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_loadu_si16
// FORCE_INLINE __m128i _mm_loadu_si16(const void *p) {}

// Load unaligned 64-bit integer from memory into the first element of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_loadu_si64
// FORCE_INLINE __m128i _mm_loadu_si64(const void *p) {}

// Allocate size bytes of memory, aligned to the alignment specified in align,
// and return a pointer to the allocated memory. _mm_free should be used to free
// memory that is allocated with _mm_malloc.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_malloc
#if !defined(SSE2RVV_ALLOC_DEFINED)
// FORCE_INLINE void *_mm_malloc(size_t size, size_t align) {}
#endif

// Conditionally store 8-bit integer elements from a into memory using mask
// (elements are not stored when the highest bit is not set in the corresponding
// element) and a non-temporal memory hint.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_maskmove_si64
// FORCE_INLINE void _mm_maskmove_si64(__m64 a, __m64 mask, char *mem_addr) {}

// Conditionally store 8-bit integer elements from a into memory using mask
// (elements are not stored when the highest bit is not set in the corresponding
// element) and a non-temporal memory hint.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_m_maskmovq
// #define _m_maskmovq(a, mask, mem_addr) _mm_maskmove_si64(a, mask, mem_addr)

// Compare packed signed 16-bit integers in a and b, and store packed maximum
// values in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_max_pi16
// FORCE_INLINE __m64 _mm_max_pi16(__m64 a, __m64 b) {}

// Compare packed single-precision (32-bit) floating-point elements in a and b,
// and store packed maximum values in dst. dst does not follow the IEEE Standard
// for Floating-Point Arithmetic (IEEE 754) maximum value when inputs are NaN or
// signed-zero values.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_max_ps
// FORCE_INLINE __m128 _mm_max_ps(__m128 a, __m128 b) {}

// Compare packed unsigned 8-bit integers in a and b, and store packed maximum
// values in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_max_pu8
// FORCE_INLINE __m64 _mm_max_pu8(__m64 a, __m64 b) {}

// Compare the lower single-precision (32-bit) floating-point elements in a and
// b, store the maximum value in the lower element of dst, and copy the upper 3
// packed elements from a to the upper element of dst. dst does not follow the
// IEEE Standard for Floating-Point Arithmetic (IEEE 754) maximum value when
// inputs are NaN or signed-zero values.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_max_ss
// FORCE_INLINE __m128 _mm_max_ss(__m128 a, __m128 b) {}

// Compare packed signed 16-bit integers in a and b, and store packed minimum
// values in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_min_pi16
// FORCE_INLINE __m64 _mm_min_pi16(__m64 a, __m64 b) {}

// Compare packed single-precision (32-bit) floating-point elements in a and b,
// and store packed minimum values in dst. dst does not follow the IEEE Standard
// for Floating-Point Arithmetic (IEEE 754) minimum value when inputs are NaN or
// signed-zero values.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_min_ps
// FORCE_INLINE __m128 _mm_min_ps(__m128 a, __m128 b) {}

// Compare packed unsigned 8-bit integers in a and b, and store packed minimum
// values in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_min_pu8
// FORCE_INLINE __m64 _mm_min_pu8(__m64 a, __m64 b) {}

// Compare the lower single-precision (32-bit) floating-point elements in a and
// b, store the minimum value in the lower element of dst, and copy the upper 3
// packed elements from a to the upper element of dst. dst does not follow the
// IEEE Standard for Floating-Point Arithmetic (IEEE 754) minimum value when
// inputs are NaN or signed-zero values.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_min_ss
// FORCE_INLINE __m128 _mm_min_ss(__m128 a, __m128 b) {}

// Move the lower single-precision (32-bit) floating-point element from b to the
// lower element of dst, and copy the upper 3 packed elements from a to the
// upper elements of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_move_ss
// FORCE_INLINE __m128 _mm_move_ss(__m128 a, __m128 b) {}

// Move the upper 2 single-precision (32-bit) floating-point elements from b to
// the lower 2 elements of dst, and copy the upper 2 elements from a to the
// upper 2 elements of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_movehl_ps
// FORCE_INLINE __m128 _mm_movehl_ps(__m128 a, __m128 b) {}

// Move the lower 2 single-precision (32-bit) floating-point elements from b to
// the upper 2 elements of dst, and copy the lower 2 elements from a to the
// lower 2 elements of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_movelh_ps
// FORCE_INLINE __m128 _mm_movelh_ps(__m128 __A, __m128 __B) {}

// Create mask from the most significant bit of each 8-bit element in a, and
// store the result in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_movemask_pi8
// FORCE_INLINE int _mm_movemask_pi8(__m64 a) {}

// Set each bit of mask dst based on the most significant bit of the
// corresponding packed single-precision (32-bit) floating-point element in a.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_movemask_ps
// FORCE_INLINE int _mm_movemask_ps(__m128 a) {}

// Multiply packed single-precision (32-bit) floating-point elements in a and b,
// and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_mul_ps
// FORCE_INLINE __m128 _mm_mul_ps(__m128 a, __m128 b) {}

// Multiply the lower single-precision (32-bit) floating-point element in a and
// b, store the result in the lower element of dst, and copy the upper 3 packed
// elements from a to the upper elements of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_mul_ss
// FORCE_INLINE __m128 _mm_mul_ss(__m128 a, __m128 b) {}

// Multiply the packed unsigned 16-bit integers in a and b, producing
// intermediate 32-bit integers, and store the high 16 bits of the intermediate
// integers in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_mulhi_pu16
// FORCE_INLINE __m64 _mm_mulhi_pu16(__m64 a, __m64 b) {}

// Compute the bitwise OR of packed single-precision (32-bit) floating-point
// elements in a and b, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_or_ps
// FORCE_INLINE __m128 _mm_or_ps(__m128 a, __m128 b) {}

// Average packed unsigned 8-bit integers in a and b, and store the results in
// dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_m_pavgb
// #define _m_pavgb(a, b) _mm_avg_pu8(a, b)

// Average packed unsigned 16-bit integers in a and b, and store the results in
// dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_m_pavgw
// #define _m_pavgw(a, b) _mm_avg_pu16(a, b)

// Extract a 16-bit integer from a, selected with imm8, and store the result in
// the lower element of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_m_pextrw
// #define _m_pextrw(a, imm) _mm_extract_pi16(a, imm)

// Copy a to dst, and insert the 16-bit integer i into dst at the location
// specified by imm8.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=m_pinsrw
// #define _m_pinsrw(a, i, imm) _mm_insert_pi16(a, i, imm)

// Compare packed signed 16-bit integers in a and b, and store packed maximum
// values in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_m_pmaxsw
// #define _m_pmaxsw(a, b) _mm_max_pi16(a, b)

// Compare packed unsigned 8-bit integers in a and b, and store packed maximum
// values in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_m_pmaxub
// #define _m_pmaxub(a, b) _mm_max_pu8(a, b)

// Compare packed signed 16-bit integers in a and b, and store packed minimum
// values in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_m_pminsw
// #define _m_pminsw(a, b) _mm_min_pi16(a, b)

// Compare packed unsigned 8-bit integers in a and b, and store packed minimum
// values in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_m_pminub
// #define _m_pminub(a, b) _mm_min_pu8(a, b)

// Create mask from the most significant bit of each 8-bit element in a, and
// store the result in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_m_pmovmskb
// #define _m_pmovmskb(a) _mm_movemask_pi8(a)

// Multiply the packed unsigned 16-bit integers in a and b, producing
// intermediate 32-bit integers, and store the high 16 bits of the intermediate
// integers in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_m_pmulhuw
// #define _m_pmulhuw(a, b) _mm_mulhi_pu16(a, b)

// Fetch the line of data from memory that contains address p to a location in
// the cache hierarchy specified by the locality hint i.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_prefetch
// FORCE_INLINE void _mm_prefetch(char const *p, int i) {}

// Compute the absolute differences of packed unsigned 8-bit integers in a and
// b, then horizontally sum each consecutive 8 differences to produce four
// unsigned 16-bit integers, and pack these unsigned 16-bit integers in the low
// 16 bits of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=m_psadbw
// #define _m_psadbw(a, b) _mm_sad_pu8(a, b)

// Shuffle 16-bit integers in a using the control in imm8, and store the results
// in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_m_pshufw
// #define _m_pshufw(a, imm) _mm_shuffle_pi16(a, imm)

// Compute the approximate reciprocal of packed single-precision (32-bit)
// floating-point elements in a, and store the results in dst. The maximum
// relative error for this approximation is less than 1.5*2^-12.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_rcp_ps
// FORCE_INLINE __m128 _mm_rcp_ps(__m128 in) {}

// Compute the approximate reciprocal of the lower single-precision (32-bit)
// floating-point element in a, store the result in the lower element of dst,
// and copy the upper 3 packed elements from a to the upper elements of dst. The
// maximum relative error for this approximation is less than 1.5*2^-12.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_rcp_ss
// FORCE_INLINE __m128 _mm_rcp_ss(__m128 a) {}

// Compute the approximate reciprocal square root of packed single-precision
// (32-bit) floating-point elements in a, and store the results in dst. The
// maximum relative error for this approximation is less than 1.5*2^-12.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_rsqrt_ps
// FORCE_INLINE __m128 _mm_rsqrt_ps(__m128 in) {}

// Compute the approximate reciprocal square root of the lower single-precision
// (32-bit) floating-point element in a, store the result in the lower element
// of dst, and copy the upper 3 packed elements from a to the upper elements of
// dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_rsqrt_ss
// FORCE_INLINE __m128 _mm_rsqrt_ss(__m128 in) {}

// Compute the absolute differences of packed unsigned 8-bit integers in a and
// b, then horizontally sum each consecutive 8 differences to produce four
// unsigned 16-bit integers, and pack these unsigned 16-bit integers in the low
// 16 bits of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_sad_pu8
// FORCE_INLINE __m64 _mm_sad_pu8(__m64 a, __m64 b) {}

// Macro: Set the flush zero bits of the MXCSR control and status register to
// the value in unsigned 32-bit integer a. The flush zero may contain any of the
// following flags: _MM_FLUSH_ZERO_ON or _MM_FLUSH_ZERO_OFF
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_MM_SET_FLUSH_ZERO_MODE
// FORCE_INLINE void _sse2rvv_mm_set_flush_zero_mode(unsigned int flag) {}

// Set packed single-precision (32-bit) floating-point elements in dst with the
// supplied values.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_set_ps
// FORCE_INLINE __m128 _mm_set_ps(float w, float z, float y, float x) {}

// Broadcast single-precision (32-bit) floating-point value a to all elements of
// dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_set_ps1
// FORCE_INLINE __m128 _mm_set_ps1(float _w) {}

// Macro: Set the rounding mode bits of the MXCSR control and status register to
// the value in unsigned 32-bit integer a. The rounding mode may contain any of
// the following flags: _MM_ROUND_NEAREST, _MM_ROUND_DOWN, _MM_ROUND_UP,
// _MM_ROUND_TOWARD_ZERO
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_MM_SET_ROUNDING_MODE
FORCE_INLINE void _MM_SET_ROUNDING_MODE(int rounding)
{
    union {
        fpcr_bitfield field;
#if defined(__aarch64__) || defined(_M_ARM64)
        uint64_t value;
#else
        uint32_t value;
#endif
    } r;

#if defined(__aarch64__) || defined(_M_ARM64)
    r.value = _sse2rvv_get_fpcr();
#else
    __asm__ __volatile__("vmrs %0, FPSCR" : "=r"(r.value)); /* read */
#endif

    switch (rounding) {
    case _MM_ROUND_TOWARD_ZERO:
        r.field.bit22 = 1;
        r.field.bit23 = 1;
        break;
    case _MM_ROUND_DOWN:
        r.field.bit22 = 0;
        r.field.bit23 = 1;
        break;
    case _MM_ROUND_UP:
        r.field.bit22 = 1;
        r.field.bit23 = 0;
        break;
    default:  //_MM_ROUND_NEAREST
        r.field.bit22 = 0;
        r.field.bit23 = 0;
    }

#if defined(__aarch64__) || defined(_M_ARM64)
    _sse2rvv_set_fpcr(r.value);
#else
    __asm__ __volatile__("vmsr FPSCR, %0" ::"r"(r)); /* write */
#endif
}

// Copy single-precision (32-bit) floating-point element a to the lower element
// of dst, and zero the upper 3 elements.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_set_ss
// FORCE_INLINE __m128 _mm_set_ss(float a) {}

// Broadcast single-precision (32-bit) floating-point value a to all elements of
// dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_set1_ps
// FORCE_INLINE __m128 _mm_set1_ps(float _w) {}

// Set the MXCSR control and status register with the value in unsigned 32-bit
// integer a.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_setcsr
// FIXME: _mm_setcsr() implementation supports changing the rounding mode only.
// FORCE_INLINE void _mm_setcsr(unsigned int a) {}

// Get the unsigned 32-bit value of the MXCSR control and status register.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_getcsr
// FIXME: _mm_getcsr() implementation supports reading the rounding mode only.
// FORCE_INLINE unsigned int _mm_getcsr(void) {}

// Set packed single-precision (32-bit) floating-point elements in dst with the
// supplied values in reverse order.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_setr_ps
// FORCE_INLINE __m128 _mm_setr_ps(float w, float z, float y, float x) {}

// Return vector of type __m128 with all elements set to zero.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_setzero_ps
// FORCE_INLINE __m128 _mm_setzero_ps(void) {}

// Shuffle 16-bit integers in a using the control in imm8, and store the results
// in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_shuffle_pi16
#ifdef _sse2rvv_shuffle
// #define _mm_shuffle_pi16(a, imm)
#else
// #define _mm_shuffle_pi16(a, imm)
#endif

// Perform a serializing operation on all store-to-memory instructions that were
// issued prior to this instruction. Guarantees that every store instruction
// that precedes, in program order, is globally visible before any store
// instruction which follows the fence in program order.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_sfence
// FORCE_INLINE void _mm_sfence(void) {}

// Perform a serializing operation on all load-from-memory and store-to-memory
// instructions that were issued prior to this instruction. Guarantees that
// every memory access that precedes, in program order, the memory fence
// instruction is globally visible before any memory instruction which follows
// the fence in program order.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_mfence
// FORCE_INLINE void _mm_mfence(void) {}

// Perform a serializing operation on all load-from-memory instructions that
// were issued prior to this instruction. Guarantees that every load instruction
// that precedes, in program order, is globally visible before any load
// instruction which follows the fence in program order.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_lfence
// FORCE_INLINE void _mm_lfence(void) {}

// FORCE_INLINE __m128 _mm_shuffle_ps(__m128 a, __m128 b, __constrange(0,255)
// int imm)
#ifdef _sse2rvv_shuffle
// #define _mm_shuffle_ps(a, b, imm)
#else  // generic
// #define _mm_shuffle_ps(a, b, imm)
#endif

// Compute the square root of packed single-precision (32-bit) floating-point
// elements in a, and store the results in dst.
// Due to ARMv7-A NEON's lack of a precise square root intrinsic, we implement
// square root by multiplying input in with its reciprocal square root before
// using the Newton-Raphson method to approximate the results.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_sqrt_ps
// FORCE_INLINE __m128 _mm_sqrt_ps(__m128 in) {}

// Compute the square root of the lower single-precision (32-bit) floating-point
// element in a, store the result in the lower element of dst, and copy the
// upper 3 packed elements from a to the upper elements of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_sqrt_ss
// FORCE_INLINE __m128 _mm_sqrt_ss(__m128 in) {}

// Store 128-bits (composed of 4 packed single-precision (32-bit) floating-point
// elements) from a into memory. mem_addr must be aligned on a 16-byte boundary
// or a general-protection exception may be generated.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_store_ps
// FORCE_INLINE void _mm_store_ps(float *p, __m128 a) {}

// Store the lower single-precision (32-bit) floating-point element from a into
// 4 contiguous elements in memory. mem_addr must be aligned on a 16-byte
// boundary or a general-protection exception may be generated.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_store_ps1
// FORCE_INLINE void _mm_store_ps1(float *p, __m128 a) {}

// Store the lower single-precision (32-bit) floating-point element from a into
// memory. mem_addr does not need to be aligned on any particular boundary.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_store_ss
// FORCE_INLINE void _mm_store_ss(float *p, __m128 a) {}

// Store the lower single-precision (32-bit) floating-point element from a into
// 4 contiguous elements in memory. mem_addr must be aligned on a 16-byte
// boundary or a general-protection exception may be generated.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_store1_ps
// #define _mm_store1_ps _mm_store_ps1

// Store the upper 2 single-precision (32-bit) floating-point elements from a
// into memory.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_storeh_pi
// FORCE_INLINE void _mm_storeh_pi(__m64 *p, __m128 a) {}

// Store the lower 2 single-precision (32-bit) floating-point elements from a
// into memory.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_storel_pi
// FORCE_INLINE void _mm_storel_pi(__m64 *p, __m128 a) {}

// Store 4 single-precision (32-bit) floating-point elements from a into memory
// in reverse order. mem_addr must be aligned on a 16-byte boundary or a
// general-protection exception may be generated.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_storer_ps
// FORCE_INLINE void _mm_storer_ps(float *p, __m128 a) {}

// Store 128-bits (composed of 4 packed single-precision (32-bit) floating-point
// elements) from a into memory. mem_addr does not need to be aligned on any
// particular boundary.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_storeu_ps
// FORCE_INLINE void _mm_storeu_ps(float *p, __m128 a) {}

// Stores 16-bits of integer data a at the address p.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_storeu_si16
// FORCE_INLINE void _mm_storeu_si16(void *p, __m128i a) {}

// Stores 64-bits of integer data a at the address p.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_storeu_si64
// FORCE_INLINE void _mm_storeu_si64(void *p, __m128i a) {}

// Store 64-bits of integer data from a into memory using a non-temporal memory
// hint.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_stream_pi
// FORCE_INLINE void _mm_stream_pi(__m64 *p, __m64 a) {}

// Store 128-bits (composed of 4 packed single-precision (32-bit) floating-
// point elements) from a into memory using a non-temporal memory hint.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_stream_ps
// FORCE_INLINE void _mm_stream_ps(float *p, __m128 a) {}

// Subtract packed single-precision (32-bit) floating-point elements in b from
// packed single-precision (32-bit) floating-point elements in a, and store the
// results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_sub_ps
// FORCE_INLINE __m128 _mm_sub_ps(__m128 a, __m128 b) {}

// Subtract the lower single-precision (32-bit) floating-point element in b from
// the lower single-precision (32-bit) floating-point element in a, store the
// result in the lower element of dst, and copy the upper 3 packed elements from
// a to the upper elements of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_sub_ss
// FORCE_INLINE __m128 _mm_sub_ss(__m128 a, __m128 b) {}

// Macro: Transpose the 4x4 matrix formed by the 4 rows of single-precision
// (32-bit) floating-point elements in row0, row1, row2, and row3, and store the
// transposed matrix in these vectors (row0 now contains column 0, etc.).
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=MM_TRANSPOSE4_PS
#define _MM_TRANSPOSE4_PS(row0, row1, row2, row3)         \
    do {                                                  \
        float32x4x2_t ROW01 = vtrnq_f32(row0, row1);      \
        float32x4x2_t ROW23 = vtrnq_f32(row2, row3);      \
        row0 = vcombine_f32(vget_low_f32(ROW01.val[0]),   \
                            vget_low_f32(ROW23.val[0]));  \
        row1 = vcombine_f32(vget_low_f32(ROW01.val[1]),   \
                            vget_low_f32(ROW23.val[1]));  \
        row2 = vcombine_f32(vget_high_f32(ROW01.val[0]),  \
                            vget_high_f32(ROW23.val[0])); \
        row3 = vcombine_f32(vget_high_f32(ROW01.val[1]),  \
                            vget_high_f32(ROW23.val[1])); \
    } while (0)

// according to the documentation, these intrinsics behave the same as the
// non-'u' versions.  We'll just alias them here.
// #define _mm_ucomieq_ss _mm_comieq_ss
// #define _mm_ucomige_ss _mm_comige_ss
// #define _mm_ucomigt_ss _mm_comigt_ss
// #define _mm_ucomile_ss _mm_comile_ss
// #define _mm_ucomilt_ss _mm_comilt_ss
// #define _mm_ucomineq_ss _mm_comineq_ss

// Return vector of type __m128i with undefined elements.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=mm_undefined_si128
// FORCE_INLINE __m128i _mm_undefined_si128(void) {}

// Return vector of type __m128 with undefined elements.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_undefined_ps
// FORCE_INLINE __m128 _mm_undefined_ps(void) {}

// Unpack and interleave single-precision (32-bit) floating-point elements from
// the high half a and b, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_unpackhi_ps
// FORCE_INLINE __m128 _mm_unpackhi_ps(__m128 a, __m128 b) {}

// Unpack and interleave single-precision (32-bit) floating-point elements from
// the low half of a and b, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_unpacklo_ps
// FORCE_INLINE __m128 _mm_unpacklo_ps(__m128 a, __m128 b) {}

// Compute the bitwise XOR of packed single-precision (32-bit) floating-point
// elements in a and b, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_xor_ps
// FORCE_INLINE __m128 _mm_xor_ps(__m128 a, __m128 b) {}

/* SSE2 */

// Add packed 16-bit integers in a and b, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_add_epi16
// FORCE_INLINE __m128i _mm_add_epi16(__m128i a, __m128i b) {}

// Add packed 32-bit integers in a and b, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_add_epi32
// FORCE_INLINE __m128i _mm_add_epi32(__m128i a, __m128i b) {}

// Add packed 64-bit integers in a and b, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_add_epi64
// FORCE_INLINE __m128i _mm_add_epi64(__m128i a, __m128i b) {}

// Add packed 8-bit integers in a and b, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_add_epi8
// FORCE_INLINE __m128i _mm_add_epi8(__m128i a, __m128i b) {}

// Add packed double-precision (64-bit) floating-point elements in a and b, and
// store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_add_pd
// FORCE_INLINE __m128d _mm_add_pd(__m128d a, __m128d b) {}

// Add the lower double-precision (64-bit) floating-point element in a and b,
// store the result in the lower element of dst, and copy the upper element from
// a to the upper element of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_add_sd
// FORCE_INLINE __m128d _mm_add_sd(__m128d a, __m128d b) {}

// Add 64-bit integers a and b, and store the result in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_add_si64
// FORCE_INLINE __m64 _mm_add_si64(__m64 a, __m64 b) {}

// Add packed signed 16-bit integers in a and b using saturation, and store the
// results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_adds_epi16
// FORCE_INLINE __m128i _mm_adds_epi16(__m128i a, __m128i b) {}

// Add packed signed 8-bit integers in a and b using saturation, and store the
// results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_adds_epi8
// FORCE_INLINE __m128i _mm_adds_epi8(__m128i a, __m128i b) {}

// Add packed unsigned 16-bit integers in a and b using saturation, and store
// the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_adds_epu16
// FORCE_INLINE __m128i _mm_adds_epu16(__m128i a, __m128i b) {}

// Add packed unsigned 8-bit integers in a and b using saturation, and store the
// results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_adds_epu8
// FORCE_INLINE __m128i _mm_adds_epu8(__m128i a, __m128i b) {}

// Compute the bitwise AND of packed double-precision (64-bit) floating-point
// elements in a and b, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_and_pd
// FORCE_INLINE __m128d _mm_and_pd(__m128d a, __m128d b) {}

// Compute the bitwise AND of 128 bits (representing integer data) in a and b,
// and store the result in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_and_si128
// FORCE_INLINE __m128i _mm_and_si128(__m128i a, __m128i b) {}

// Compute the bitwise NOT of packed double-precision (64-bit) floating-point
// elements in a and then AND with b, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_andnot_pd
// FORCE_INLINE __m128d _mm_andnot_pd(__m128d a, __m128d b) {}

// Compute the bitwise NOT of 128 bits (representing integer data) in a and then
// AND with b, and store the result in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_andnot_si128
// FORCE_INLINE __m128i _mm_andnot_si128(__m128i a, __m128i b) {}

// Average packed unsigned 16-bit integers in a and b, and store the results in
// dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_avg_epu16
// FORCE_INLINE __m128i _mm_avg_epu16(__m128i a, __m128i b) {}

// Average packed unsigned 8-bit integers in a and b, and store the results in
// dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_avg_epu8
// FORCE_INLINE __m128i _mm_avg_epu8(__m128i a, __m128i b) {}

// Shift a left by imm8 bytes while shifting in zeros, and store the results in
// dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_bslli_si128
// #define _mm_bslli_si128(a, imm) _mm_slli_si128(a, imm)

// Shift a right by imm8 bytes while shifting in zeros, and store the results in
// dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_bsrli_si128
// #define _mm_bsrli_si128(a, imm) _mm_srli_si128(a, imm)

// Cast vector of type __m128d to type __m128. This intrinsic is only used for
// compilation and does not generate any instructions, thus it has zero latency.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_castpd_ps
// FORCE_INLINE __m128 _mm_castpd_ps(__m128d a) {}

// Cast vector of type __m128d to type __m128i. This intrinsic is only used for
// compilation and does not generate any instructions, thus it has zero latency.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_castpd_si128
// FORCE_INLINE __m128i _mm_castpd_si128(__m128d a) {}

// Cast vector of type __m128 to type __m128d. This intrinsic is only used for
// compilation and does not generate any instructions, thus it has zero latency.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_castps_pd
// FORCE_INLINE __m128d _mm_castps_pd(__m128 a) {}

// Cast vector of type __m128 to type __m128i. This intrinsic is only used for
// compilation and does not generate any instructions, thus it has zero latency.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_castps_si128
// FORCE_INLINE __m128i _mm_castps_si128(__m128 a) {}

// Cast vector of type __m128i to type __m128d. This intrinsic is only used for
// compilation and does not generate any instructions, thus it has zero latency.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_castsi128_pd
// FORCE_INLINE __m128d _mm_castsi128_pd(__m128i a) {}

// Cast vector of type __m128i to type __m128. This intrinsic is only used for
// compilation and does not generate any instructions, thus it has zero latency.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_castsi128_ps
// FORCE_INLINE __m128 _mm_castsi128_ps(__m128i a) {}

// Invalidate and flush the cache line that contains p from all levels of the
// cache hierarchy.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_clflush
#if defined(__APPLE__)
#include <libkern/OSCacheControl.h>
#endif
// FORCE_INLINE void _mm_clflush(void const *p) {}

// Compare packed 16-bit integers in a and b for equality, and store the results
// in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpeq_epi16
// FORCE_INLINE __m128i _mm_cmpeq_epi16(__m128i a, __m128i b) {}

// Compare packed 32-bit integers in a and b for equality, and store the results
// in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpeq_epi32
// FORCE_INLINE __m128i _mm_cmpeq_epi32(__m128i a, __m128i b) {}

// Compare packed 8-bit integers in a and b for equality, and store the results
// in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpeq_epi8
// FORCE_INLINE __m128i _mm_cmpeq_epi8(__m128i a, __m128i b) {}

// Compare packed double-precision (64-bit) floating-point elements in a and b
// for equality, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpeq_pd
// FORCE_INLINE __m128d _mm_cmpeq_pd(__m128d a, __m128d b) {}

// Compare the lower double-precision (64-bit) floating-point elements in a and
// b for equality, store the result in the lower element of dst, and copy the
// upper element from a to the upper element of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpeq_sd
// FORCE_INLINE __m128d _mm_cmpeq_sd(__m128d a, __m128d b) {}

// Compare packed double-precision (64-bit) floating-point elements in a and b
// for greater-than-or-equal, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpge_pd
// FORCE_INLINE __m128d _mm_cmpge_pd(__m128d a, __m128d b) {}

// Compare the lower double-precision (64-bit) floating-point elements in a and
// b for greater-than-or-equal, store the result in the lower element of dst,
// and copy the upper element from a to the upper element of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpge_sd
// FORCE_INLINE __m128d _mm_cmpge_sd(__m128d a, __m128d b) {}

// Compare packed signed 16-bit integers in a and b for greater-than, and store
// the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpgt_epi16
// FORCE_INLINE __m128i _mm_cmpgt_epi16(__m128i a, __m128i b) {}

// Compare packed signed 32-bit integers in a and b for greater-than, and store
// the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpgt_epi32
// FORCE_INLINE __m128i _mm_cmpgt_epi32(__m128i a, __m128i b) {}

// Compare packed signed 8-bit integers in a and b for greater-than, and store
// the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpgt_epi8
// FORCE_INLINE __m128i _mm_cmpgt_epi8(__m128i a, __m128i b) {}

// Compare packed double-precision (64-bit) floating-point elements in a and b
// for greater-than, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpgt_pd
// FORCE_INLINE __m128d _mm_cmpgt_pd(__m128d a, __m128d b) {}

// Compare the lower double-precision (64-bit) floating-point elements in a and
// b for greater-than, store the result in the lower element of dst, and copy
// the upper element from a to the upper element of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpgt_sd
// FORCE_INLINE __m128d _mm_cmpgt_sd(__m128d a, __m128d b) {}

// Compare packed double-precision (64-bit) floating-point elements in a and b
// for less-than-or-equal, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmple_pd
// FORCE_INLINE __m128d _mm_cmple_pd(__m128d a, __m128d b) {}

// Compare the lower double-precision (64-bit) floating-point elements in a and
// b for less-than-or-equal, store the result in the lower element of dst, and
// copy the upper element from a to the upper element of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmple_sd
// FORCE_INLINE __m128d _mm_cmple_sd(__m128d a, __m128d b) {}

// Compare packed signed 16-bit integers in a and b for less-than, and store the
// results in dst. Note: This intrinsic emits the pcmpgtw instruction with the
// order of the operands switched.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmplt_epi16
// FORCE_INLINE __m128i _mm_cmplt_epi16(__m128i a, __m128i b) {}

// Compare packed signed 32-bit integers in a and b for less-than, and store the
// results in dst. Note: This intrinsic emits the pcmpgtd instruction with the
// order of the operands switched.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmplt_epi32
// FORCE_INLINE __m128i _mm_cmplt_epi32(__m128i a, __m128i b) {}

// Compare packed signed 8-bit integers in a and b for less-than, and store the
// results in dst. Note: This intrinsic emits the pcmpgtb instruction with the
// order of the operands switched.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmplt_epi8
// FORCE_INLINE __m128i _mm_cmplt_epi8(__m128i a, __m128i b) {}

// Compare packed double-precision (64-bit) floating-point elements in a and b
// for less-than, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmplt_pd
// FORCE_INLINE __m128d _mm_cmplt_pd(__m128d a, __m128d b) {}

// Compare the lower double-precision (64-bit) floating-point elements in a and
// b for less-than, store the result in the lower element of dst, and copy the
// upper element from a to the upper element of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmplt_sd
// FORCE_INLINE __m128d _mm_cmplt_sd(__m128d a, __m128d b) {}

// Compare packed double-precision (64-bit) floating-point elements in a and b
// for not-equal, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpneq_pd
// FORCE_INLINE __m128d _mm_cmpneq_pd(__m128d a, __m128d b) {}

// Compare the lower double-precision (64-bit) floating-point elements in a and
// b for not-equal, store the result in the lower element of dst, and copy the
// upper element from a to the upper element of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpneq_sd
// FORCE_INLINE __m128d _mm_cmpneq_sd(__m128d a, __m128d b) {}

// Compare packed double-precision (64-bit) floating-point elements in a and b
// for not-greater-than-or-equal, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpnge_pd
// FORCE_INLINE __m128d _mm_cmpnge_pd(__m128d a, __m128d b) {}

// Compare the lower double-precision (64-bit) floating-point elements in a and
// b for not-greater-than-or-equal, store the result in the lower element of
// dst, and copy the upper element from a to the upper element of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpnge_sd
// FORCE_INLINE __m128d _mm_cmpnge_sd(__m128d a, __m128d b) {}

// Compare packed double-precision (64-bit) floating-point elements in a and b
// for not-greater-than, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_cmpngt_pd
// FORCE_INLINE __m128d _mm_cmpngt_pd(__m128d a, __m128d b) {}

// Compare the lower double-precision (64-bit) floating-point elements in a and
// b for not-greater-than, store the result in the lower element of dst, and
// copy the upper element from a to the upper element of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpngt_sd
// FORCE_INLINE __m128d _mm_cmpngt_sd(__m128d a, __m128d b) {}

// Compare packed double-precision (64-bit) floating-point elements in a and b
// for not-less-than-or-equal, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpnle_pd
// FORCE_INLINE __m128d _mm_cmpnle_pd(__m128d a, __m128d b) {}

// Compare the lower double-precision (64-bit) floating-point elements in a and
// b for not-less-than-or-equal, store the result in the lower element of dst,
// and copy the upper element from a to the upper element of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpnle_sd
// FORCE_INLINE __m128d _mm_cmpnle_sd(__m128d a, __m128d b) {}

// Compare packed double-precision (64-bit) floating-point elements in a and b
// for not-less-than, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpnlt_pd
// FORCE_INLINE __m128d _mm_cmpnlt_pd(__m128d a, __m128d b) {}

// Compare the lower double-precision (64-bit) floating-point elements in a and
// b for not-less-than, store the result in the lower element of dst, and copy
// the upper element from a to the upper element of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpnlt_sd
// FORCE_INLINE __m128d _mm_cmpnlt_sd(__m128d a, __m128d b) {}

// Compare packed double-precision (64-bit) floating-point elements in a and b
// to see if neither is NaN, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpord_pd
// FORCE_INLINE __m128d _mm_cmpord_pd(__m128d a, __m128d b) {}

// Compare the lower double-precision (64-bit) floating-point elements in a and
// b to see if neither is NaN, store the result in the lower element of dst, and
// copy the upper element from a to the upper element of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpord_sd
// FORCE_INLINE __m128d _mm_cmpord_sd(__m128d a, __m128d b) {}

// Compare packed double-precision (64-bit) floating-point elements in a and b
// to see if either is NaN, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpunord_pd
// FORCE_INLINE __m128d _mm_cmpunord_pd(__m128d a, __m128d b) {}

// Compare the lower double-precision (64-bit) floating-point elements in a and
// b to see if either is NaN, store the result in the lower element of dst, and
// copy the upper element from a to the upper element of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpunord_sd
// FORCE_INLINE __m128d _mm_cmpunord_sd(__m128d a, __m128d b) {}

// Compare the lower double-precision (64-bit) floating-point element in a and b
// for greater-than-or-equal, and return the boolean result (0 or 1).
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_comige_sd
// FORCE_INLINE int _mm_comige_sd(__m128d a, __m128d b) {}

// Compare the lower double-precision (64-bit) floating-point element in a and b
// for greater-than, and return the boolean result (0 or 1).
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_comigt_sd
// FORCE_INLINE int _mm_comigt_sd(__m128d a, __m128d b) {}

// Compare the lower double-precision (64-bit) floating-point element in a and b
// for less-than-or-equal, and return the boolean result (0 or 1).
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_comile_sd
// FORCE_INLINE int _mm_comile_sd(__m128d a, __m128d b) {}

// Compare the lower double-precision (64-bit) floating-point element in a and b
// for less-than, and return the boolean result (0 or 1).
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_comilt_sd
// FORCE_INLINE int _mm_comilt_sd(__m128d a, __m128d b) {}

// Compare the lower double-precision (64-bit) floating-point element in a and b
// for equality, and return the boolean result (0 or 1).
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_comieq_sd
// FORCE_INLINE int _mm_comieq_sd(__m128d a, __m128d b) {}

// Compare the lower double-precision (64-bit) floating-point element in a and b
// for not-equal, and return the boolean result (0 or 1).
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_comineq_sd
// FORCE_INLINE int _mm_comineq_sd(__m128d a, __m128d b) {}

// Convert packed signed 32-bit integers in a to packed double-precision
// (64-bit) floating-point elements, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtepi32_pd
// FORCE_INLINE __m128d _mm_cvtepi32_pd(__m128i a) {}

// Convert packed signed 32-bit integers in a to packed single-precision
// (32-bit) floating-point elements, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtepi32_ps
// FORCE_INLINE __m128 _mm_cvtepi32_ps(__m128i a) {}

// Convert packed double-precision (64-bit) floating-point elements in a to
// packed 32-bit integers, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtpd_epi32
// FORCE_INLINE __m128i _mm_cvtpd_epi32(__m128d a) {}

// Convert packed double-precision (64-bit) floating-point elements in a to
// packed 32-bit integers, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtpd_pi32
// FORCE_INLINE __m64 _mm_cvtpd_pi32(__m128d a) {}

// Convert packed double-precision (64-bit) floating-point elements in a to
// packed single-precision (32-bit) floating-point elements, and store the
// results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtpd_ps
// FORCE_INLINE __m128 _mm_cvtpd_ps(__m128d a) {}

// Convert packed signed 32-bit integers in a to packed double-precision
// (64-bit) floating-point elements, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtpi32_pd
// FORCE_INLINE __m128d _mm_cvtpi32_pd(__m64 a) {}

// Convert packed single-precision (32-bit) floating-point elements in a to
// packed 32-bit integers, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtps_epi32
// FORCE_INLINE __m128i _mm_cvtps_epi32(__m128 a) {}

// Convert packed single-precision (32-bit) floating-point elements in a to
// packed double-precision (64-bit) floating-point elements, and store the
// results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtps_pd
// FORCE_INLINE __m128d _mm_cvtps_pd(__m128 a) {}

// Copy the lower double-precision (64-bit) floating-point element of a to dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtsd_f64
// FORCE_INLINE double _mm_cvtsd_f64(__m128d a) {}

// Convert the lower double-precision (64-bit) floating-point element in a to a
// 32-bit integer, and store the result in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtsd_si32
// FORCE_INLINE int32_t _mm_cvtsd_si32(__m128d a) {}

// Convert the lower double-precision (64-bit) floating-point element in a to a
// 64-bit integer, and store the result in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtsd_si64
// FORCE_INLINE int64_t _mm_cvtsd_si64(__m128d a) {}

// Convert the lower double-precision (64-bit) floating-point element in a to a
// 64-bit integer, and store the result in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtsd_si64x
// #define _mm_cvtsd_si64x _mm_cvtsd_si64

// Convert the lower double-precision (64-bit) floating-point element in b to a
// single-precision (32-bit) floating-point element, store the result in the
// lower element of dst, and copy the upper 3 packed elements from a to the
// upper elements of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtsd_ss
// FORCE_INLINE __m128 _mm_cvtsd_ss(__m128 a, __m128d b) {}

// Copy the lower 32-bit integer in a to dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtsi128_si32
// FORCE_INLINE int _mm_cvtsi128_si32(__m128i a) {}

// Copy the lower 64-bit integer in a to dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtsi128_si64
// FORCE_INLINE int64_t _mm_cvtsi128_si64(__m128i a) {}

// Copy the lower 64-bit integer in a to dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtsi128_si64x
// #define _mm_cvtsi128_si64x(a) _mm_cvtsi128_si64(a)

// Convert the signed 32-bit integer b to a double-precision (64-bit)
// floating-point element, store the result in the lower element of dst, and
// copy the upper element from a to the upper element of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtsi32_sd
// FORCE_INLINE __m128d _mm_cvtsi32_sd(__m128d a, int32_t b) {}

// Copy the lower 64-bit integer in a to dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtsi128_si64x
// #define _mm_cvtsi128_si64x(a) _mm_cvtsi128_si64(a)

// Copy 32-bit integer a to the lower elements of dst, and zero the upper
// elements of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtsi32_si128
// FORCE_INLINE __m128i _mm_cvtsi32_si128(int a) {}

// Convert the signed 64-bit integer b to a double-precision (64-bit)
// floating-point element, store the result in the lower element of dst, and
// copy the upper element from a to the upper element of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtsi64_sd
// FORCE_INLINE __m128d _mm_cvtsi64_sd(__m128d a, int64_t b) {}

// Copy 64-bit integer a to the lower element of dst, and zero the upper
// element.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtsi64_si128
// FORCE_INLINE __m128i _mm_cvtsi64_si128(int64_t a) {}

// Copy 64-bit integer a to the lower element of dst, and zero the upper
// element.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtsi64x_si128
// #define _mm_cvtsi64x_si128(a) _mm_cvtsi64_si128(a)

// Convert the signed 64-bit integer b to a double-precision (64-bit)
// floating-point element, store the result in the lower element of dst, and
// copy the upper element from a to the upper element of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtsi64x_sd
// #define _mm_cvtsi64x_sd(a, b) _mm_cvtsi64_sd(a, b)

// Convert the lower single-precision (32-bit) floating-point element in b to a
// double-precision (64-bit) floating-point element, store the result in the
// lower element of dst, and copy the upper element from a to the upper element
// of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtss_sd
// FORCE_INLINE __m128d _mm_cvtss_sd(__m128d a, __m128 b) {}

// Convert packed double-precision (64-bit) floating-point elements in a to
// packed 32-bit integers with truncation, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvttpd_epi32
// FORCE_INLINE __m128i _mm_cvttpd_epi32(__m128d a) {}

// Convert packed double-precision (64-bit) floating-point elements in a to
// packed 32-bit integers with truncation, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvttpd_pi32
// FORCE_INLINE __m64 _mm_cvttpd_pi32(__m128d a) {}

// Convert packed single-precision (32-bit) floating-point elements in a to
// packed 32-bit integers with truncation, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvttps_epi32
// FORCE_INLINE __m128i _mm_cvttps_epi32(__m128 a) {}

// Convert the lower double-precision (64-bit) floating-point element in a to a
// 32-bit integer with truncation, and store the result in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvttsd_si32
// FORCE_INLINE int32_t _mm_cvttsd_si32(__m128d a) {}

// Convert the lower double-precision (64-bit) floating-point element in a to a
// 64-bit integer with truncation, and store the result in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvttsd_si64
// FORCE_INLINE int64_t _mm_cvttsd_si64(__m128d a) {}

// Convert the lower double-precision (64-bit) floating-point element in a to a
// 64-bit integer with truncation, and store the result in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvttsd_si64x
// #define _mm_cvttsd_si64x(a) _mm_cvttsd_si64(a)

// Divide packed double-precision (64-bit) floating-point elements in a by
// packed elements in b, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_div_pd
// FORCE_INLINE __m128d _mm_div_pd(__m128d a, __m128d b) {}

// Divide the lower double-precision (64-bit) floating-point element in a by the
// lower double-precision (64-bit) floating-point element in b, store the result
// in the lower element of dst, and copy the upper element from a to the upper
// element of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_div_sd
// FORCE_INLINE __m128d _mm_div_sd(__m128d a, __m128d b) {}

// Extract a 16-bit integer from a, selected with imm8, and store the result in
// the lower element of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_extract_epi16
// FORCE_INLINE int _mm_extract_epi16(__m128i a, __constrange(0,8) int imm)
// #define _mm_extract_epi16(a, imm)

// Copy a to dst, and insert the 16-bit integer i into dst at the location
// specified by imm8.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_insert_epi16
// FORCE_INLINE __m128i _mm_insert_epi16(__m128i a, int b,
//                                       __constrange(0,8) int imm)
// #define _mm_insert_epi16(a, b, imm)

// Load 128-bits (composed of 2 packed double-precision (64-bit) floating-point
// elements) from memory into dst. mem_addr must be aligned on a 16-byte
// boundary or a general-protection exception may be generated.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_load_pd
// FORCE_INLINE __m128d _mm_load_pd(const double *p) {}

// Load a double-precision (64-bit) floating-point element from memory into both
// elements of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_load_pd1
// #define _mm_load_pd1 _mm_load1_pd

// Load a double-precision (64-bit) floating-point element from memory into the
// lower of dst, and zero the upper element. mem_addr does not need to be
// aligned on any particular boundary.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_load_sd
// FORCE_INLINE __m128d _mm_load_sd(const double *p) {}

// Load 128-bits of integer data from memory into dst. mem_addr must be aligned
// on a 16-byte boundary or a general-protection exception may be generated.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_load_si128
// FORCE_INLINE __m128i _mm_load_si128(const __m128i *p) {}

// Load a double-precision (64-bit) floating-point element from memory into both
// elements of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_load1_pd
// FORCE_INLINE __m128d _mm_load1_pd(const double *p) {}

// Load a double-precision (64-bit) floating-point element from memory into the
// upper element of dst, and copy the lower element from a to dst. mem_addr does
// not need to be aligned on any particular boundary.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_loadh_pd
// FORCE_INLINE __m128d _mm_loadh_pd(__m128d a, const double *p) {}

// Load 64-bit integer from memory into the first element of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_loadl_epi64
// FORCE_INLINE __m128i _mm_loadl_epi64(__m128i const *p) {}

// Load a double-precision (64-bit) floating-point element from memory into the
// lower element of dst, and copy the upper element from a to dst. mem_addr does
// not need to be aligned on any particular boundary.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_loadl_pd
// FORCE_INLINE __m128d _mm_loadl_pd(__m128d a, const double *p) {}

// Load 2 double-precision (64-bit) floating-point elements from memory into dst
// in reverse order. mem_addr must be aligned on a 16-byte boundary or a
// general-protection exception may be generated.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_loadr_pd
// FORCE_INLINE __m128d _mm_loadr_pd(const double *p) {}

// Loads two double-precision from unaligned memory, floating-point values.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_loadu_pd
// FORCE_INLINE __m128d _mm_loadu_pd(const double *p) {}

// Load 128-bits of integer data from memory into dst. mem_addr does not need to
// be aligned on any particular boundary.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_loadu_si128
// FORCE_INLINE __m128i _mm_loadu_si128(const __m128i *p) {}

// Load unaligned 32-bit integer from memory into the first element of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_loadu_si32
// FORCE_INLINE __m128i _mm_loadu_si32(const void *p) {}

// Multiply packed signed 16-bit integers in a and b, producing intermediate
// signed 32-bit integers. Horizontally add adjacent pairs of intermediate
// 32-bit integers, and pack the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_madd_epi16
// FORCE_INLINE __m128i _mm_madd_epi16(__m128i a, __m128i b) {}

// Conditionally store 8-bit integer elements from a into memory using mask
// (elements are not stored when the highest bit is not set in the corresponding
// element) and a non-temporal memory hint. mem_addr does not need to be aligned
// on any particular boundary.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_maskmoveu_si128
// FORCE_INLINE void _mm_maskmoveu_si128(__m128i a, __m128i mask, char *mem_addr) {}

// Compare packed signed 16-bit integers in a and b, and store packed maximum
// values in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_max_epi16
// FORCE_INLINE __m128i _mm_max_epi16(__m128i a, __m128i b) {}

// Compare packed unsigned 8-bit integers in a and b, and store packed maximum
// values in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_max_epu8
// FORCE_INLINE __m128i _mm_max_epu8(__m128i a, __m128i b) {}

// Compare packed double-precision (64-bit) floating-point elements in a and b,
// and store packed maximum values in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_max_pd
// FORCE_INLINE __m128d _mm_max_pd(__m128d a, __m128d b) {}

// Compare the lower double-precision (64-bit) floating-point elements in a and
// b, store the maximum value in the lower element of dst, and copy the upper
// element from a to the upper element of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_max_sd
// FORCE_INLINE __m128d _mm_max_sd(__m128d a, __m128d b) {}

// Compare packed signed 16-bit integers in a and b, and store packed minimum
// values in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_min_epi16
// FORCE_INLINE __m128i _mm_min_epi16(__m128i a, __m128i b) {}

// Compare packed unsigned 8-bit integers in a and b, and store packed minimum
// values in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_min_epu8
// FORCE_INLINE __m128i _mm_min_epu8(__m128i a, __m128i b) {}

// Compare packed double-precision (64-bit) floating-point elements in a and b,
// and store packed minimum values in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_min_pd
// FORCE_INLINE __m128d _mm_min_pd(__m128d a, __m128d b) {}

// Compare the lower double-precision (64-bit) floating-point elements in a and
// b, store the minimum value in the lower element of dst, and copy the upper
// element from a to the upper element of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_min_sd
// FORCE_INLINE __m128d _mm_min_sd(__m128d a, __m128d b) {}

// Copy the lower 64-bit integer in a to the lower element of dst, and zero the
// upper element.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_move_epi64
// FORCE_INLINE __m128i _mm_move_epi64(__m128i a) {}

// Move the lower double-precision (64-bit) floating-point element from b to the
// lower element of dst, and copy the upper element from a to the upper element
// of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_move_sd
// FORCE_INLINE __m128d _mm_move_sd(__m128d a, __m128d b) {}

// Create mask from the most significant bit of each 8-bit element in a, and
// store the result in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_movemask_epi8
// FORCE_INLINE int _mm_movemask_epi8(__m128i a) {}

// Set each bit of mask dst based on the most significant bit of the
// corresponding packed double-precision (64-bit) floating-point element in a.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_movemask_pd
// FORCE_INLINE int _mm_movemask_pd(__m128d a) {}

// Copy the lower 64-bit integer in a to dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_movepi64_pi64
// FORCE_INLINE __m64 _mm_movepi64_pi64(__m128i a) {}

// Copy the 64-bit integer a to the lower element of dst, and zero the upper
// element.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_movpi64_epi64
// FORCE_INLINE __m128i _mm_movpi64_epi64(__m64 a) {}

// Multiply the low unsigned 32-bit integers from each packed 64-bit element in
// a and b, and store the unsigned 64-bit results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_mul_epu32
// FORCE_INLINE __m128i _mm_mul_epu32(__m128i a, __m128i b) {}

// Multiply packed double-precision (64-bit) floating-point elements in a and b,
// and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_mul_pd
// FORCE_INLINE __m128d _mm_mul_pd(__m128d a, __m128d b) {}

// Multiply the lower double-precision (64-bit) floating-point element in a and
// b, store the result in the lower element of dst, and copy the upper element
// from a to the upper element of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=mm_mul_sd
// FORCE_INLINE __m128d _mm_mul_sd(__m128d a, __m128d b) {}

// Multiply the low unsigned 32-bit integers from a and b, and store the
// unsigned 64-bit result in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_mul_su32
// FORCE_INLINE __m64 _mm_mul_su32(__m64 a, __m64 b) {}

// Multiply the packed signed 16-bit integers in a and b, producing intermediate
// 32-bit integers, and store the high 16 bits of the intermediate integers in
// dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_mulhi_epi16
// FORCE_INLINE __m128i _mm_mulhi_epi16(__m128i a, __m128i b) {}

// Multiply the packed unsigned 16-bit integers in a and b, producing
// intermediate 32-bit integers, and store the high 16 bits of the intermediate
// integers in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_mulhi_epu16
// FORCE_INLINE __m128i _mm_mulhi_epu16(__m128i a, __m128i b) {}

// Multiply the packed 16-bit integers in a and b, producing intermediate 32-bit
// integers, and store the low 16 bits of the intermediate integers in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_mullo_epi16
// FORCE_INLINE __m128i _mm_mullo_epi16(__m128i a, __m128i b) {}

// Compute the bitwise OR of packed double-precision (64-bit) floating-point
// elements in a and b, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=mm_or_pd
// FORCE_INLINE __m128d _mm_or_pd(__m128d a, __m128d b) {}

// Compute the bitwise OR of 128 bits (representing integer data) in a and b,
// and store the result in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_or_si128
// FORCE_INLINE __m128i _mm_or_si128(__m128i a, __m128i b) {}

// Convert packed signed 16-bit integers from a and b to packed 8-bit integers
// using signed saturation, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_packs_epi16
// FORCE_INLINE __m128i _mm_packs_epi16(__m128i a, __m128i b) {}

// Convert packed signed 32-bit integers from a and b to packed 16-bit integers
// using signed saturation, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_packs_epi32
// FORCE_INLINE __m128i _mm_packs_epi32(__m128i a, __m128i b) {}

// Convert packed signed 16-bit integers from a and b to packed 8-bit integers
// using unsigned saturation, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_packus_epi16
// FORCE_INLINE __m128i _mm_packus_epi16(const __m128i a, const __m128i b) {}

// Pause the processor. This is typically used in spin-wait loops and depending
// on the x86 processor typical values are in the 40-100 cycle range. The
// 'yield' instruction isn't a good fit because it's effectively a nop on most
// Arm cores. Experience with several databases has shown has shown an 'isb' is
// a reasonable approximation.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_pause
// FORCE_INLINE void _mm_pause(void) {}

// Compute the absolute differences of packed unsigned 8-bit integers in a and
// b, then horizontally sum each consecutive 8 differences to produce two
// unsigned 16-bit integers, and pack these unsigned 16-bit integers in the low
// 16 bits of 64-bit elements in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_sad_epu8
// FORCE_INLINE __m128i _mm_sad_epu8(__m128i a, __m128i b) {}

// Set packed 16-bit integers in dst with the supplied values.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_set_epi16
// FORCE_INLINE __m128i _mm_set_epi16(short i7,
//                                    short i6,
//                                    short i5,
//                                    short i4,
//                                    short i3,
//                                    short i2,
//                                    short i1,
//                                    short i0) {}

// Set packed 32-bit integers in dst with the supplied values.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_set_epi32
// FORCE_INLINE __m128i _mm_set_epi32(int i3, int i2, int i1, int i0) {}

// Set packed 64-bit integers in dst with the supplied values.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_set_epi64
// FORCE_INLINE __m128i _mm_set_epi64(__m64 i1, __m64 i2) {}

// Set packed 64-bit integers in dst with the supplied values.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_set_epi64x
// FORCE_INLINE __m128i _mm_set_epi64x(int64_t i1, int64_t i2) {}

// Set packed 8-bit integers in dst with the supplied values.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_set_epi8
// FORCE_INLINE __m128i _mm_set_epi8(signed char b15,
//                                   signed char b14,
//                                   signed char b13,
//                                   signed char b12,
//                                   signed char b11,
//                                   signed char b10,
//                                   signed char b9,
//                                   signed char b8,
//                                   signed char b7,
//                                   signed char b6,
//                                   signed char b5,
//                                   signed char b4,
//                                   signed char b3,
//                                   signed char b2,
//                                   signed char b1,
//                                   signed char b0) {}

// Set packed double-precision (64-bit) floating-point elements in dst with the
// supplied values.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_set_pd
// FORCE_INLINE __m128d _mm_set_pd(double e1, double e0) {}

// Broadcast double-precision (64-bit) floating-point value a to all elements of
// dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_set_pd1
// #define _mm_set_pd1 _mm_set1_pd

// Copy double-precision (64-bit) floating-point element a to the lower element
// of dst, and zero the upper element.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_set_sd
// FORCE_INLINE __m128d _mm_set_sd(double a) {}

// Broadcast 16-bit integer a to all elements of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_set1_epi16
// FORCE_INLINE __m128i _mm_set1_epi16(short w) {}

// Broadcast 32-bit integer a to all elements of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_set1_epi32
// FORCE_INLINE __m128i _mm_set1_epi32(int _i) {}

// Broadcast 64-bit integer a to all elements of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_set1_epi64
// FORCE_INLINE __m128i _mm_set1_epi64(__m64 _i) {}

// Broadcast 64-bit integer a to all elements of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_set1_epi64x
// FORCE_INLINE __m128i _mm_set1_epi64x(int64_t _i) {}

// Broadcast 8-bit integer a to all elements of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_set1_epi8
// FORCE_INLINE __m128i _mm_set1_epi8(signed char w) {}

// Broadcast double-precision (64-bit) floating-point value a to all elements of
// dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_set1_pd
// FORCE_INLINE __m128d _mm_set1_pd(double d) {}

// Set packed 16-bit integers in dst with the supplied values in reverse order.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_setr_epi16
// FORCE_INLINE __m128i _mm_setr_epi16(short w0,
//                                     short w1,
//                                     short w2,
//                                     short w3,
//                                     short w4,
//                                     short w5,
//                                     short w6,
//                                     short w7) {}

// Set packed 32-bit integers in dst with the supplied values in reverse order.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_setr_epi32
// FORCE_INLINE __m128i _mm_setr_epi32(int i3, int i2, int i1, int i0) {}

// Set packed 64-bit integers in dst with the supplied values in reverse order.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_setr_epi64
// FORCE_INLINE __m128i _mm_setr_epi64(__m64 e1, __m64 e0) {}

// Set packed 8-bit integers in dst with the supplied values in reverse order.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_setr_epi8
// FORCE_INLINE __m128i _mm_setr_epi8(signed char b0,
//                                    signed char b1,
//                                    signed char b2,
//                                    signed char b3,
//                                    signed char b4,
//                                    signed char b5,
//                                    signed char b6,
//                                    signed char b7,
//                                    signed char b8,
//                                    signed char b9,
//                                    signed char b10,
//                                    signed char b11,
//                                    signed char b12,
//                                    signed char b13,
//                                    signed char b14,
//                                    signed char b15) {}

// Set packed double-precision (64-bit) floating-point elements in dst with the
// supplied values in reverse order.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_setr_pd
// FORCE_INLINE __m128d _mm_setr_pd(double e1, double e0) {}

// Return vector of type __m128d with all elements set to zero.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_setzero_pd
// FORCE_INLINE __m128d _mm_setzero_pd(void) {}

// Return vector of type __m128i with all elements set to zero.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_setzero_si128
// FORCE_INLINE __m128i _mm_setzero_si128(void) {}

// Shuffle 32-bit integers in a using the control in imm8, and store the results
// in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_shuffle_epi32
// FORCE_INLINE __m128i _mm_shuffle_epi32(__m128i a,
//                                        __constrange(0,255) int imm)
#if defined(_sse2rvv_shuffle)
// #define _mm_shuffle_epi32(a, imm)
#else  // generic
// #define _mm_shuffle_epi32(a, imm)
#endif

// Shuffle double-precision (64-bit) floating-point elements using the control
// in imm8, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_shuffle_pd
#ifdef _sse2rvv_shuffle
// #define _mm_shuffle_pd(a, b, imm8)
#else
// #define _mm_shuffle_pd(a, b, imm8)
#endif

// FORCE_INLINE __m128i _mm_shufflehi_epi16(__m128i a,
//                                          __constrange(0,255) int imm)
#if defined(_sse2rvv_shuffle)
// #define _mm_shufflehi_epi16(a, imm)
#else  // generic
// #define _mm_shufflehi_epi16(a, imm) _mm_shufflehi_epi16_function((a), (imm))
#endif

// FORCE_INLINE __m128i _mm_shufflelo_epi16(__m128i a,
//                                          __constrange(0,255) int imm)
#if defined(_sse2rvv_shuffle)
// #define _mm_shufflelo_epi16(a, imm)
#else  // generic
// #define _mm_shufflelo_epi16(a, imm) _mm_shufflelo_epi16_function((a), (imm))
#endif

// Shift packed 16-bit integers in a left by count while shifting in zeros, and
// store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_sll_epi16
// FORCE_INLINE __m128i _mm_sll_epi16(__m128i a, __m128i count) {}

// Shift packed 32-bit integers in a left by count while shifting in zeros, and
// store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_sll_epi32
// FORCE_INLINE __m128i _mm_sll_epi32(__m128i a, __m128i count) {}

// Shift packed 64-bit integers in a left by count while shifting in zeros, and
// store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_sll_epi64
// FORCE_INLINE __m128i _mm_sll_epi64(__m128i a, __m128i count) {}

// Shift packed 16-bit integers in a left by imm8 while shifting in zeros, and
// store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_slli_epi16
// FORCE_INLINE __m128i _mm_slli_epi16(__m128i a, int imm) {}

// Shift packed 32-bit integers in a left by imm8 while shifting in zeros, and
// store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_slli_epi32
// FORCE_INLINE __m128i _mm_slli_epi32(__m128i a, int imm) {}

// Shift packed 64-bit integers in a left by imm8 while shifting in zeros, and
// store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_slli_epi64
// FORCE_INLINE __m128i _mm_slli_epi64(__m128i a, int imm) {}

// Shift a left by imm8 bytes while shifting in zeros, and store the results in
// dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_slli_si128
// #define _mm_slli_si128(a, imm)

// Compute the square root of packed double-precision (64-bit) floating-point
// elements in a, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_sqrt_pd
// FORCE_INLINE __m128d _mm_sqrt_pd(__m128d a) {}

// Compute the square root of the lower double-precision (64-bit) floating-point
// element in b, store the result in the lower element of dst, and copy the
// upper element from a to the upper element of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_sqrt_sd
// FORCE_INLINE __m128d _mm_sqrt_sd(__m128d a, __m128d b) {}

// Shift packed 16-bit integers in a right by count while shifting in sign bits,
// and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_sra_epi16
// FORCE_INLINE __m128i _mm_sra_epi16(__m128i a, __m128i count) {}

// Shift packed 32-bit integers in a right by count while shifting in sign bits,
// and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_sra_epi32
// FORCE_INLINE __m128i _mm_sra_epi32(__m128i a, __m128i count) {}

// Shift packed 16-bit integers in a right by imm8 while shifting in sign
// bits, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_srai_epi16
// FORCE_INLINE __m128i _mm_srai_epi16(__m128i a, int imm) {}

// Shift packed 32-bit integers in a right by imm8 while shifting in sign bits,
// and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_srai_epi32
// FORCE_INLINE __m128i _mm_srai_epi32(__m128i a, __constrange(0,255) int imm)
// #define _mm_srai_epi32(a, imm)

// Shift packed 16-bit integers in a right by count while shifting in zeros, and
// store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_srl_epi16
// FORCE_INLINE __m128i _mm_srl_epi16(__m128i a, __m128i count) {}

// Shift packed 32-bit integers in a right by count while shifting in zeros, and
// store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_srl_epi32
// FORCE_INLINE __m128i _mm_srl_epi32(__m128i a, __m128i count) {}

// Shift packed 64-bit integers in a right by count while shifting in zeros, and
// store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_srl_epi64
// FORCE_INLINE __m128i _mm_srl_epi64(__m128i a, __m128i count) {}

// Shift packed 16-bit integers in a right by imm8 while shifting in zeros, and
// store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_srli_epi16
// #define _mm_srli_epi16(a, imm)

// Shift packed 32-bit integers in a right by imm8 while shifting in zeros, and
// store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_srli_epi32
// FORCE_INLINE __m128i _mm_srli_epi32(__m128i a, __constrange(0,255) int imm)
// #define _mm_srli_epi32(a, imm)

// Shift packed 64-bit integers in a right by imm8 while shifting in zeros, and
// store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_srli_epi64
// #define _mm_srli_epi64(a, imm)

// Shift a right by imm8 bytes while shifting in zeros, and store the results in
// dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_srli_si128
// #define _mm_srli_si128(a, imm)

// Store 128-bits (composed of 2 packed double-precision (64-bit) floating-point
// elements) from a into memory. mem_addr must be aligned on a 16-byte boundary
// or a general-protection exception may be generated.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_store_pd
// FORCE_INLINE void _mm_store_pd(double *mem_addr, __m128d a) {}

// Store the lower double-precision (64-bit) floating-point element from a into
// 2 contiguous elements in memory. mem_addr must be aligned on a 16-byte
// boundary or a general-protection exception may be generated.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_store_pd1
// FORCE_INLINE void _mm_store_pd1(double *mem_addr, __m128d a) {}

// Store the lower double-precision (64-bit) floating-point element from a into
// memory. mem_addr does not need to be aligned on any particular boundary.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=mm_store_sd
// FORCE_INLINE void _mm_store_sd(double *mem_addr, __m128d a) {}

// Store 128-bits of integer data from a into memory. mem_addr must be aligned
// on a 16-byte boundary or a general-protection exception may be generated.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_store_si128
// FORCE_INLINE void _mm_store_si128(__m128i *p, __m128i a) {}

// Store the lower double-precision (64-bit) floating-point element from a into
// 2 contiguous elements in memory. mem_addr must be aligned on a 16-byte
// boundary or a general-protection exception may be generated.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#expand=9,526,5601&text=_mm_store1_pd
// #define _mm_store1_pd _mm_store_pd1

// Store the upper double-precision (64-bit) floating-point element from a into
// memory.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_storeh_pd
// FORCE_INLINE void _mm_storeh_pd(double *mem_addr, __m128d a) {}

// Store 64-bit integer from the first element of a into memory.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_storel_epi64
// FORCE_INLINE void _mm_storel_epi64(__m128i *a, __m128i b) {}

// Store the lower double-precision (64-bit) floating-point element from a into
// memory.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_storel_pd
// FORCE_INLINE void _mm_storel_pd(double *mem_addr, __m128d a) {}

// Store 2 double-precision (64-bit) floating-point elements from a into memory
// in reverse order. mem_addr must be aligned on a 16-byte boundary or a
// general-protection exception may be generated.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_storer_pd
// FORCE_INLINE void _mm_storer_pd(double *mem_addr, __m128d a) {}

// Store 128-bits (composed of 2 packed double-precision (64-bit) floating-point
// elements) from a into memory. mem_addr does not need to be aligned on any
// particular boundary.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_storeu_pd
// FORCE_INLINE void _mm_storeu_pd(double *mem_addr, __m128d a) {}

// Store 128-bits of integer data from a into memory. mem_addr does not need to
// be aligned on any particular boundary.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_storeu_si128
// FORCE_INLINE void _mm_storeu_si128(__m128i *p, __m128i a) {}

// Store 32-bit integer from the first element of a into memory. mem_addr does
// not need to be aligned on any particular boundary.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_storeu_si32
// FORCE_INLINE void _mm_storeu_si32(void *p, __m128i a) {}

// Store 128-bits (composed of 2 packed double-precision (64-bit) floating-point
// elements) from a into memory using a non-temporal memory hint. mem_addr must
// be aligned on a 16-byte boundary or a general-protection exception may be
// generated.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_stream_pd
// FORCE_INLINE void _mm_stream_pd(double *p, __m128d a) {}

// Store 128-bits of integer data from a into memory using a non-temporal memory
// hint. mem_addr must be aligned on a 16-byte boundary or a general-protection
// exception may be generated.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_stream_si128
// FORCE_INLINE void _mm_stream_si128(__m128i *p, __m128i a) {}

// Store 32-bit integer a into memory using a non-temporal hint to minimize
// cache pollution. If the cache line containing address mem_addr is already in
// the cache, the cache will be updated.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_stream_si32
// FORCE_INLINE void _mm_stream_si32(int *p, int a) {}

// Store 64-bit integer a into memory using a non-temporal hint to minimize
// cache pollution. If the cache line containing address mem_addr is already in
// the cache, the cache will be updated.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_stream_si64
// FORCE_INLINE void _mm_stream_si64(__int64 *p, __int64 a) {}

// Subtract packed 16-bit integers in b from packed 16-bit integers in a, and
// store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_sub_epi16
// FORCE_INLINE __m128i _mm_sub_epi16(__m128i a, __m128i b) {}

// Subtract packed 32-bit integers in b from packed 32-bit integers in a, and
// store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_sub_epi32
// FORCE_INLINE __m128i _mm_sub_epi32(__m128i a, __m128i b) {}

// Subtract packed 64-bit integers in b from packed 64-bit integers in a, and
// store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_sub_epi64
// FORCE_INLINE __m128i _mm_sub_epi64(__m128i a, __m128i b) {}

// Subtract packed 8-bit integers in b from packed 8-bit integers in a, and
// store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_sub_epi8
// FORCE_INLINE __m128i _mm_sub_epi8(__m128i a, __m128i b) {}

// Subtract packed double-precision (64-bit) floating-point elements in b from
// packed double-precision (64-bit) floating-point elements in a, and store the
// results in dst.
//  https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=mm_sub_pd
// FORCE_INLINE __m128d _mm_sub_pd(__m128d a, __m128d b) {}

// Subtract the lower double-precision (64-bit) floating-point element in b from
// the lower double-precision (64-bit) floating-point element in a, store the
// result in the lower element of dst, and copy the upper element from a to the
// upper element of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_sub_sd
// FORCE_INLINE __m128d _mm_sub_sd(__m128d a, __m128d b) {}

// Subtract 64-bit integer b from 64-bit integer a, and store the result in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_sub_si64
// FORCE_INLINE __m64 _mm_sub_si64(__m64 a, __m64 b) {}

// Subtract packed signed 16-bit integers in b from packed 16-bit integers in a
// using saturation, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_subs_epi16
// FORCE_INLINE __m128i _mm_subs_epi16(__m128i a, __m128i b) {}

// Subtract packed signed 8-bit integers in b from packed 8-bit integers in a
// using saturation, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_subs_epi8
// FORCE_INLINE __m128i _mm_subs_epi8(__m128i a, __m128i b) {}

// Subtract packed unsigned 16-bit integers in b from packed unsigned 16-bit
// integers in a using saturation, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_subs_epu16
// FORCE_INLINE __m128i _mm_subs_epu16(__m128i a, __m128i b) {}

// Subtract packed unsigned 8-bit integers in b from packed unsigned 8-bit
// integers in a using saturation, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_subs_epu8
// FORCE_INLINE __m128i _mm_subs_epu8(__m128i a, __m128i b) {}

// #define _mm_ucomieq_sd _mm_comieq_sd
// #define _mm_ucomige_sd _mm_comige_sd
// #define _mm_ucomigt_sd _mm_comigt_sd
// #define _mm_ucomile_sd _mm_comile_sd
// #define _mm_ucomilt_sd _mm_comilt_sd
// #define _mm_ucomineq_sd _mm_comineq_sd

// Return vector of type __m128d with undefined elements.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_undefined_pd
// FORCE_INLINE __m128d _mm_undefined_pd(void) {}

// Unpack and interleave 16-bit integers from the high half of a and b, and
// store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_unpackhi_epi16
// FORCE_INLINE __m128i _mm_unpackhi_epi16(__m128i a, __m128i b) {}

// Unpack and interleave 32-bit integers from the high half of a and b, and
// store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_unpackhi_epi32
// FORCE_INLINE __m128i _mm_unpackhi_epi32(__m128i a, __m128i b) {}

// Unpack and interleave 64-bit integers from the high half of a and b, and
// store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_unpackhi_epi64
// FORCE_INLINE __m128i _mm_unpackhi_epi64(__m128i a, __m128i b) {}

// Unpack and interleave 8-bit integers from the high half of a and b, and store
// the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_unpackhi_epi8
// FORCE_INLINE __m128i _mm_unpackhi_epi8(__m128i a, __m128i b) {}

// Unpack and interleave double-precision (64-bit) floating-point elements from
// the high half of a and b, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_unpackhi_pd
// FORCE_INLINE __m128d _mm_unpackhi_pd(__m128d a, __m128d b) {}

// Unpack and interleave 16-bit integers from the low half of a and b, and store
// the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_unpacklo_epi16
// FORCE_INLINE __m128i _mm_unpacklo_epi16(__m128i a, __m128i b) {}

// Unpack and interleave 32-bit integers from the low half of a and b, and store
// the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_unpacklo_epi32
// FORCE_INLINE __m128i _mm_unpacklo_epi32(__m128i a, __m128i b) {}

// Unpack and interleave 64-bit integers from the low half of a and b, and store
// the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_unpacklo_epi64
// FORCE_INLINE __m128i _mm_unpacklo_epi64(__m128i a, __m128i b) {}

// Unpack and interleave 8-bit integers from the low half of a and b, and store
// the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_unpacklo_epi8
// FORCE_INLINE __m128i _mm_unpacklo_epi8(__m128i a, __m128i b) {}

// Unpack and interleave double-precision (64-bit) floating-point elements from
// the low half of a and b, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_unpacklo_pd
// FORCE_INLINE __m128d _mm_unpacklo_pd(__m128d a, __m128d b) {}

// Compute the bitwise XOR of packed double-precision (64-bit) floating-point
// elements in a and b, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_xor_pd
// FORCE_INLINE __m128d _mm_xor_pd(__m128d a, __m128d b) {}

// Compute the bitwise XOR of 128 bits (representing integer data) in a and b,
// and store the result in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_xor_si128
// FORCE_INLINE __m128i _mm_xor_si128(__m128i a, __m128i b) {}

/* SSE3 */

// Alternatively add and subtract packed double-precision (64-bit)
// floating-point elements in a to/from packed elements in b, and store the
// results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_addsub_pd
// FORCE_INLINE __m128d _mm_addsub_pd(__m128d a, __m128d b) {}

// Alternatively add and subtract packed single-precision (32-bit)
// floating-point elements in a to/from packed elements in b, and store the
// results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=addsub_ps
// FORCE_INLINE __m128 _mm_addsub_ps(__m128 a, __m128 b) {}

// Horizontally add adjacent pairs of double-precision (64-bit) floating-point
// elements in a and b, and pack the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_hadd_pd
// FORCE_INLINE __m128d _mm_hadd_pd(__m128d a, __m128d b) {}

// Horizontally add adjacent pairs of single-precision (32-bit) floating-point
// elements in a and b, and pack the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_hadd_ps
// FORCE_INLINE __m128 _mm_hadd_ps(__m128 a, __m128 b) {}

// Horizontally subtract adjacent pairs of double-precision (64-bit)
// floating-point elements in a and b, and pack the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_hsub_pd
// FORCE_INLINE __m128d _mm_hsub_pd(__m128d _a, __m128d _b) {}

// Horizontally subtract adjacent pairs of single-precision (32-bit)
// floating-point elements in a and b, and pack the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_hsub_ps
// FORCE_INLINE __m128 _mm_hsub_ps(__m128 _a, __m128 _b) {}

// Load 128-bits of integer data from unaligned memory into dst. This intrinsic
// may perform better than _mm_loadu_si128 when the data crosses a cache line
// boundary.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_lddqu_si128
// #define _mm_lddqu_si128 _mm_loadu_si128

// Load a double-precision (64-bit) floating-point element from memory into both
// elements of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_loaddup_pd
// #define _mm_loaddup_pd _mm_load1_pd

// Duplicate the low double-precision (64-bit) floating-point element from a,
// and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_movedup_pd
// FORCE_INLINE __m128d _mm_movedup_pd(__m128d a) {}

// Duplicate odd-indexed single-precision (32-bit) floating-point elements
// from a, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_movehdup_ps
// FORCE_INLINE __m128 _mm_movehdup_ps(__m128 a) {}

// Duplicate even-indexed single-precision (32-bit) floating-point elements
// from a, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_moveldup_ps
// FORCE_INLINE __m128 _mm_moveldup_ps(__m128 a) {}

/* SSSE3 */

// Compute the absolute value of packed signed 16-bit integers in a, and store
// the unsigned results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_abs_epi16
// FORCE_INLINE __m128i _mm_abs_epi16(__m128i a) {}

// Compute the absolute value of packed signed 32-bit integers in a, and store
// the unsigned results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_abs_epi32
// FORCE_INLINE __m128i _mm_abs_epi32(__m128i a) {}

// Compute the absolute value of packed signed 8-bit integers in a, and store
// the unsigned results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_abs_epi8
// FORCE_INLINE __m128i _mm_abs_epi8(__m128i a) {}

// Compute the absolute value of packed signed 16-bit integers in a, and store
// the unsigned results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_abs_pi16
// FORCE_INLINE __m64 _mm_abs_pi16(__m64 a) {}

// Compute the absolute value of packed signed 32-bit integers in a, and store
// the unsigned results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_abs_pi32
// FORCE_INLINE __m64 _mm_abs_pi32(__m64 a) {}

// Compute the absolute value of packed signed 8-bit integers in a, and store
// the unsigned results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_abs_pi8
// FORCE_INLINE __m64 _mm_abs_pi8(__m64 a) {}

// Concatenate 16-byte blocks in a and b into a 32-byte temporary result, shift
// the result right by imm8 bytes, and store the low 16 bytes in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_alignr_epi8
#if defined(__GNUC__) && !defined(__clang__)
// #define _mm_alignr_epi8(a, b, imm)

#else
// #define _mm_alignr_epi8(a, b, imm)

#endif

// Concatenate 8-byte blocks in a and b into a 16-byte temporary result, shift
// the result right by imm8 bytes, and store the low 8 bytes in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_alignr_pi8
// #define _mm_alignr_pi8(a, b, imm)

// Horizontally add adjacent pairs of 16-bit integers in a and b, and pack the
// signed 16-bit results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_hadd_epi16
// FORCE_INLINE __m128i _mm_hadd_epi16(__m128i _a, __m128i _b) {}

// Horizontally add adjacent pairs of 32-bit integers in a and b, and pack the
// signed 32-bit results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_hadd_epi32
// FORCE_INLINE __m128i _mm_hadd_epi32(__m128i _a, __m128i _b) {}

// Horizontally add adjacent pairs of 16-bit integers in a and b, and pack the
// signed 16-bit results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_hadd_pi16
// FORCE_INLINE __m64 _mm_hadd_pi16(__m64 a, __m64 b) {}

// Horizontally add adjacent pairs of 32-bit integers in a and b, and pack the
// signed 32-bit results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_hadd_pi32
// FORCE_INLINE __m64 _mm_hadd_pi32(__m64 a, __m64 b) {}

// Horizontally add adjacent pairs of signed 16-bit integers in a and b using
// saturation, and pack the signed 16-bit results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_hadds_epi16
// FORCE_INLINE __m128i _mm_hadds_epi16(__m128i _a, __m128i _b) {}

// Horizontally add adjacent pairs of signed 16-bit integers in a and b using
// saturation, and pack the signed 16-bit results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_hadds_pi16
// FORCE_INLINE __m64 _mm_hadds_pi16(__m64 _a, __m64 _b) {}

// Horizontally subtract adjacent pairs of 16-bit integers in a and b, and pack
// the signed 16-bit results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_hsub_epi16
// FORCE_INLINE __m128i _mm_hsub_epi16(__m128i _a, __m128i _b) {}

// Horizontally subtract adjacent pairs of 32-bit integers in a and b, and pack
// the signed 32-bit results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_hsub_epi32
// FORCE_INLINE __m128i _mm_hsub_epi32(__m128i _a, __m128i _b) {}

// Horizontally subtract adjacent pairs of 16-bit integers in a and b, and pack
// the signed 16-bit results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_hsub_pi16
// FORCE_INLINE __m64 _mm_hsub_pi16(__m64 _a, __m64 _b) {}

// Horizontally subtract adjacent pairs of 32-bit integers in a and b, and pack
// the signed 32-bit results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=mm_hsub_pi32
// FORCE_INLINE __m64 _mm_hsub_pi32(__m64 _a, __m64 _b) {}

// Horizontally subtract adjacent pairs of signed 16-bit integers in a and b
// using saturation, and pack the signed 16-bit results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_hsubs_epi16
// FORCE_INLINE __m128i _mm_hsubs_epi16(__m128i _a, __m128i _b) {}

// Horizontally subtract adjacent pairs of signed 16-bit integers in a and b
// using saturation, and pack the signed 16-bit results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_hsubs_pi16
// FORCE_INLINE __m64 _mm_hsubs_pi16(__m64 _a, __m64 _b) {}

// Vertically multiply each unsigned 8-bit integer from a with the corresponding
// signed 8-bit integer from b, producing intermediate signed 16-bit integers.
// Horizontally add adjacent pairs of intermediate signed 16-bit integers,
// and pack the saturated results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_maddubs_epi16
// FORCE_INLINE __m128i _mm_maddubs_epi16(__m128i _a, __m128i _b) {}

// Vertically multiply each unsigned 8-bit integer from a with the corresponding
// signed 8-bit integer from b, producing intermediate signed 16-bit integers.
// Horizontally add adjacent pairs of intermediate signed 16-bit integers, and
// pack the saturated results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_maddubs_pi16
// FORCE_INLINE __m64 _mm_maddubs_pi16(__m64 _a, __m64 _b) {}

// Multiply packed signed 16-bit integers in a and b, producing intermediate
// signed 32-bit integers. Shift right by 15 bits while rounding up, and store
// the packed 16-bit integers in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_mulhrs_epi16
// FORCE_INLINE __m128i _mm_mulhrs_epi16(__m128i a, __m128i b) {}

// Multiply packed signed 16-bit integers in a and b, producing intermediate
// signed 32-bit integers. Truncate each intermediate integer to the 18 most
// significant bits, round by adding 1, and store bits [16:1] to dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_mulhrs_pi16
// FORCE_INLINE __m64 _mm_mulhrs_pi16(__m64 a, __m64 b) {}

// Shuffle packed 8-bit integers in a according to shuffle control mask in the
// corresponding 8-bit element of b, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_shuffle_epi8
// FORCE_INLINE __m128i _mm_shuffle_epi8(__m128i a, __m128i b) {}

// Shuffle packed 8-bit integers in a according to shuffle control mask in the
// corresponding 8-bit element of b, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_shuffle_pi8
// FORCE_INLINE __m64 _mm_shuffle_pi8(__m64 a, __m64 b) {}

// Negate packed 16-bit integers in a when the corresponding signed
// 16-bit integer in b is negative, and store the results in dst.
// Element in dst are zeroed out when the corresponding element
// in b is zero.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_sign_epi16
// FORCE_INLINE __m128i _mm_sign_epi16(__m128i _a, __m128i _b) {}

// Negate packed 32-bit integers in a when the corresponding signed
// 32-bit integer in b is negative, and store the results in dst.
// Element in dst are zeroed out when the corresponding element
// in b is zero.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_sign_epi32
// FORCE_INLINE __m128i _mm_sign_epi32(__m128i _a, __m128i _b) {}

// Negate packed 8-bit integers in a when the corresponding signed
// 8-bit integer in b is negative, and store the results in dst.
// Element in dst are zeroed out when the corresponding element
// in b is zero.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_sign_epi8
// FORCE_INLINE __m128i _mm_sign_epi8(__m128i _a, __m128i _b) {}

// Negate packed 16-bit integers in a when the corresponding signed 16-bit
// integer in b is negative, and store the results in dst. Element in dst are
// zeroed out when the corresponding element in b is zero.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_sign_pi16
// FORCE_INLINE __m64 _mm_sign_pi16(__m64 _a, __m64 _b) {}

// Negate packed 32-bit integers in a when the corresponding signed 32-bit
// integer in b is negative, and store the results in dst. Element in dst are
// zeroed out when the corresponding element in b is zero.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_sign_pi32
// FORCE_INLINE __m64 _mm_sign_pi32(__m64 _a, __m64 _b) {}

// Negate packed 8-bit integers in a when the corresponding signed 8-bit integer
// in b is negative, and store the results in dst. Element in dst are zeroed out
// when the corresponding element in b is zero.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_sign_pi8
// FORCE_INLINE __m64 _mm_sign_pi8(__m64 _a, __m64 _b) {}

/* SSE4.1 */

// Blend packed 16-bit integers from a and b using control mask imm8, and store
// the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_blend_epi16
// FORCE_INLINE __m128i _mm_blend_epi16(__m128i a, __m128i b,
//                                      __constrange(0,255) int imm)
// #define _mm_blend_epi16(a, b, imm)

// Blend packed double-precision (64-bit) floating-point elements from a and b
// using control mask imm8, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_blend_pd
// #define _mm_blend_pd(a, b, imm)

// Blend packed single-precision (32-bit) floating-point elements from a and b
// using mask, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_blend_ps
// FORCE_INLINE __m128 _mm_blend_ps(__m128 _a, __m128 _b, const char imm8) {}

// Blend packed 8-bit integers from a and b using mask, and store the results in
// dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_blendv_epi8
// FORCE_INLINE __m128i _mm_blendv_epi8(__m128i _a, __m128i _b, __m128i _mask) {}

// Blend packed double-precision (64-bit) floating-point elements from a and b
// using mask, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_blendv_pd
// FORCE_INLINE __m128d _mm_blendv_pd(__m128d _a, __m128d _b, __m128d _mask) {}

// Blend packed single-precision (32-bit) floating-point elements from a and b
// using mask, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_blendv_ps
// FORCE_INLINE __m128 _mm_blendv_ps(__m128 _a, __m128 _b, __m128 _mask) {}

// Round the packed double-precision (64-bit) floating-point elements in a up
// to an integer value, and store the results as packed double-precision
// floating-point elements in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_ceil_pd
// FORCE_INLINE __m128d _mm_ceil_pd(__m128d a) {}

// Round the packed single-precision (32-bit) floating-point elements in a up to
// an integer value, and store the results as packed single-precision
// floating-point elements in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_ceil_ps
// FORCE_INLINE __m128 _mm_ceil_ps(__m128 a) {}

// Round the lower double-precision (64-bit) floating-point element in b up to
// an integer value, store the result as a double-precision floating-point
// element in the lower element of dst, and copy the upper element from a to the
// upper element of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_ceil_sd
// FORCE_INLINE __m128d _mm_ceil_sd(__m128d a, __m128d b) {}

// Round the lower single-precision (32-bit) floating-point element in b up to
// an integer value, store the result as a single-precision floating-point
// element in the lower element of dst, and copy the upper 3 packed elements
// from a to the upper elements of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_ceil_ss
// FORCE_INLINE __m128 _mm_ceil_ss(__m128 a, __m128 b) {}

// Compare packed 64-bit integers in a and b for equality, and store the results
// in dst
// FORCE_INLINE __m128i _mm_cmpeq_epi64(__m128i a, __m128i b) {}

// Sign extend packed 16-bit integers in a to packed 32-bit integers, and store
// the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtepi16_epi32
// FORCE_INLINE __m128i _mm_cvtepi16_epi32(__m128i a) {}

// Sign extend packed 16-bit integers in a to packed 64-bit integers, and store
// the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtepi16_epi64
// FORCE_INLINE __m128i _mm_cvtepi16_epi64(__m128i a) {}

// Sign extend packed 32-bit integers in a to packed 64-bit integers, and store
// the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtepi32_epi64
// FORCE_INLINE __m128i _mm_cvtepi32_epi64(__m128i a) {}

// Sign extend packed 8-bit integers in a to packed 16-bit integers, and store
// the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtepi8_epi16
// FORCE_INLINE __m128i _mm_cvtepi8_epi16(__m128i a) {}

// Sign extend packed 8-bit integers in a to packed 32-bit integers, and store
// the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtepi8_epi32
// FORCE_INLINE __m128i _mm_cvtepi8_epi32(__m128i a) {}

// Sign extend packed 8-bit integers in the low 8 bytes of a to packed 64-bit
// integers, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtepi8_epi64
// FORCE_INLINE __m128i _mm_cvtepi8_epi64(__m128i a) {}

// Zero extend packed unsigned 16-bit integers in a to packed 32-bit integers,
// and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtepu16_epi32
// FORCE_INLINE __m128i _mm_cvtepu16_epi32(__m128i a) {}

// Zero extend packed unsigned 16-bit integers in a to packed 64-bit integers,
// and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtepu16_epi64
// FORCE_INLINE __m128i _mm_cvtepu16_epi64(__m128i a) {}

// Zero extend packed unsigned 32-bit integers in a to packed 64-bit integers,
// and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtepu32_epi64
// FORCE_INLINE __m128i _mm_cvtepu32_epi64(__m128i a) {}

// Zero extend packed unsigned 8-bit integers in a to packed 16-bit integers,
// and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtepu8_epi16
// FORCE_INLINE __m128i _mm_cvtepu8_epi16(__m128i a) {}

// Zero extend packed unsigned 8-bit integers in a to packed 32-bit integers,
// and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtepu8_epi32
// FORCE_INLINE __m128i _mm_cvtepu8_epi32(__m128i a) {}

// Zero extend packed unsigned 8-bit integers in the low 8 bytes of a to packed
// 64-bit integers, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cvtepu8_epi64
// FORCE_INLINE __m128i _mm_cvtepu8_epi64(__m128i a) {}

// Conditionally multiply the packed double-precision (64-bit) floating-point
// elements in a and b using the high 4 bits in imm8, sum the four products, and
// conditionally store the sum in dst using the low 4 bits of imm8.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_dp_pd
// FORCE_INLINE __m128d _mm_dp_pd(__m128d a, __m128d b, const int imm) {}

// Conditionally multiply the packed single-precision (32-bit) floating-point
// elements in a and b using the high 4 bits in imm8, sum the four products,
// and conditionally store the sum in dst using the low 4 bits of imm.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_dp_ps
// FORCE_INLINE __m128 _mm_dp_ps(__m128 a, __m128 b, const int imm) {}

// Extract a 32-bit integer from a, selected with imm8, and store the result in
// dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_extract_epi32
// FORCE_INLINE int _mm_extract_epi32(__m128i a, __constrange(0,4) int imm)
// #define _mm_extract_epi32(a, imm)

// Extract a 64-bit integer from a, selected with imm8, and store the result in
// dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_extract_epi64
// FORCE_INLINE __int64 _mm_extract_epi64(__m128i a, __constrange(0,2) int imm)
// #define _mm_extract_epi64(a, imm)

// Extract an 8-bit integer from a, selected with imm8, and store the result in
// the lower element of dst. FORCE_INLINE int _mm_extract_epi8(__m128i a,
// __constrange(0,16) int imm)
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_extract_epi8
// #define _mm_extract_epi8(a, imm) vgetq_lane_u8(vreinterpretq_u8_m128i(a), (imm))

// Extracts the selected single-precision (32-bit) floating-point from a.
// FORCE_INLINE int _mm_extract_ps(__m128 a, __constrange(0,4) int imm)
// #define _mm_extract_ps(a, imm) vgetq_lane_s32(vreinterpretq_s32_m128(a), (imm))

// Round the packed double-precision (64-bit) floating-point elements in a down
// to an integer value, and store the results as packed double-precision
// floating-point elements in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_floor_pd
// FORCE_INLINE __m128d _mm_floor_pd(__m128d a) {}

// Round the packed single-precision (32-bit) floating-point elements in a down
// to an integer value, and store the results as packed single-precision
// floating-point elements in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_floor_ps
// FORCE_INLINE __m128 _mm_floor_ps(__m128 a) {}

// Round the lower double-precision (64-bit) floating-point element in b down to
// an integer value, store the result as a double-precision floating-point
// element in the lower element of dst, and copy the upper element from a to the
// upper element of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_floor_sd
// FORCE_INLINE __m128d _mm_floor_sd(__m128d a, __m128d b) {}

// Round the lower single-precision (32-bit) floating-point element in b down to
// an integer value, store the result as a single-precision floating-point
// element in the lower element of dst, and copy the upper 3 packed elements
// from a to the upper elements of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_floor_ss
// FORCE_INLINE __m128 _mm_floor_ss(__m128 a, __m128 b) {}

// Copy a to dst, and insert the 32-bit integer i into dst at the location
// specified by imm8.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_insert_epi32
// FORCE_INLINE __m128i _mm_insert_epi32(__m128i a, int b,
//                                       __constrange(0,4) int imm)
// #define _mm_insert_epi32(a, b, imm)

// Copy a to dst, and insert the 64-bit integer i into dst at the location
// specified by imm8.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_insert_epi64
// FORCE_INLINE __m128i _mm_insert_epi64(__m128i a, __int64 b,
//                                       __constrange(0,2) int imm)
// #define _mm_insert_epi64(a, b, imm)

// Copy a to dst, and insert the lower 8-bit integer from i into dst at the
// location specified by imm8.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_insert_epi8
// FORCE_INLINE __m128i _mm_insert_epi8(__m128i a, int b,
//                                      __constrange(0,16) int imm)
// #define _mm_insert_epi8(a, b, imm)

// Copy a to tmp, then insert a single-precision (32-bit) floating-point
// element from b into tmp using the control in imm8. Store tmp to dst using
// the mask in imm8 (elements are zeroed out when the corresponding bit is set).
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=insert_ps
// #define _mm_insert_ps(a, b, imm8)

// Compare packed signed 32-bit integers in a and b, and store packed maximum
// values in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_max_epi32
// FORCE_INLINE __m128i _mm_max_epi32(__m128i a, __m128i b) {}

// Compare packed signed 8-bit integers in a and b, and store packed maximum
// values in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_max_epi8
// FORCE_INLINE __m128i _mm_max_epi8(__m128i a, __m128i b) {}

// Compare packed unsigned 16-bit integers in a and b, and store packed maximum
// values in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_max_epu16
// FORCE_INLINE __m128i _mm_max_epu16(__m128i a, __m128i b) {}

// Compare packed unsigned 32-bit integers in a and b, and store packed maximum
// values in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_max_epu32
// FORCE_INLINE __m128i _mm_max_epu32(__m128i a, __m128i b) {}

// Compare packed signed 32-bit integers in a and b, and store packed minimum
// values in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_min_epi32
// FORCE_INLINE __m128i _mm_min_epi32(__m128i a, __m128i b) {}

// Compare packed signed 8-bit integers in a and b, and store packed minimum
// values in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_min_epi8
// FORCE_INLINE __m128i _mm_min_epi8(__m128i a, __m128i b) {}

// Compare packed unsigned 16-bit integers in a and b, and store packed minimum
// values in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_min_epu16
// FORCE_INLINE __m128i _mm_min_epu16(__m128i a, __m128i b) {}

// Compare packed unsigned 32-bit integers in a and b, and store packed minimum
// values in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_max_epu32
// FORCE_INLINE __m128i _mm_min_epu32(__m128i a, __m128i b) {}

// Horizontally compute the minimum amongst the packed unsigned 16-bit integers
// in a, store the minimum and index in dst, and zero the remaining bits in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_minpos_epu16
// FORCE_INLINE __m128i _mm_minpos_epu16(__m128i a) {}

// Compute the sum of absolute differences (SADs) of quadruplets of unsigned
// 8-bit integers in a compared to those in b, and store the 16-bit results in
// dst. Eight SADs are performed using one quadruplet from b and eight
// quadruplets from a. One quadruplet is selected from b starting at on the
// offset specified in imm8. Eight quadruplets are formed from sequential 8-bit
// integers selected from a starting at the offset specified in imm8.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_mpsadbw_epu8
// FORCE_INLINE __m128i _mm_mpsadbw_epu8(__m128i a, __m128i b, const int imm) {}

// Multiply the low signed 32-bit integers from each packed 64-bit element in
// a and b, and store the signed 64-bit results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_mul_epi32
// FORCE_INLINE __m128i _mm_mul_epi32(__m128i a, __m128i b) {}

// Multiply the packed 32-bit integers in a and b, producing intermediate 64-bit
// integers, and store the low 32 bits of the intermediate integers in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_mullo_epi32
// FORCE_INLINE __m128i _mm_mullo_epi32(__m128i a, __m128i b) {}

// Convert packed signed 32-bit integers from a and b to packed 16-bit integers
// using unsigned saturation, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_packus_epi32
// FORCE_INLINE __m128i _mm_packus_epi32(__m128i a, __m128i b) {}

// Round the packed double-precision (64-bit) floating-point elements in a using
// the rounding parameter, and store the results as packed double-precision
// floating-point elements in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_round_pd
// FORCE_INLINE __m128d _mm_round_pd(__m128d a, int rounding) {}

// Round the packed single-precision (32-bit) floating-point elements in a using
// the rounding parameter, and store the results as packed single-precision
// floating-point elements in dst.
// software.intel.com/sites/landingpage/IntrinsicsGuide/#text=_mm_round_ps
// FORCE_INLINE __m128 _mm_round_ps(__m128 a, int rounding) {}

// Round the lower double-precision (64-bit) floating-point element in b using
// the rounding parameter, store the result as a double-precision floating-point
// element in the lower element of dst, and copy the upper element from a to the
// upper element of dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_round_sd
// FORCE_INLINE __m128d _mm_round_sd(__m128d a, __m128d b, int rounding) {}

// Round the lower single-precision (32-bit) floating-point element in b using
// the rounding parameter, store the result as a single-precision floating-point
// element in the lower element of dst, and copy the upper 3 packed elements
// from a to the upper elements of dst. Rounding is done according to the
// rounding[3:0] parameter, which can be one of:
//     (_MM_FROUND_TO_NEAREST_INT |_MM_FROUND_NO_EXC) // round to nearest, and
//     suppress exceptions
//     (_MM_FROUND_TO_NEG_INF |_MM_FROUND_NO_EXC)     // round down, and
//     suppress exceptions
//     (_MM_FROUND_TO_POS_INF |_MM_FROUND_NO_EXC)     // round up, and suppress
//     exceptions
//     (_MM_FROUND_TO_ZERO |_MM_FROUND_NO_EXC)        // truncate, and suppress
//     exceptions _MM_FROUND_CUR_DIRECTION // use MXCSR.RC; see
//     _MM_SET_ROUNDING_MODE
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_round_ss
// FORCE_INLINE __m128 _mm_round_ss(__m128 a, __m128 b, int rounding) {}

// Load 128-bits of integer data from memory into dst using a non-temporal
// memory hint. mem_addr must be aligned on a 16-byte boundary or a
// general-protection exception may be generated.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_stream_load_si128
// FORCE_INLINE __m128i _mm_stream_load_si128(__m128i *p) {}

// Compute the bitwise NOT of a and then AND with a 128-bit vector containing
// all 1's, and return 1 if the result is zero, otherwise return 0.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_test_all_ones
// FORCE_INLINE int _mm_test_all_ones(__m128i a) {}

// Compute the bitwise AND of 128 bits (representing integer data) in a and
// mask, and return 1 if the result is zero, otherwise return 0.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_test_all_zeros
// FORCE_INLINE int _mm_test_all_zeros(__m128i a, __m128i mask) {}

// Compute the bitwise AND of 128 bits (representing integer data) in a and
// mask, and set ZF to 1 if the result is zero, otherwise set ZF to 0. Compute
// the bitwise NOT of a and then AND with mask, and set CF to 1 if the result is
// zero, otherwise set CF to 0. Return 1 if both the ZF and CF values are zero,
// otherwise return 0.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=mm_test_mix_ones_zero
// FORCE_INLINE int _mm_test_mix_ones_zeros(__m128i a, __m128i mask) {}

// Compute the bitwise AND of 128 bits (representing integer data) in a and b,
// and set ZF to 1 if the result is zero, otherwise set ZF to 0. Compute the
// bitwise NOT of a and then AND with b, and set CF to 1 if the result is zero,
// otherwise set CF to 0. Return the CF value.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_testc_si128
// FORCE_INLINE int _mm_testc_si128(__m128i a, __m128i b) {}

// Compute the bitwise AND of 128 bits (representing integer data) in a and b,
// and set ZF to 1 if the result is zero, otherwise set ZF to 0. Compute the
// bitwise NOT of a and then AND with b, and set CF to 1 if the result is zero,
// otherwise set CF to 0. Return 1 if both the ZF and CF values are zero,
// otherwise return 0.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_testnzc_si128
// #define _mm_testnzc_si128(a, b) _mm_test_mix_ones_zeros(a, b)

// Compute the bitwise AND of 128 bits (representing integer data) in a and b,
// and set ZF to 1 if the result is zero, otherwise set ZF to 0. Compute the
// bitwise NOT of a and then AND with b, and set CF to 1 if the result is zero,
// otherwise set CF to 0. Return the ZF value.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_testz_si128
// FORCE_INLINE int _mm_testz_si128(__m128i a, __m128i b) {}

/* SSE4.2 */

static const uint16_t ALIGN_STRUCT(16) _sse2rvv_cmpestr_mask16b[8] = {
    0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80,
};
static const uint8_t ALIGN_STRUCT(16) _sse2rvv_cmpestr_mask8b[16] = {
    0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80,
    0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80,
};

/* specify the source data format */
#define _SIDD_UBYTE_OPS 0x00 /* unsigned 8-bit characters */
#define _SIDD_UWORD_OPS 0x01 /* unsigned 16-bit characters */
#define _SIDD_SBYTE_OPS 0x02 /* signed 8-bit characters */
#define _SIDD_SWORD_OPS 0x03 /* signed 16-bit characters */

/* specify the comparison operation */
#define _SIDD_CMP_EQUAL_ANY 0x00     /* compare equal any: strchr */
#define _SIDD_CMP_RANGES 0x04        /* compare ranges */
#define _SIDD_CMP_EQUAL_EACH 0x08    /* compare equal each: strcmp */
#define _SIDD_CMP_EQUAL_ORDERED 0x0C /* compare equal ordered */

/* specify the polarity */
#define _SIDD_POSITIVE_POLARITY 0x00
#define _SIDD_MASKED_POSITIVE_POLARITY 0x20
#define _SIDD_NEGATIVE_POLARITY 0x10 /* negate results */
#define _SIDD_MASKED_NEGATIVE_POLARITY \
    0x30 /* negate results only before end of string */

/* specify the output selection in _mm_cmpXstri */
#define _SIDD_LEAST_SIGNIFICANT 0x00
#define _SIDD_MOST_SIGNIFICANT 0x40

/* specify the output selection in _mm_cmpXstrm */
#define _SIDD_BIT_MASK 0x00
#define _SIDD_UNIT_MASK 0x40

/* Pattern Matching for C macros.
 * https://github.com/pfultz2/Cloak/wiki/C-Preprocessor-tricks,-tips,-and-idioms
 */

/* catenate */
#define SSE2RVV_PRIMITIVE_CAT(a, ...) a##__VA_ARGS__
#define SSE2RVV_CAT(a, b) SSE2RVV_PRIMITIVE_CAT(a, b)

#define SSE2RVV_IIF(c) SSE2RVV_PRIMITIVE_CAT(SSE2RVV_IIF_, c)
/* run the 2nd parameter */
#define SSE2RVV_IIF_0(t, ...) __VA_ARGS__
/* run the 1st parameter */
#define SSE2RVV_IIF_1(t, ...) t

#define SSE2RVV_COMPL(b) SSE2RVV_PRIMITIVE_CAT(SSE2RVV_COMPL_, b)
#define SSE2RVV_COMPL_0 1
#define SSE2RVV_COMPL_1 0

#define SSE2RVV_DEC(x) SSE2RVV_PRIMITIVE_CAT(SSE2RVV_DEC_, x)
#define SSE2RVV_DEC_1 0
#define SSE2RVV_DEC_2 1
#define SSE2RVV_DEC_3 2
#define SSE2RVV_DEC_4 3
#define SSE2RVV_DEC_5 4
#define SSE2RVV_DEC_6 5
#define SSE2RVV_DEC_7 6
#define SSE2RVV_DEC_8 7
#define SSE2RVV_DEC_9 8
#define SSE2RVV_DEC_10 9
#define SSE2RVV_DEC_11 10
#define SSE2RVV_DEC_12 11
#define SSE2RVV_DEC_13 12
#define SSE2RVV_DEC_14 13
#define SSE2RVV_DEC_15 14
#define SSE2RVV_DEC_16 15

/* detection */
#define SSE2RVV_CHECK_N(x, n, ...) n
#define SSE2RVV_CHECK(...) SSE2RVV_CHECK_N(__VA_ARGS__, 0, )
#define SSE2RVV_PROBE(x) x, 1,

#define SSE2RVV_NOT(x) SSE2RVV_CHECK(SSE2RVV_PRIMITIVE_CAT(SSE2RVV_NOT_, x))
#define SSE2RVV_NOT_0 SSE2RVV_PROBE(~)

#define SSE2RVV_BOOL(x) SSE2RVV_COMPL(SSE2RVV_NOT(x))
#define SSE2RVV_IF(c) SSE2RVV_IIF(SSE2RVV_BOOL(c))

#define SSE2RVV_EAT(...)
#define SSE2RVV_EXPAND(...) __VA_ARGS__
#define SSE2RVV_WHEN(c) SSE2RVV_IF(c)(SSE2RVV_EXPAND, SSE2RVV_EAT)

/* recursion */
/* deferred expression */
#define SSE2RVV_EMPTY()
#define SSE2RVV_DEFER(id) id SSE2RVV_EMPTY()
#define SSE2RVV_OBSTRUCT(...) __VA_ARGS__ SSE2RVV_DEFER(SSE2RVV_EMPTY)()
#define SSE2RVV_EXPAND(...) __VA_ARGS__

#define SSE2RVV_EVAL(...) \
    SSE2RVV_EVAL1(SSE2RVV_EVAL1(SSE2RVV_EVAL1(__VA_ARGS__)))
#define SSE2RVV_EVAL1(...) \
    SSE2RVV_EVAL2(SSE2RVV_EVAL2(SSE2RVV_EVAL2(__VA_ARGS__)))
#define SSE2RVV_EVAL2(...) \
    SSE2RVV_EVAL3(SSE2RVV_EVAL3(SSE2RVV_EVAL3(__VA_ARGS__)))
#define SSE2RVV_EVAL3(...) __VA_ARGS__

#define SSE2RVV_REPEAT(count, macro, ...)         \
    SSE2RVV_WHEN(count)                           \
    (SSE2RVV_OBSTRUCT(SSE2RVV_REPEAT_INDIRECT)()( \
        SSE2RVV_DEC(count), macro,                \
        __VA_ARGS__) SSE2RVV_OBSTRUCT(macro)(SSE2RVV_DEC(count), __VA_ARGS__))
#define SSE2RVV_REPEAT_INDIRECT() SSE2RVV_REPEAT

#define SSE2RVV_SIZE_OF_byte 8
#define SSE2RVV_NUMBER_OF_LANES_byte 16
#define SSE2RVV_SIZE_OF_word 16
#define SSE2RVV_NUMBER_OF_LANES_word 8

#define SSE2RVV_COMPARE_EQUAL_THEN_FILL_LANE(i, type)                          \
    mtx[i] = vreinterpretq_m128i_##type(vceqq_##type(                          \
        vdupq_n_##type(vgetq_lane_##type(vreinterpretq_##type##_m128i(b), i)), \
        vreinterpretq_##type##_m128i(a)));

#define SSE2RVV_FILL_LANE(i, type) \
    vec_b[i] =                     \
        vdupq_n_##type(vgetq_lane_##type(vreinterpretq_##type##_m128i(b), i));

#define PCMPSTR_RANGES(a, b, mtx, data_type_prefix, type_prefix, size,         \
                       number_of_lanes, byte_or_word)                          \
    do {                                                                       \
        SSE2RVV_CAT(                                                           \
            data_type_prefix,                                                  \
            SSE2RVV_CAT(size,                                                  \
                        SSE2RVV_CAT(x, SSE2RVV_CAT(number_of_lanes, _t))))     \
        vec_b[number_of_lanes];                                                \
        __m128i mask = SSE2RVV_IIF(byte_or_word)(                              \
            vreinterpretq_m128i_u16(vdupq_n_u16(0xff)),                        \
            vreinterpretq_m128i_u32(vdupq_n_u32(0xffff)));                     \
        SSE2RVV_EVAL(SSE2RVV_REPEAT(number_of_lanes, SSE2RVV_FILL_LANE,        \
                                    SSE2RVV_CAT(type_prefix, size)))           \
        for (int i = 0; i < number_of_lanes; i++) {                            \
            mtx[i] = SSE2RVV_CAT(vreinterpretq_m128i_u,                        \
                                 size)(SSE2RVV_CAT(vbslq_u, size)(             \
                SSE2RVV_CAT(vreinterpretq_u, SSE2RVV_CAT(size, _m128i))(mask), \
                SSE2RVV_CAT(vcgeq_, SSE2RVV_CAT(type_prefix, size))(           \
                    vec_b[i],                                                  \
                    SSE2RVV_CAT(vreinterpretq_,                                \
                                SSE2RVV_CAT(type_prefix,                       \
                                            SSE2RVV_CAT(size, _m128i(a))))),   \
                SSE2RVV_CAT(vcleq_, SSE2RVV_CAT(type_prefix, size))(           \
                    vec_b[i],                                                  \
                    SSE2RVV_CAT(vreinterpretq_,                                \
                                SSE2RVV_CAT(type_prefix,                       \
                                            SSE2RVV_CAT(size, _m128i(a))))))); \
        }                                                                      \
    } while (0)

#define PCMPSTR_EQ(a, b, mtx, size, number_of_lanes)                      \
    do {                                                                  \
        SSE2RVV_EVAL(SSE2RVV_REPEAT(number_of_lanes,                      \
                                    SSE2RVV_COMPARE_EQUAL_THEN_FILL_LANE, \
                                    SSE2RVV_CAT(u, size)))                \
    } while (0)

#define SSE2RVV_CMP_EQUAL_ANY_IMPL(type)                                     \
    static int _sse2rvv_cmp_##type##_equal_any(__m128i a, int la, __m128i b, \
                                               int lb)                       \
    {                                                                        \
        __m128i mtx[16];                                                     \
        PCMPSTR_EQ(a, b, mtx, SSE2RVV_CAT(SSE2RVV_SIZE_OF_, type),           \
                   SSE2RVV_CAT(SSE2RVV_NUMBER_OF_LANES_, type));             \
        return SSE2RVV_CAT(                                                  \
            _sse2rvv_aggregate_equal_any_,                                   \
            SSE2RVV_CAT(SSE2RVV_CAT(SSE2RVV_SIZE_OF_, type),                 \
                        SSE2RVV_CAT(x, SSE2RVV_CAT(SSE2RVV_NUMBER_OF_LANES_, \
                                                   type))))(la, lb, mtx);    \
    }

#define SSE2RVV_CMP_RANGES_IMPL(type, data_type, us, byte_or_word)            \
    static int _sse2rvv_cmp_##us##type##_ranges(__m128i a, int la, __m128i b, \
                                                int lb)                       \
    {                                                                         \
        __m128i mtx[16];                                                      \
        PCMPSTR_RANGES(                                                       \
            a, b, mtx, data_type, us, SSE2RVV_CAT(SSE2RVV_SIZE_OF_, type),    \
            SSE2RVV_CAT(SSE2RVV_NUMBER_OF_LANES_, type), byte_or_word);       \
        return SSE2RVV_CAT(                                                   \
            _sse2rvv_aggregate_ranges_,                                       \
            SSE2RVV_CAT(SSE2RVV_CAT(SSE2RVV_SIZE_OF_, type),                  \
                        SSE2RVV_CAT(x, SSE2RVV_CAT(SSE2RVV_NUMBER_OF_LANES_,  \
                                                   type))))(la, lb, mtx);     \
    }

#define SSE2RVV_CMP_EQUAL_ORDERED_IMPL(type)                                   \
    static int _sse2rvv_cmp_##type##_equal_ordered(__m128i a, int la,          \
                                                   __m128i b, int lb)          \
    {                                                                          \
        __m128i mtx[16];                                                       \
        PCMPSTR_EQ(a, b, mtx, SSE2RVV_CAT(SSE2RVV_SIZE_OF_, type),             \
                   SSE2RVV_CAT(SSE2RVV_NUMBER_OF_LANES_, type));               \
        return SSE2RVV_CAT(                                                    \
            _sse2rvv_aggregate_equal_ordered_,                                 \
            SSE2RVV_CAT(                                                       \
                SSE2RVV_CAT(SSE2RVV_SIZE_OF_, type),                           \
                SSE2RVV_CAT(x, SSE2RVV_CAT(SSE2RVV_NUMBER_OF_LANES_, type))))( \
            SSE2RVV_CAT(SSE2RVV_NUMBER_OF_LANES_, type), la, lb, mtx);         \
    }

static int _sse2rvv_aggregate_equal_any_8x16(int la, int lb, __m128i mtx[16])
{
    int res = 0;
    int m = (1 << la) - 1;
    uint8x8_t vec_mask = vld1_u8(_sse2rvv_cmpestr_mask8b);
    uint8x8_t t_lo = vtst_u8(vdup_n_u8(m & 0xff), vec_mask);
    uint8x8_t t_hi = vtst_u8(vdup_n_u8(m >> 8), vec_mask);
    uint8x16_t vec = vcombine_u8(t_lo, t_hi);
    for (int j = 0; j < lb; j++) {
        mtx[j] = vreinterpretq_m128i_u8(
            vandq_u8(vec, vreinterpretq_u8_m128i(mtx[j])));
        mtx[j] = vreinterpretq_m128i_u8(
            vshrq_n_u8(vreinterpretq_u8_m128i(mtx[j]), 7));
        int tmp = _sse2rvv_vaddvq_u8(vreinterpretq_u8_m128i(mtx[j])) ? 1 : 0;
        res |= (tmp << j);
    }
    return res;
}

static int _sse2rvv_aggregate_equal_any_16x8(int la, int lb, __m128i mtx[16])
{
    int res = 0;
    int m = (1 << la) - 1;
    uint16x8_t vec =
        vtstq_u16(vdupq_n_u16(m), vld1q_u16(_sse2rvv_cmpestr_mask16b));
    for (int j = 0; j < lb; j++) {
        mtx[j] = vreinterpretq_m128i_u16(
            vandq_u16(vec, vreinterpretq_u16_m128i(mtx[j])));
        mtx[j] = vreinterpretq_m128i_u16(
            vshrq_n_u16(vreinterpretq_u16_m128i(mtx[j]), 15));
        int tmp = _sse2rvv_vaddvq_u16(vreinterpretq_u16_m128i(mtx[j])) ? 1 : 0;
        res |= (tmp << j);
    }
    return res;
}

/* clang-format off */
#define SSE2RVV_GENERATE_CMP_EQUAL_ANY(prefix) \
    prefix##IMPL(byte) \
    prefix##IMPL(word)
/* clang-format on */

SSE2RVV_GENERATE_CMP_EQUAL_ANY(SSE2RVV_CMP_EQUAL_ANY_)

static int _sse2rvv_aggregate_ranges_16x8(int la, int lb, __m128i mtx[16])
{
    int res = 0;
    int m = (1 << la) - 1;
    uint16x8_t vec =
        vtstq_u16(vdupq_n_u16(m), vld1q_u16(_sse2rvv_cmpestr_mask16b));
    for (int j = 0; j < lb; j++) {
        mtx[j] = vreinterpretq_m128i_u16(
            vandq_u16(vec, vreinterpretq_u16_m128i(mtx[j])));
        mtx[j] = vreinterpretq_m128i_u16(
            vshrq_n_u16(vreinterpretq_u16_m128i(mtx[j]), 15));
        __m128i tmp = vreinterpretq_m128i_u32(
            vshrq_n_u32(vreinterpretq_u32_m128i(mtx[j]), 16));
        uint32x4_t vec_res = vandq_u32(vreinterpretq_u32_m128i(mtx[j]),
                                       vreinterpretq_u32_m128i(tmp));
#if defined(__aarch64__) || defined(_M_ARM64)
        int t = vaddvq_u32(vec_res) ? 1 : 0;
#else
        uint64x2_t sumh = vpaddlq_u32(vec_res);
        int t = vgetq_lane_u64(sumh, 0) + vgetq_lane_u64(sumh, 1);
#endif
        res |= (t << j);
    }
    return res;
}

static int _sse2rvv_aggregate_ranges_8x16(int la, int lb, __m128i mtx[16])
{
    int res = 0;
    int m = (1 << la) - 1;
    uint8x8_t vec_mask = vld1_u8(_sse2rvv_cmpestr_mask8b);
    uint8x8_t t_lo = vtst_u8(vdup_n_u8(m & 0xff), vec_mask);
    uint8x8_t t_hi = vtst_u8(vdup_n_u8(m >> 8), vec_mask);
    uint8x16_t vec = vcombine_u8(t_lo, t_hi);
    for (int j = 0; j < lb; j++) {
        mtx[j] = vreinterpretq_m128i_u8(
            vandq_u8(vec, vreinterpretq_u8_m128i(mtx[j])));
        mtx[j] = vreinterpretq_m128i_u8(
            vshrq_n_u8(vreinterpretq_u8_m128i(mtx[j]), 7));
        __m128i tmp = vreinterpretq_m128i_u16(
            vshrq_n_u16(vreinterpretq_u16_m128i(mtx[j]), 8));
        uint16x8_t vec_res = vandq_u16(vreinterpretq_u16_m128i(mtx[j]),
                                       vreinterpretq_u16_m128i(tmp));
        int t = _sse2rvv_vaddvq_u16(vec_res) ? 1 : 0;
        res |= (t << j);
    }
    return res;
}

#define SSE2RVV_CMP_RANGES_IS_BYTE 1
#define SSE2RVV_CMP_RANGES_IS_WORD 0

/* clang-format off */
#define SSE2RVV_GENERATE_CMP_RANGES(prefix)             \
    prefix##IMPL(byte, uint, u, prefix##IS_BYTE)         \
    prefix##IMPL(byte, int, s, prefix##IS_BYTE)          \
    prefix##IMPL(word, uint, u, prefix##IS_WORD)         \
    prefix##IMPL(word, int, s, prefix##IS_WORD)
/* clang-format on */

SSE2RVV_GENERATE_CMP_RANGES(SSE2RVV_CMP_RANGES_)

#undef SSE2RVV_CMP_RANGES_IS_BYTE
#undef SSE2RVV_CMP_RANGES_IS_WORD

static int _sse2rvv_cmp_byte_equal_each(__m128i a, int la, __m128i b, int lb)
{
    uint8x16_t mtx =
        vceqq_u8(vreinterpretq_u8_m128i(a), vreinterpretq_u8_m128i(b));
    int m0 = (la < lb) ? 0 : ((1 << la) - (1 << lb));
    int m1 = 0x10000 - (1 << la);
    int tb = 0x10000 - (1 << lb);
    uint8x8_t vec_mask, vec0_lo, vec0_hi, vec1_lo, vec1_hi;
    uint8x8_t tmp_lo, tmp_hi, res_lo, res_hi;
    vec_mask = vld1_u8(_sse2rvv_cmpestr_mask8b);
    vec0_lo = vtst_u8(vdup_n_u8(m0), vec_mask);
    vec0_hi = vtst_u8(vdup_n_u8(m0 >> 8), vec_mask);
    vec1_lo = vtst_u8(vdup_n_u8(m1), vec_mask);
    vec1_hi = vtst_u8(vdup_n_u8(m1 >> 8), vec_mask);
    tmp_lo = vtst_u8(vdup_n_u8(tb), vec_mask);
    tmp_hi = vtst_u8(vdup_n_u8(tb >> 8), vec_mask);

    res_lo = vbsl_u8(vec0_lo, vdup_n_u8(0), vget_low_u8(mtx));
    res_hi = vbsl_u8(vec0_hi, vdup_n_u8(0), vget_high_u8(mtx));
    res_lo = vbsl_u8(vec1_lo, tmp_lo, res_lo);
    res_hi = vbsl_u8(vec1_hi, tmp_hi, res_hi);
    res_lo = vand_u8(res_lo, vec_mask);
    res_hi = vand_u8(res_hi, vec_mask);

    int res = _sse2rvv_vaddv_u8(res_lo) + (_sse2rvv_vaddv_u8(res_hi) << 8);
    return res;
}

static int _sse2rvv_cmp_word_equal_each(__m128i a, int la, __m128i b, int lb)
{
    uint16x8_t mtx =
        vceqq_u16(vreinterpretq_u16_m128i(a), vreinterpretq_u16_m128i(b));
    int m0 = (la < lb) ? 0 : ((1 << la) - (1 << lb));
    int m1 = 0x100 - (1 << la);
    int tb = 0x100 - (1 << lb);
    uint16x8_t vec_mask = vld1q_u16(_sse2rvv_cmpestr_mask16b);
    uint16x8_t vec0 = vtstq_u16(vdupq_n_u16(m0), vec_mask);
    uint16x8_t vec1 = vtstq_u16(vdupq_n_u16(m1), vec_mask);
    uint16x8_t tmp = vtstq_u16(vdupq_n_u16(tb), vec_mask);
    mtx = vbslq_u16(vec0, vdupq_n_u16(0), mtx);
    mtx = vbslq_u16(vec1, tmp, mtx);
    mtx = vandq_u16(mtx, vec_mask);
    return _sse2rvv_vaddvq_u16(mtx);
}

#define SSE2RVV_AGGREGATE_EQUAL_ORDER_IS_UBYTE 1
#define SSE2RVV_AGGREGATE_EQUAL_ORDER_IS_UWORD 0

#define SSE2RVV_AGGREGATE_EQUAL_ORDER_IMPL(size, number_of_lanes, data_type)   \
    static int _sse2rvv_aggregate_equal_ordered_##size##x##number_of_lanes(    \
        int bound, int la, int lb, __m128i mtx[16])                            \
    {                                                                          \
        int res = 0;                                                           \
        int m1 = SSE2RVV_IIF(data_type)(0x10000, 0x100) - (1 << la);           \
        uint##size##x8_t vec_mask = SSE2RVV_IIF(data_type)(                    \
            vld1_u##size(_sse2rvv_cmpestr_mask##size##b),                      \
            vld1q_u##size(_sse2rvv_cmpestr_mask##size##b));                    \
        uint##size##x##number_of_lanes##_t vec1 = SSE2RVV_IIF(data_type)(      \
            vcombine_u##size(vtst_u##size(vdup_n_u##size(m1), vec_mask),       \
                             vtst_u##size(vdup_n_u##size(m1 >> 8), vec_mask)), \
            vtstq_u##size(vdupq_n_u##size(m1), vec_mask));                     \
        uint##size##x##number_of_lanes##_t vec_minusone = vdupq_n_u##size(-1); \
        uint##size##x##number_of_lanes##_t vec_zero = vdupq_n_u##size(0);      \
        for (int j = 0; j < lb; j++) {                                         \
            mtx[j] = vreinterpretq_m128i_u##size(vbslq_u##size(                \
                vec1, vec_minusone, vreinterpretq_u##size##_m128i(mtx[j])));   \
        }                                                                      \
        for (int j = lb; j < bound; j++) {                                     \
            mtx[j] = vreinterpretq_m128i_u##size(                              \
                vbslq_u##size(vec1, vec_minusone, vec_zero));                  \
        }                                                                      \
        unsigned SSE2RVV_IIF(data_type)(char, short) *ptr =                    \
            (unsigned SSE2RVV_IIF(data_type)(char, short) *) mtx;              \
        for (int i = 0; i < bound; i++) {                                      \
            int val = 1;                                                       \
            for (int j = 0, k = i; j < bound - i && k < bound; j++, k++)       \
                val &= ptr[k * bound + j];                                     \
            res += val << i;                                                   \
        }                                                                      \
        return res;                                                            \
    }

/* clang-format off */
#define SSE2RVV_GENERATE_AGGREGATE_EQUAL_ORDER(prefix) \
    prefix##IMPL(8, 16, prefix##IS_UBYTE)               \
    prefix##IMPL(16, 8, prefix##IS_UWORD)
/* clang-format on */

SSE2RVV_GENERATE_AGGREGATE_EQUAL_ORDER(SSE2RVV_AGGREGATE_EQUAL_ORDER_)

#undef SSE2RVV_AGGREGATE_EQUAL_ORDER_IS_UBYTE
#undef SSE2RVV_AGGREGATE_EQUAL_ORDER_IS_UWORD

/* clang-format off */
#define SSE2RVV_GENERATE_CMP_EQUAL_ORDERED(prefix) \
    prefix##IMPL(byte)                              \
    prefix##IMPL(word)
/* clang-format on */

SSE2RVV_GENERATE_CMP_EQUAL_ORDERED(SSE2RVV_CMP_EQUAL_ORDERED_)

#define SSE2RVV_CMPESTR_LIST                           \
    _(CMP_UBYTE_EQUAL_ANY, cmp_byte_equal_any)         \
    _(CMP_UWORD_EQUAL_ANY, cmp_word_equal_any)         \
    _(CMP_SBYTE_EQUAL_ANY, cmp_byte_equal_any)         \
    _(CMP_SWORD_EQUAL_ANY, cmp_word_equal_any)         \
    _(CMP_UBYTE_RANGES, cmp_ubyte_ranges)              \
    _(CMP_UWORD_RANGES, cmp_uword_ranges)              \
    _(CMP_SBYTE_RANGES, cmp_sbyte_ranges)              \
    _(CMP_SWORD_RANGES, cmp_sword_ranges)              \
    _(CMP_UBYTE_EQUAL_EACH, cmp_byte_equal_each)       \
    _(CMP_UWORD_EQUAL_EACH, cmp_word_equal_each)       \
    _(CMP_SBYTE_EQUAL_EACH, cmp_byte_equal_each)       \
    _(CMP_SWORD_EQUAL_EACH, cmp_word_equal_each)       \
    _(CMP_UBYTE_EQUAL_ORDERED, cmp_byte_equal_ordered) \
    _(CMP_UWORD_EQUAL_ORDERED, cmp_word_equal_ordered) \
    _(CMP_SBYTE_EQUAL_ORDERED, cmp_byte_equal_ordered) \
    _(CMP_SWORD_EQUAL_ORDERED, cmp_word_equal_ordered)

enum {
#define _(name, func_suffix) name,
    SSE2RVV_CMPESTR_LIST
#undef _
};
typedef int (*cmpestr_func_t)(__m128i a, int la, __m128i b, int lb);
static cmpestr_func_t _sse2rvv_cmpfunc_table[] = {
#define _(name, func_suffix) _sse2rvv_##func_suffix,
    SSE2RVV_CMPESTR_LIST
#undef _
};

FORCE_INLINE int _sse2rvv_sido_negative(int res, int lb, int imm8, int bound)
{
    switch (imm8 & 0x30) {
    case _SIDD_NEGATIVE_POLARITY:
        res ^= 0xffffffff;
        break;
    case _SIDD_MASKED_NEGATIVE_POLARITY:
        res ^= (1 << lb) - 1;
        break;
    default:
        break;
    }

    return res & ((bound == 8) ? 0xFF : 0xFFFF);
}

FORCE_INLINE int _sse2rvv_clz(unsigned int x)
{
#ifdef _MSC_VER
    unsigned long cnt = 0;
    if (_BitScanReverse(&cnt, x))
        return 31 - cnt;
    return 32;
#else
    return x != 0 ? __builtin_clz(x) : 32;
#endif
}

FORCE_INLINE int _sse2rvv_ctz(unsigned int x)
{
#ifdef _MSC_VER
    unsigned long cnt = 0;
    if (_BitScanForward(&cnt, x))
        return cnt;
    return 32;
#else
    return x != 0 ? __builtin_ctz(x) : 32;
#endif
}

FORCE_INLINE int _sse2rvv_ctzll(unsigned long long x)
{
#ifdef _MSC_VER
    unsigned long cnt;
#if defined(SSE2RVV_HAS_BITSCAN64)
    if (_BitScanForward64(&cnt, x))
        return (int) (cnt);
#else
    if (_BitScanForward(&cnt, (unsigned long) (x)))
        return (int) cnt;
    if (_BitScanForward(&cnt, (unsigned long) (x >> 32)))
        return (int) (cnt + 32);
#endif /* SSE2RVV_HAS_BITSCAN64 */
    return 64;
#else /* assume GNU compatible compilers */
    return x != 0 ? __builtin_ctzll(x) : 64;
#endif
}

#define SSE2RVV_MIN(x, y) (x) < (y) ? (x) : (y)

#define SSE2RVV_CMPSTR_SET_UPPER(var, imm) const int var = (imm & 0x01) ? 8 : 16

#define SSE2RVV_CMPESTRX_LEN_PAIR(a, b, la, lb) \
    int tmp1 = la ^ (la >> 31);                 \
    la = tmp1 - (la >> 31);                     \
    int tmp2 = lb ^ (lb >> 31);                 \
    lb = tmp2 - (lb >> 31);                     \
    la = SSE2RVV_MIN(la, bound);                \
    lb = SSE2RVV_MIN(lb, bound)

// Compare all pairs of character in string a and b,
// then aggregate the result.
// As the only difference of PCMPESTR* and PCMPISTR* is the way to calculate the
// length of string, we use SSE2RVV_CMP{I,E}STRX_GET_LEN to get the length of
// string a and b.
#define SSE2RVV_COMP_AGG(a, b, la, lb, imm8, IE)                  \
    SSE2RVV_CMPSTR_SET_UPPER(bound, imm8);                        \
    SSE2RVV_##IE##_LEN_PAIR(a, b, la, lb);                        \
    int r2 = (_sse2rvv_cmpfunc_table[imm8 & 0x0f])(a, la, b, lb); \
    r2 = _sse2rvv_sido_negative(r2, lb, imm8, bound)

#define SSE2RVV_CMPSTR_GENERATE_INDEX(r2, bound, imm8) \
    return (r2 == 0)                                   \
               ? bound                                 \
               : ((imm8 & 0x40) ? (31 - _sse2rvv_clz(r2)) : _sse2rvv_ctz(r2))

#define SSE2RVV_CMPSTR_GENERATE_MASK(dst)                                      \
    __m128i dst = vreinterpretq_m128i_u8(vdupq_n_u8(0));                       \
    if (imm8 & 0x40) {                                                         \
        if (bound == 8) {                                                      \
            uint16x8_t tmp = vtstq_u16(vdupq_n_u16(r2),                        \
                                       vld1q_u16(_sse2rvv_cmpestr_mask16b));   \
            dst = vreinterpretq_m128i_u16(vbslq_u16(                           \
                tmp, vdupq_n_u16(-1), vreinterpretq_u16_m128i(dst)));          \
        } else {                                                               \
            uint8x16_t vec_r2 =                                                \
                vcombine_u8(vdup_n_u8(r2), vdup_n_u8(r2 >> 8));                \
            uint8x16_t tmp =                                                   \
                vtstq_u8(vec_r2, vld1q_u8(_sse2rvv_cmpestr_mask8b));           \
            dst = vreinterpretq_m128i_u8(                                      \
                vbslq_u8(tmp, vdupq_n_u8(-1), vreinterpretq_u8_m128i(dst)));   \
        }                                                                      \
    } else {                                                                   \
        if (bound == 16) {                                                     \
            dst = vreinterpretq_m128i_u16(                                     \
                vsetq_lane_u16(r2 & 0xffff, vreinterpretq_u16_m128i(dst), 0)); \
        } else {                                                               \
            dst = vreinterpretq_m128i_u8(                                      \
                vsetq_lane_u8(r2 & 0xff, vreinterpretq_u8_m128i(dst), 0));     \
        }                                                                      \
    }                                                                          \
    return dst

// Compare packed strings in a and b with lengths la and lb using the control
// in imm8, and returns 1 if b did not contain a null character and the
// resulting mask was zero, and 0 otherwise.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpestra
// FORCE_INLINE int _mm_cmpestra(__m128i a,
//                               int la,
//                               __m128i b,
//                               int lb,
//                               const int imm8) {}

// Compare packed strings in a and b with lengths la and lb using the control in
// imm8, and returns 1 if the resulting mask was non-zero, and 0 otherwise.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpestrc
// FORCE_INLINE int _mm_cmpestrc(__m128i a,
//                               int la,
//                               __m128i b,
//                               int lb,
//                               const int imm8) {}

// Compare packed strings in a and b with lengths la and lb using the control
// in imm8, and store the generated index in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpestri
// FORCE_INLINE int _mm_cmpestri(__m128i a,
//                               int la,
//                               __m128i b,
//                               int lb,
//                               const int imm8) {}

// Compare packed strings in a and b with lengths la and lb using the control
// in imm8, and store the generated mask in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpestrm
FORCE_INLINE __m128i
_mm_cmpestrm(__m128i a, int la, __m128i b, int lb, const int imm8)
{
    SSE2RVV_COMP_AGG(a, b, la, lb, imm8, CMPESTRX);
    SSE2RVV_CMPSTR_GENERATE_MASK(dst);
}

// Compare packed strings in a and b with lengths la and lb using the control in
// imm8, and returns bit 0 of the resulting bit mask.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpestro
// FORCE_INLINE int _mm_cmpestro(__m128i a,
//                               int la,
//                               __m128i b,
//                               int lb,
//                               const int imm8) {}

// Compare packed strings in a and b with lengths la and lb using the control in
// imm8, and returns 1 if any character in a was null, and 0 otherwise.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpestrs
// FORCE_INLINE int _mm_cmpestrs(__m128i a,
//                               int la,
//                               __m128i b,
//                               int lb,
//                               const int imm8) {}

// Compare packed strings in a and b with lengths la and lb using the control in
// imm8, and returns 1 if any character in b was null, and 0 otherwise.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpestrz
// FORCE_INLINE int _mm_cmpestrz(__m128i a,
//                               int la,
//                               __m128i b,
//                               int lb,
//                               const int imm8) {}

#define SSE2RVV_CMPISTRX_LENGTH(str, len, imm8)                          \
    do {                                                                 \
        if (imm8 & 0x01) {                                               \
            uint16x8_t equal_mask_##str =                                \
                vceqq_u16(vreinterpretq_u16_m128i(str), vdupq_n_u16(0)); \
            uint8x8_t res_##str = vshrn_n_u16(equal_mask_##str, 4);      \
            uint64_t matches_##str =                                     \
                vget_lane_u64(vreinterpret_u64_u8(res_##str), 0);        \
            len = _sse2rvv_ctzll(matches_##str) >> 3;                    \
        } else {                                                         \
            uint16x8_t equal_mask_##str = vreinterpretq_u16_u8(          \
                vceqq_u8(vreinterpretq_u8_m128i(str), vdupq_n_u8(0)));   \
            uint8x8_t res_##str = vshrn_n_u16(equal_mask_##str, 4);      \
            uint64_t matches_##str =                                     \
                vget_lane_u64(vreinterpret_u64_u8(res_##str), 0);        \
            len = _sse2rvv_ctzll(matches_##str) >> 2;                    \
        }                                                                \
    } while (0)

#define SSE2RVV_CMPISTRX_LEN_PAIR(a, b, la, lb) \
    int la, lb;                                 \
    do {                                        \
        SSE2RVV_CMPISTRX_LENGTH(a, la, imm8);   \
        SSE2RVV_CMPISTRX_LENGTH(b, lb, imm8);   \
    } while (0)

// Compare packed strings with implicit lengths in a and b using the control in
// imm8, and returns 1 if b did not contain a null character and the resulting
// mask was zero, and 0 otherwise.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpistra
// FORCE_INLINE int _mm_cmpistra(__m128i a, __m128i b, const int imm8) {}

// Compare packed strings with implicit lengths in a and b using the control in
// imm8, and returns 1 if the resulting mask was non-zero, and 0 otherwise.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpistrc
// FORCE_INLINE int _mm_cmpistrc(__m128i a, __m128i b, const int imm8) {}

// Compare packed strings with implicit lengths in a and b using the control in
// imm8, and store the generated index in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpistri
// FORCE_INLINE int _mm_cmpistri(__m128i a, __m128i b, const int imm8) {}

// Compare packed strings with implicit lengths in a and b using the control in
// imm8, and store the generated mask in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpistrm
// FORCE_INLINE __m128i _mm_cmpistrm(__m128i a, __m128i b, const int imm8) {}

// Compare packed strings with implicit lengths in a and b using the control in
// imm8, and returns bit 0 of the resulting bit mask.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpistro
// FORCE_INLINE int _mm_cmpistro(__m128i a, __m128i b, const int imm8) {}

// Compare packed strings with implicit lengths in a and b using the control in
// imm8, and returns 1 if any character in a was null, and 0 otherwise.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpistrs
// FORCE_INLINE int _mm_cmpistrs(__m128i a, __m128i b, const int imm8) {}

// Compare packed strings with implicit lengths in a and b using the control in
// imm8, and returns 1 if any character in b was null, and 0 otherwise.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_cmpistrz
// FORCE_INLINE int _mm_cmpistrz(__m128i a, __m128i b, const int imm8) {}

// Compares the 2 signed 64-bit integers in a and the 2 signed 64-bit integers
// in b for greater than.
// FORCE_INLINE __m128i _mm_cmpgt_epi64(__m128i a, __m128i b) {}

// Starting with the initial value in crc, accumulates a CRC32 value for
// unsigned 16-bit integer v, and stores the result in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_crc32_u16
// FORCE_INLINE uint32_t _mm_crc32_u16(uint32_t crc, uint16_t v) {}

// Starting with the initial value in crc, accumulates a CRC32 value for
// unsigned 32-bit integer v, and stores the result in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_crc32_u32
// FORCE_INLINE uint32_t _mm_crc32_u32(uint32_t crc, uint32_t v) {}

// Starting with the initial value in crc, accumulates a CRC32 value for
// unsigned 64-bit integer v, and stores the result in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_crc32_u64
// FORCE_INLINE uint64_t _mm_crc32_u64(uint64_t crc, uint64_t v) {}

// Starting with the initial value in crc, accumulates a CRC32 value for
// unsigned 8-bit integer v, and stores the result in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_crc32_u8
// FORCE_INLINE uint32_t _mm_crc32_u8(uint32_t crc, uint8_t v) {}

/* AES */

#if !defined(__ARM_FEATURE_CRYPTO) && (!defined(_M_ARM64) || defined(__clang__))
/* clang-format off */
#define SSE2RVV_AES_SBOX(w)                                           \
    {                                                                  \
        w(0x63), w(0x7c), w(0x77), w(0x7b), w(0xf2), w(0x6b), w(0x6f), \
        w(0xc5), w(0x30), w(0x01), w(0x67), w(0x2b), w(0xfe), w(0xd7), \
        w(0xab), w(0x76), w(0xca), w(0x82), w(0xc9), w(0x7d), w(0xfa), \
        w(0x59), w(0x47), w(0xf0), w(0xad), w(0xd4), w(0xa2), w(0xaf), \
        w(0x9c), w(0xa4), w(0x72), w(0xc0), w(0xb7), w(0xfd), w(0x93), \
        w(0x26), w(0x36), w(0x3f), w(0xf7), w(0xcc), w(0x34), w(0xa5), \
        w(0xe5), w(0xf1), w(0x71), w(0xd8), w(0x31), w(0x15), w(0x04), \
        w(0xc7), w(0x23), w(0xc3), w(0x18), w(0x96), w(0x05), w(0x9a), \
        w(0x07), w(0x12), w(0x80), w(0xe2), w(0xeb), w(0x27), w(0xb2), \
        w(0x75), w(0x09), w(0x83), w(0x2c), w(0x1a), w(0x1b), w(0x6e), \
        w(0x5a), w(0xa0), w(0x52), w(0x3b), w(0xd6), w(0xb3), w(0x29), \
        w(0xe3), w(0x2f), w(0x84), w(0x53), w(0xd1), w(0x00), w(0xed), \
        w(0x20), w(0xfc), w(0xb1), w(0x5b), w(0x6a), w(0xcb), w(0xbe), \
        w(0x39), w(0x4a), w(0x4c), w(0x58), w(0xcf), w(0xd0), w(0xef), \
        w(0xaa), w(0xfb), w(0x43), w(0x4d), w(0x33), w(0x85), w(0x45), \
        w(0xf9), w(0x02), w(0x7f), w(0x50), w(0x3c), w(0x9f), w(0xa8), \
        w(0x51), w(0xa3), w(0x40), w(0x8f), w(0x92), w(0x9d), w(0x38), \
        w(0xf5), w(0xbc), w(0xb6), w(0xda), w(0x21), w(0x10), w(0xff), \
        w(0xf3), w(0xd2), w(0xcd), w(0x0c), w(0x13), w(0xec), w(0x5f), \
        w(0x97), w(0x44), w(0x17), w(0xc4), w(0xa7), w(0x7e), w(0x3d), \
        w(0x64), w(0x5d), w(0x19), w(0x73), w(0x60), w(0x81), w(0x4f), \
        w(0xdc), w(0x22), w(0x2a), w(0x90), w(0x88), w(0x46), w(0xee), \
        w(0xb8), w(0x14), w(0xde), w(0x5e), w(0x0b), w(0xdb), w(0xe0), \
        w(0x32), w(0x3a), w(0x0a), w(0x49), w(0x06), w(0x24), w(0x5c), \
        w(0xc2), w(0xd3), w(0xac), w(0x62), w(0x91), w(0x95), w(0xe4), \
        w(0x79), w(0xe7), w(0xc8), w(0x37), w(0x6d), w(0x8d), w(0xd5), \
        w(0x4e), w(0xa9), w(0x6c), w(0x56), w(0xf4), w(0xea), w(0x65), \
        w(0x7a), w(0xae), w(0x08), w(0xba), w(0x78), w(0x25), w(0x2e), \
        w(0x1c), w(0xa6), w(0xb4), w(0xc6), w(0xe8), w(0xdd), w(0x74), \
        w(0x1f), w(0x4b), w(0xbd), w(0x8b), w(0x8a), w(0x70), w(0x3e), \
        w(0xb5), w(0x66), w(0x48), w(0x03), w(0xf6), w(0x0e), w(0x61), \
        w(0x35), w(0x57), w(0xb9), w(0x86), w(0xc1), w(0x1d), w(0x9e), \
        w(0xe1), w(0xf8), w(0x98), w(0x11), w(0x69), w(0xd9), w(0x8e), \
        w(0x94), w(0x9b), w(0x1e), w(0x87), w(0xe9), w(0xce), w(0x55), \
        w(0x28), w(0xdf), w(0x8c), w(0xa1), w(0x89), w(0x0d), w(0xbf), \
        w(0xe6), w(0x42), w(0x68), w(0x41), w(0x99), w(0x2d), w(0x0f), \
        w(0xb0), w(0x54), w(0xbb), w(0x16)                             \
    }
#define SSE2RVV_AES_RSBOX(w)                                          \
    {                                                                  \
        w(0x52), w(0x09), w(0x6a), w(0xd5), w(0x30), w(0x36), w(0xa5), \
        w(0x38), w(0xbf), w(0x40), w(0xa3), w(0x9e), w(0x81), w(0xf3), \
        w(0xd7), w(0xfb), w(0x7c), w(0xe3), w(0x39), w(0x82), w(0x9b), \
        w(0x2f), w(0xff), w(0x87), w(0x34), w(0x8e), w(0x43), w(0x44), \
        w(0xc4), w(0xde), w(0xe9), w(0xcb), w(0x54), w(0x7b), w(0x94), \
        w(0x32), w(0xa6), w(0xc2), w(0x23), w(0x3d), w(0xee), w(0x4c), \
        w(0x95), w(0x0b), w(0x42), w(0xfa), w(0xc3), w(0x4e), w(0x08), \
        w(0x2e), w(0xa1), w(0x66), w(0x28), w(0xd9), w(0x24), w(0xb2), \
        w(0x76), w(0x5b), w(0xa2), w(0x49), w(0x6d), w(0x8b), w(0xd1), \
        w(0x25), w(0x72), w(0xf8), w(0xf6), w(0x64), w(0x86), w(0x68), \
        w(0x98), w(0x16), w(0xd4), w(0xa4), w(0x5c), w(0xcc), w(0x5d), \
        w(0x65), w(0xb6), w(0x92), w(0x6c), w(0x70), w(0x48), w(0x50), \
        w(0xfd), w(0xed), w(0xb9), w(0xda), w(0x5e), w(0x15), w(0x46), \
        w(0x57), w(0xa7), w(0x8d), w(0x9d), w(0x84), w(0x90), w(0xd8), \
        w(0xab), w(0x00), w(0x8c), w(0xbc), w(0xd3), w(0x0a), w(0xf7), \
        w(0xe4), w(0x58), w(0x05), w(0xb8), w(0xb3), w(0x45), w(0x06), \
        w(0xd0), w(0x2c), w(0x1e), w(0x8f), w(0xca), w(0x3f), w(0x0f), \
        w(0x02), w(0xc1), w(0xaf), w(0xbd), w(0x03), w(0x01), w(0x13), \
        w(0x8a), w(0x6b), w(0x3a), w(0x91), w(0x11), w(0x41), w(0x4f), \
        w(0x67), w(0xdc), w(0xea), w(0x97), w(0xf2), w(0xcf), w(0xce), \
        w(0xf0), w(0xb4), w(0xe6), w(0x73), w(0x96), w(0xac), w(0x74), \
        w(0x22), w(0xe7), w(0xad), w(0x35), w(0x85), w(0xe2), w(0xf9), \
        w(0x37), w(0xe8), w(0x1c), w(0x75), w(0xdf), w(0x6e), w(0x47), \
        w(0xf1), w(0x1a), w(0x71), w(0x1d), w(0x29), w(0xc5), w(0x89), \
        w(0x6f), w(0xb7), w(0x62), w(0x0e), w(0xaa), w(0x18), w(0xbe), \
        w(0x1b), w(0xfc), w(0x56), w(0x3e), w(0x4b), w(0xc6), w(0xd2), \
        w(0x79), w(0x20), w(0x9a), w(0xdb), w(0xc0), w(0xfe), w(0x78), \
        w(0xcd), w(0x5a), w(0xf4), w(0x1f), w(0xdd), w(0xa8), w(0x33), \
        w(0x88), w(0x07), w(0xc7), w(0x31), w(0xb1), w(0x12), w(0x10), \
        w(0x59), w(0x27), w(0x80), w(0xec), w(0x5f), w(0x60), w(0x51), \
        w(0x7f), w(0xa9), w(0x19), w(0xb5), w(0x4a), w(0x0d), w(0x2d), \
        w(0xe5), w(0x7a), w(0x9f), w(0x93), w(0xc9), w(0x9c), w(0xef), \
        w(0xa0), w(0xe0), w(0x3b), w(0x4d), w(0xae), w(0x2a), w(0xf5), \
        w(0xb0), w(0xc8), w(0xeb), w(0xbb), w(0x3c), w(0x83), w(0x53), \
        w(0x99), w(0x61), w(0x17), w(0x2b), w(0x04), w(0x7e), w(0xba), \
        w(0x77), w(0xd6), w(0x26), w(0xe1), w(0x69), w(0x14), w(0x63), \
        w(0x55), w(0x21), w(0x0c), w(0x7d)                             \
    }
/* clang-format on */

/* X Macro trick. See https://en.wikipedia.org/wiki/X_Macro */
#define SSE2RVV_AES_H0(x) (x)
static const uint8_t _sse2rvv_sbox[256] = SSE2RVV_AES_SBOX(SSE2RVV_AES_H0);
static const uint8_t _sse2rvv_rsbox[256] = SSE2RVV_AES_RSBOX(SSE2RVV_AES_H0);
#undef SSE2RVV_AES_H0

/* x_time function and matrix multiply function */
#if !defined(__aarch64__) && !defined(_M_ARM64)
#define SSE2RVV_XT(x) (((x) << 1) ^ ((((x) >> 7) & 1) * 0x1b))
#define SSE2RVV_MULTIPLY(x, y)                                \
    (((y & 1) * x) ^ ((y >> 1 & 1) * SSE2RVV_XT(x)) ^         \
     ((y >> 2 & 1) * SSE2RVV_XT(SSE2RVV_XT(x))) ^             \
     ((y >> 3 & 1) * SSE2RVV_XT(SSE2RVV_XT(SSE2RVV_XT(x)))) ^ \
     ((y >> 4 & 1) * SSE2RVV_XT(SSE2RVV_XT(SSE2RVV_XT(SSE2RVV_XT(x))))))
#endif

// In the absence of crypto extensions, implement aesenc using regular NEON
// intrinsics instead. See:
// https://www.workofard.com/2017/01/accelerated-aes-for-the-arm64-linux-kernel/
// https://www.workofard.com/2017/07/ghash-for-low-end-cores/ and
// for more information.
// FORCE_INLINE __m128i _mm_aesenc_si128(__m128i a, __m128i RoundKey) {}

// Perform one round of an AES decryption flow on data (state) in a using the
// round key in RoundKey, and store the result in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_aesdec_si128
// FORCE_INLINE __m128i _mm_aesdec_si128(__m128i a, __m128i RoundKey) {}

// Perform the last round of an AES encryption flow on data (state) in a using
// the round key in RoundKey, and store the result in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_aesenclast_si128
// FORCE_INLINE __m128i _mm_aesenclast_si128(__m128i a, __m128i RoundKey) {}

// Perform the last round of an AES decryption flow on data (state) in a using
// the round key in RoundKey, and store the result in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_aesdeclast_si128
// FORCE_INLINE __m128i _mm_aesdeclast_si128(__m128i a, __m128i RoundKey) {}

// Perform the InvMixColumns transformation on a and store the result in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_aesimc_si128
// FORCE_INLINE __m128i _mm_aesimc_si128(__m128i a) {}

// Assist in expanding the AES cipher key by computing steps towards generating
// a round key for encryption cipher using data from a and an 8-bit round
// constant specified in imm8, and store the result in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_aeskeygenassist_si128
//
// Emits the Advanced Encryption Standard (AES) instruction aeskeygenassist.
// This instruction generates a round key for AES encryption. See
// https://kazakov.life/2017/11/01/cryptocurrency-mining-on-ios-devices/
// for details.
// FORCE_INLINE __m128i _mm_aeskeygenassist_si128(__m128i a, const int rcon) {}
#undef SSE2RVV_AES_SBOX
#undef SSE2RVV_AES_RSBOX

#if defined(__aarch64__)
#undef SSE2RVV_XT
#undef SSE2RVV_MULTIPLY
#endif

#else /* __ARM_FEATURE_CRYPTO */
// Implements equivalent of 'aesenc' by combining AESE (with an empty key) and
// AESMC and then manually applying the real key as an xor operation. This
// unfortunately means an additional xor op; the compiler should be able to
// optimize this away for repeated calls however. See
// https://blog.michaelbrase.com/2018/05/08/emulating-x86-aes-intrinsics-on-armv8-a
// for more details.
// FORCE_INLINE __m128i _mm_aesenc_si128(__m128i a, __m128i b) {}

// Perform one round of an AES decryption flow on data (state) in a using the
// round key in RoundKey, and store the result in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_aesdec_si128
// FORCE_INLINE __m128i _mm_aesdec_si128(__m128i a, __m128i RoundKey) {}

// Perform the last round of an AES encryption flow on data (state) in a using
// the round key in RoundKey, and store the result in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_aesenclast_si128
// FORCE_INLINE __m128i _mm_aesenclast_si128(__m128i a, __m128i RoundKey) {}

// Perform the last round of an AES decryption flow on data (state) in a using
// the round key in RoundKey, and store the result in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_aesdeclast_si128
// FORCE_INLINE __m128i _mm_aesdeclast_si128(__m128i a, __m128i RoundKey) {}

// Perform the InvMixColumns transformation on a and store the result in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_aesimc_si128
// FORCE_INLINE __m128i _mm_aesimc_si128(__m128i a) {}

// Assist in expanding the AES cipher key by computing steps towards generating
// a round key for encryption cipher using data from a and an 8-bit round
// constant specified in imm8, and store the result in dst."
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_aeskeygenassist_si128
// FORCE_INLINE __m128i _mm_aeskeygenassist_si128(__m128i a, const int rcon) {}
#endif

/* Others */

// Perform a carry-less multiplication of two 64-bit integers, selected from a
// and b according to imm8, and store the results in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_clmulepi64_si128
// FORCE_INLINE __m128i _mm_clmulepi64_si128(__m128i _a, __m128i _b, const int imm) {}

// FORCE_INLINE unsigned int _sse2rvv_mm_get_denormals_zero_mode(void) {}

// Count the number of bits set to 1 in unsigned 32-bit integer a, and
// return that count in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_popcnt_u32
// FORCE_INLINE int _mm_popcnt_u32(unsigned int a) {}

// Count the number of bits set to 1 in unsigned 64-bit integer a, and
// return that count in dst.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_popcnt_u64
// FORCE_INLINE int64_t _mm_popcnt_u64(uint64_t a) {}

// FORCE_INLINE void _sse2rvv_mm_set_denormals_zero_mode(unsigned int flag) {}

// Return the current 64-bit value of the processor's time-stamp counter.
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=rdtsc
FORCE_INLINE uint64_t _rdtsc(void)
{
#if defined(__aarch64__) || defined(_M_ARM64)
    uint64_t val;

    /* According to ARM DDI 0487F.c, from Armv8.0 to Armv8.5 inclusive, the
     * system counter is at least 56 bits wide; from Armv8.6, the counter
     * must be 64 bits wide.  So the system counter could be less than 64
     * bits wide and it is attributed with the flag 'cap_user_time_short'
     * is true.
     */
#if defined(_MSC_VER)
    val = _ReadStatusReg(ARM64_SYSREG(3, 3, 14, 0, 2));
#else
    __asm__ __volatile__("mrs %0, cntvct_el0" : "=r"(val));
#endif

    return val;
#else
    uint32_t pmccntr, pmuseren, pmcntenset;
    // Read the user mode Performance Monitoring Unit (PMU)
    // User Enable Register (PMUSERENR) access permissions.
    __asm__ __volatile__("mrc p15, 0, %0, c9, c14, 0" : "=r"(pmuseren));
    if (pmuseren & 1) {  // Allows reading PMUSERENR for user mode code.
        __asm__ __volatile__("mrc p15, 0, %0, c9, c12, 1" : "=r"(pmcntenset));
        if (pmcntenset & 0x80000000UL) {  // Is it counting?
            __asm__ __volatile__("mrc p15, 0, %0, c9, c13, 0" : "=r"(pmccntr));
            // The counter is set up to count every 64th cycle
            return (uint64_t) (pmccntr) << 6;
        }
    }

    // Fallback to syscall as we can't enable PMUSERENR in user mode.
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (uint64_t) (tv.tv_sec) * 1000000 + tv.tv_usec;
#endif
}

#if defined(__GNUC__) || defined(__clang__)
#pragma pop_macro("ALIGN_STRUCT")
#pragma pop_macro("FORCE_INLINE")
#endif

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC pop_options
#endif

#endif