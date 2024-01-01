#include <assert.h>
#include <float.h>
#include <inttypes.h>
#include <math.h>
#include <stdalign.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <utility>

#include "binding.h"
#include "impl.h"

// Try 10,000 random floating point values for each test we run
#define MAX_TEST_VALUE 10000

/* Pattern Matching for C macros.
 * https://github.com/pfultz2/Cloak/wiki/C-Preprocessor-tricks,-tips,-and-idioms
 */

/* catenate */
#define PRIMITIVE_CAT(a, ...) a##__VA_ARGS__

#define IIF(c) PRIMITIVE_CAT(IIF_, c)
/* run the 2nd parameter */
#define IIF_0(t, ...) __VA_ARGS__
/* run the 1st parameter */
#define IIF_1(t, ...) t

// This program a set of unit tests to ensure that each SSE call provide the
// output we expect.  If this fires an assert, then something didn't match up.
//
// Functions with "test_" prefix will be called in run_single_test.
namespace SSE2RVV {
// Forward declaration
class SSE2RVV_TEST_IMPL : public SSE2RVV_TEST {
public:
  SSE2RVV_TEST_IMPL(void);
  result_t load_test_float_pointers(uint32_t i);
  result_t load_test_int_pointers(uint32_t i);
  result_t run_single_test(INSTRUCTION_TEST test, uint32_t i);

  float *test_cases_float_pointer1;
  float *test_cases_float_pointer2;
  int32_t *test_cases_int_pointer1;
  int32_t *test_cases_int_pointer2;
  float test_cases_floats[MAX_TEST_VALUE];
  int32_t test_cases_ints[MAX_TEST_VALUE];

  virtual ~SSE2RVV_TEST_IMPL(void) {
    platform_aligned_free(test_cases_float_pointer1);
    platform_aligned_free(test_cases_float_pointer2);
    platform_aligned_free(test_cases_int_pointer1);
    platform_aligned_free(test_cases_int_pointer2);
  }
  virtual void release(void) { delete this; }
  virtual result_t run_test(INSTRUCTION_TEST test) {
    result_t ret = TEST_SUCCESS;

    // Test a whole bunch of values
    for (uint32_t i = 0; i < (MAX_TEST_VALUE - 8); i++) {
      ret = load_test_float_pointers(i); // Load some random float values
      if (ret == TEST_FAIL)
        break;                         // load test float failed??
      ret = load_test_int_pointers(i); // load some random int values
      if (ret == TEST_FAIL)
        break; // load test float failed??
      // If we are testing the reciprocal, then invert the input data
      // (easier for debugging)
      if (test == it_mm_rcp_ps) {
        test_cases_float_pointer1[0] = 1.0f / test_cases_float_pointer1[0];
        test_cases_float_pointer1[1] = 1.0f / test_cases_float_pointer1[1];
        test_cases_float_pointer1[2] = 1.0f / test_cases_float_pointer1[2];
        test_cases_float_pointer1[3] = 1.0f / test_cases_float_pointer1[3];
      }
      if (test == it_mm_rcp_ps || test == it_mm_rcp_ss ||
          test == it_mm_rsqrt_ps || test == it_mm_rsqrt_ss) {
        if ((rand() & 3) == 0) {
          uint32_t r1 = rand() & 3;
          uint32_t r2 = rand() & 3;
          uint32_t r3 = rand() & 3;
          uint32_t r4 = rand() & 3;
          uint32_t r5 = rand() & 3;
          uint32_t r6 = rand() & 3;
          uint32_t r7 = rand() & 3;
          uint32_t r8 = rand() & 3;
          test_cases_float_pointer1[r1] = 0.0f;
          test_cases_float_pointer1[r2] = 0.0f;
          test_cases_float_pointer1[r3] = 0.0f;
          test_cases_float_pointer1[r4] = 0.0f;
          test_cases_float_pointer1[r5] = -0.0f;
          test_cases_float_pointer1[r6] = -0.0f;
          test_cases_float_pointer1[r7] = -0.0f;
          test_cases_float_pointer1[r8] = -0.0f;
        }
      }
      if (test == it_mm_cmpge_ps || test == it_mm_cmpge_ss ||
          test == it_mm_cmple_ps || test == it_mm_cmple_ss ||
          test == it_mm_cmpeq_ps || test == it_mm_cmpeq_ss) {
        // Make sure at least one value is the same.
        test_cases_float_pointer1[3] = test_cases_float_pointer2[3];
      }

      if (test == it_mm_cmpord_ps || test == it_mm_cmpord_ss ||
          test == it_mm_cmpunord_ps || test == it_mm_cmpunord_ss ||
          test == it_mm_cmpeq_ps || test == it_mm_cmpeq_ss ||
          test == it_mm_cmpge_ps || test == it_mm_cmpge_ss ||
          test == it_mm_cmpgt_ps || test == it_mm_cmpgt_ss ||
          test == it_mm_cmple_ps || test == it_mm_cmple_ss ||
          test == it_mm_cmplt_ps || test == it_mm_cmplt_ss ||
          test == it_mm_cmpneq_ps || test == it_mm_cmpneq_ss ||
          test == it_mm_cmpnge_ps || test == it_mm_cmpnge_ss ||
          test == it_mm_cmpngt_ps || test == it_mm_cmpngt_ss ||
          test == it_mm_cmpnle_ps || test == it_mm_cmpnle_ss ||
          test == it_mm_cmpnlt_ps || test == it_mm_cmpnlt_ss ||
          test == it_mm_comieq_ss || test == it_mm_ucomieq_ss ||
          test == it_mm_comige_ss || test == it_mm_ucomige_ss ||
          test == it_mm_comigt_ss || test == it_mm_ucomigt_ss ||
          test == it_mm_comile_ss || test == it_mm_ucomile_ss ||
          test == it_mm_comilt_ss || test == it_mm_ucomilt_ss ||
          test == it_mm_comineq_ss || test == it_mm_ucomineq_ss) {
        // Make sure the NaN values are included in the testing
        // one out of four times.
        if ((rand() & 3) == 0) {
          uint32_t r1 = rand() & 3;
          uint32_t r2 = rand() & 3;
          test_cases_float_pointer1[r1] = nanf("");
          test_cases_float_pointer2[r2] = nanf("");
        }
      }

      if (test == it_mm_cmpord_pd || test == it_mm_cmpord_sd ||
          test == it_mm_cmpunord_pd || test == it_mm_cmpunord_sd ||
          test == it_mm_cmpeq_pd || test == it_mm_cmpeq_sd ||
          test == it_mm_cmpge_pd || test == it_mm_cmpge_sd ||
          test == it_mm_cmpgt_pd || test == it_mm_cmpgt_sd ||
          test == it_mm_cmple_pd || test == it_mm_cmple_sd ||
          test == it_mm_cmplt_pd || test == it_mm_cmplt_sd ||
          test == it_mm_cmpneq_pd || test == it_mm_cmpneq_sd ||
          test == it_mm_cmpnge_pd || test == it_mm_cmpnge_sd ||
          test == it_mm_cmpngt_pd || test == it_mm_cmpngt_sd ||
          test == it_mm_cmpnle_pd || test == it_mm_cmpnle_sd ||
          test == it_mm_cmpnlt_pd || test == it_mm_cmpnlt_sd ||
          test == it_mm_comieq_sd || test == it_mm_ucomieq_sd ||
          test == it_mm_comige_sd || test == it_mm_ucomige_sd ||
          test == it_mm_comigt_sd || test == it_mm_ucomigt_sd ||
          test == it_mm_comile_sd || test == it_mm_ucomile_sd ||
          test == it_mm_comilt_sd || test == it_mm_ucomilt_sd ||
          test == it_mm_comineq_sd || test == it_mm_ucomineq_sd) {
        // Make sure the NaN values are included in the testing
        // one out of four times.
        if ((rand() & 3) == 0) {
          // FIXME:
          // The argument "0xFFFFFFFFFFFF" is a tricky workaround to
          // set the NaN value for doubles. The code is not intuitive
          // and should be fixed in the future.
          uint32_t r1 = ((rand() & 1) << 1) + 1;
          uint32_t r2 = ((rand() & 1) << 1) + 1;
          test_cases_float_pointer1[r1] = nanf("0xFFFFFFFFFFFF");
          test_cases_float_pointer2[r2] = nanf("0xFFFFFFFFFFFF");
        }
      }

      if (test == it_mm_max_pd || test == it_mm_max_sd ||
          test == it_mm_min_pd || test == it_mm_min_sd) {
        // Make sure the positive/negative infinity values are included
        // in the testing one out of four times.
        if ((rand() & 3) == 0) {
          uint32_t r1 = ((rand() & 1) << 1) + 1;
          uint32_t r2 = ((rand() & 1) << 1) + 1;
          uint32_t r3 = ((rand() & 1) << 1) + 1;
          uint32_t r4 = ((rand() & 1) << 1) + 1;
          test_cases_float_pointer1[r1] = INFINITY;
          test_cases_float_pointer2[r2] = INFINITY;
          test_cases_float_pointer1[r3] = -INFINITY;
          test_cases_float_pointer1[r4] = -INFINITY;
        }
      }

#if SSE2RVV_PRECISE_MINMAX
      if (test == it_mm_max_ps || test == it_mm_max_ss ||
          test == it_mm_min_ps || test == it_mm_min_ss) {
        // Make sure the NaN values are included in the testing
        // one out of four times.
        if ((rand() & 3) == 0) {
          uint32_t r1 = rand() & 3;
          uint32_t r2 = rand() & 3;
          test_cases_float_pointer1[r1] = nanf("");
          test_cases_float_pointer2[r2] = nanf("");
        }
      }

      if (test == it_mm_max_pd || test == it_mm_max_sd ||
          test == it_mm_min_pd || test == it_mm_min_sd) {
        // Make sure the NaN values are included in the testing
        // one out of four times.
        if ((rand() & 3) == 0) {
          // FIXME:
          // The argument "0xFFFFFFFFFFFF" is a tricky workaround to
          // set the NaN value for doubles. The code is not intuitive
          // and should be fixed in the future.
          uint32_t r1 = ((rand() & 1) << 1) + 1;
          uint32_t r2 = ((rand() & 1) << 1) + 1;
          test_cases_float_pointer1[r1] = nanf("0xFFFFFFFFFFFF");
          test_cases_float_pointer2[r2] = nanf("0xFFFFFFFFFFFF");
        }
      }
#endif

      // one out of every random 64 times or so, mix up the test floats to
      // contain some integer values
      if ((rand() & 63) == 0) {
        uint32_t option = rand() & 3;
        switch (option) {
        // All integers..
        case 0:
          test_cases_float_pointer1[0] = float(test_cases_int_pointer1[0]);
          test_cases_float_pointer1[1] = float(test_cases_int_pointer1[1]);
          test_cases_float_pointer1[2] = float(test_cases_int_pointer1[2]);
          test_cases_float_pointer1[3] = float(test_cases_int_pointer1[3]);

          test_cases_float_pointer2[0] = float(test_cases_int_pointer2[0]);
          test_cases_float_pointer2[1] = float(test_cases_int_pointer2[1]);
          test_cases_float_pointer2[2] = float(test_cases_int_pointer2[2]);
          test_cases_float_pointer2[3] = float(test_cases_int_pointer2[3]);

          break;
        case 1: {
          uint32_t index = rand() & 3;
          test_cases_float_pointer1[index] =
              float(test_cases_int_pointer1[index]);
          index = rand() & 3;
          test_cases_float_pointer2[index] =
              float(test_cases_int_pointer2[index]);
        } break;
        case 2: {
          uint32_t index1 = rand() & 3;
          uint32_t index2 = rand() & 3;
          test_cases_float_pointer1[index1] =
              float(test_cases_int_pointer1[index1]);
          test_cases_float_pointer1[index2] =
              float(test_cases_int_pointer1[index2]);
          index1 = rand() & 3;
          index2 = rand() & 3;
          test_cases_float_pointer1[index1] =
              float(test_cases_int_pointer1[index1]);
          test_cases_float_pointer1[index2] =
              float(test_cases_int_pointer1[index2]);
        } break;
        case 3:
          test_cases_float_pointer1[0] = float(test_cases_int_pointer1[0]);
          test_cases_float_pointer1[1] = float(test_cases_int_pointer1[1]);
          test_cases_float_pointer1[2] = float(test_cases_int_pointer1[2]);
          test_cases_float_pointer1[3] = float(test_cases_int_pointer1[3]);
          break;
        }
        if ((rand() & 3) == 0) { // one out of 4 times, make halves
          for (uint32_t j = 0; j < 4; j++) {
            test_cases_float_pointer1[j] *= 0.5f;
            test_cases_float_pointer2[j] *= 0.5f;
          }
        }
      }

      ret = run_single_test(test, i);
      if (ret == TEST_FAIL) // the test failed...
      {
        // Set a breakpoint here if you want to step through the failure
        // case in the debugger
        ret = run_single_test(test, i);
        break;
      }
    }
    return ret;
  }
};

const char *instructionString[] = {
#define _(x) #x,
    INTRIN_LIST
#undef _
};

// Produce rounding which is the same as SSE instructions with _MM_ROUND_NEAREST
// rounding mode
static inline float bankers_rounding(float val) {
  if (val < 0)
    return -bankers_rounding(-val);

  float ret;
  float roundDown = floorf(val); // Round down value
  float roundUp = ceilf(val);    // Round up value
  float diffDown = val - roundDown;
  float diffUp = roundUp - val;

  if (diffDown < diffUp) {
    /* If it's closer to the round down value, then use it */
    ret = roundDown;
  } else if (diffDown > diffUp) {
    /* If it's closer to the round up value, then use it */
    ret = roundUp;
  } else {
    /* If it's equidistant between round up and round down value, pick the
     * one which is an even number */
    float half = roundDown / 2;
    if (half != floorf(half)) {
      /* If the round down value is odd, return the round up value */
      ret = roundUp;
    } else {
      /* If the round up value is odd, return the round down value */
      ret = roundDown;
    }
  }
  return ret;
}

static inline double bankers_rounding(double val) {
  if (val < 0)
    return -bankers_rounding(-val);

  double ret;
  double roundDown = floor(val); // Round down value
  double roundUp = ceil(val);    // Round up value
  double diffDown = val - roundDown;
  double diffUp = roundUp - val;

  if (diffDown < diffUp) {
    /* If it's closer to the round down value, then use it */
    ret = roundDown;
  } else if (diffDown > diffUp) {
    /* If it's closer to the round up value, then use it */
    ret = roundUp;
  } else {
    /* If it's equidistant between round up and round down value, pick the
     * one which is an even number */
    double half = roundDown / 2;
    if (half != floor(half)) {
      /* If the round down value is odd, return the round up value */
      ret = roundUp;
    } else {
      /* If the round up value is odd, return the round down value */
      ret = roundDown;
    }
  }
  return ret;
}

// SplitMix64 PRNG by Sebastiano Vigna, see:
// <https://xoshiro.di.unimi.it/splitmix64.c>
static uint64_t state; // the state of SplitMix64 PRNG
const double TWOPOWER64 = pow(2, 64);

#define SSE2RVV_INIT_RNG(seed)                                                 \
  do {                                                                         \
    state = seed;                                                              \
  } while (0)

static double next() {
  uint64_t z = (state += 0x9e3779b97f4a7c15);
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
  z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
  return z ^ (z >> 31);
}

static float ranf() { return next() / TWOPOWER64; }

static float ranf(float low, float high) { return ranf() * (high - low) + low; }

// Enable the tests which are using the macro of another tests
result_t test_mm_slli_si128(const SSE2RVV_TEST_IMPL &impl, uint32_t iter);
result_t test_mm_srli_si128(const SSE2RVV_TEST_IMPL &impl, uint32_t iter);
result_t test_mm_shuffle_pi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter);

// This function is not called from "run_single_test", but for other intrinsic
// tests that might need to call "_mm_set_epi32".
__m128i do_mm_set_epi32(int32_t x, int32_t y, int32_t z, int32_t w) {
  __m128i a = _mm_set_epi32(x, y, z, w);
  validate_int32(a, w, z, y, x);
  return a;
}

// This function is not called from "run_single_test", but for other intrinsic
// tests that might need to load __m64 data.
template <class T> __m64 load_m64(const T *p) { return *((const __m64 *)p); }

// This function is not called from "run_single_test", but for other intrinsic
// tests that might need to call "_mm_load_ps".
template <class T> __m128 load_m128(const T *p) {
  return _mm_loadu_ps((const float *)p);
}

// This function is not called from "run_single_test", but for other intrinsic
// tests that might need to call "_mm_load_ps".
template <class T> __m128i load_m128i(const T *p) {
  __m128 a = _mm_loadu_ps((const float *)p);
  __m128i ia = *(const __m128i *)&a;
  return ia;
}

// This function is not called from "run_single_test", but for other intrinsic
// tests that might need to call "_mm_load_pd".
template <class T> __m128d load_m128d(const T *p) {
  return _mm_loadu_pd((const double *)p);
}

// This function is not called from "run_single_test", but for other intrinsic
// tests that might need to call "_mm_store_ps".
result_t do_mm_store_ps(float *p, float x, float y, float z, float w) {
  __m128 a = _mm_set_ps(x, y, z, w);
  _mm_store_ps(p, a);
  ASSERT_RETURN(p[0] == w);
  ASSERT_RETURN(p[1] == z);
  ASSERT_RETURN(p[2] == y);
  ASSERT_RETURN(p[3] == x);
  return TEST_SUCCESS;
}

// This function is not called from "run_single_test", but for other intrinsic
// tests that might need to call "_mm_store_ps".
result_t do_mm_store_ps(int32_t *p, int32_t x, int32_t y, int32_t z,
                        int32_t w) {
  __m128i a = _mm_set_epi32(x, y, z, w);
  _mm_store_ps((float *)p, *(const __m128 *)&a);
  ASSERT_RETURN(p[0] == w);
  ASSERT_RETURN(p[1] == z);
  ASSERT_RETURN(p[2] == y);
  ASSERT_RETURN(p[3] == x);
  return TEST_SUCCESS;
}

float cmp_noNaN(float a, float b) {
  return (!isnan(a) && !isnan(b)) ? ALL_BIT_1_32 : 0.0f;
}

double cmp_noNaN(double a, double b) {
  return (!isnan(a) && !isnan(b)) ? ALL_BIT_1_64 : 0.0f;
}

float cmp_hasNaN(float a, float b) {
  return (isnan(a) || isnan(b)) ? ALL_BIT_1_32 : 0.0f;
}

double cmp_hasNaN(double a, double b) {
  return (isnan(a) || isnan(b)) ? ALL_BIT_1_64 : 0.0f;
}

int32_t comilt_ss(float a, float b) {
  if (isnan(a) || isnan(b))
    return 0;
  return (a < b);
}

int32_t comigt_ss(float a, float b) {
  if (isnan(a) || isnan(b))
    return 0;
  return (a > b);
}

int32_t comile_ss(float a, float b) {
  if (isnan(a) || isnan(b))
    return 0;
  return (a <= b);
}

int32_t comige_ss(float a, float b) {
  if (isnan(a) || isnan(b))
    return 0;
  return (a >= b);
}

int32_t comieq_ss(float a, float b) {
  if (isnan(a) || isnan(b))
    return 0;
  return (a == b);
}

int32_t comineq_ss(float a, float b) {
  if (isnan(a) || isnan(b))
    return 1;
  return (a != b);
}

static inline int16_t saturate_16(int32_t a) {
  int32_t max = (1 << 15) - 1;
  int32_t min = -(1 << 15);
  if (a > max)
    return max;
  if (a < min)
    return min;
  return a;
}

uint32_t canonical_crc32_u8(uint32_t crc, uint8_t v) {
  crc ^= v;
  for (int bit = 0; bit < 8; bit++) {
    if (crc & 1)
      crc = (crc >> 1) ^ uint32_t(0x82f63b78);
    else
      crc = (crc >> 1);
  }
  return crc;
}

uint32_t canonical_crc32_u16(uint32_t crc, uint16_t v) {
  crc = canonical_crc32_u8(crc, v & 0xff);
  crc = canonical_crc32_u8(crc, (v >> 8) & 0xff);
  return crc;
}

uint32_t canonical_crc32_u32(uint32_t crc, uint32_t v) {
  crc = canonical_crc32_u16(crc, v & 0xffff);
  crc = canonical_crc32_u16(crc, (v >> 16) & 0xffff);
  return crc;
}

uint64_t canonical_crc32_u64(uint64_t crc, uint64_t v) {
  crc = canonical_crc32_u32((uint32_t)(crc), v & 0xffffffff);
  crc = canonical_crc32_u32((uint32_t)(crc), (v >> 32) & 0xffffffff);
  return crc;
}

static const uint8_t crypto_aes_sbox[256] = {
    0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5, 0x30, 0x01, 0x67, 0x2b,
    0xfe, 0xd7, 0xab, 0x76, 0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0,
    0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0, 0xb7, 0xfd, 0x93, 0x26,
    0x36, 0x3f, 0xf7, 0xcc, 0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15,
    0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a, 0x07, 0x12, 0x80, 0xe2,
    0xeb, 0x27, 0xb2, 0x75, 0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0,
    0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84, 0x53, 0xd1, 0x00, 0xed,
    0x20, 0xfc, 0xb1, 0x5b, 0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf,
    0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85, 0x45, 0xf9, 0x02, 0x7f,
    0x50, 0x3c, 0x9f, 0xa8, 0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5,
    0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2, 0xcd, 0x0c, 0x13, 0xec,
    0x5f, 0x97, 0x44, 0x17, 0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73,
    0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88, 0x46, 0xee, 0xb8, 0x14,
    0xde, 0x5e, 0x0b, 0xdb, 0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c,
    0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79, 0xe7, 0xc8, 0x37, 0x6d,
    0x8d, 0xd5, 0x4e, 0xa9, 0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08,
    0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6, 0xe8, 0xdd, 0x74, 0x1f,
    0x4b, 0xbd, 0x8b, 0x8a, 0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e,
    0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e, 0xe1, 0xf8, 0x98, 0x11,
    0x69, 0xd9, 0x8e, 0x94, 0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf,
    0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68, 0x41, 0x99, 0x2d, 0x0f,
    0xb0, 0x54, 0xbb, 0x16,
};

static const uint8_t crypto_aes_rsbox[256] = {
    0x52, 0x09, 0x6a, 0xd5, 0x30, 0x36, 0xa5, 0x38, 0xbf, 0x40, 0xa3, 0x9e,
    0x81, 0xf3, 0xd7, 0xfb, 0x7c, 0xe3, 0x39, 0x82, 0x9b, 0x2f, 0xff, 0x87,
    0x34, 0x8e, 0x43, 0x44, 0xc4, 0xde, 0xe9, 0xcb, 0x54, 0x7b, 0x94, 0x32,
    0xa6, 0xc2, 0x23, 0x3d, 0xee, 0x4c, 0x95, 0x0b, 0x42, 0xfa, 0xc3, 0x4e,
    0x08, 0x2e, 0xa1, 0x66, 0x28, 0xd9, 0x24, 0xb2, 0x76, 0x5b, 0xa2, 0x49,
    0x6d, 0x8b, 0xd1, 0x25, 0x72, 0xf8, 0xf6, 0x64, 0x86, 0x68, 0x98, 0x16,
    0xd4, 0xa4, 0x5c, 0xcc, 0x5d, 0x65, 0xb6, 0x92, 0x6c, 0x70, 0x48, 0x50,
    0xfd, 0xed, 0xb9, 0xda, 0x5e, 0x15, 0x46, 0x57, 0xa7, 0x8d, 0x9d, 0x84,
    0x90, 0xd8, 0xab, 0x00, 0x8c, 0xbc, 0xd3, 0x0a, 0xf7, 0xe4, 0x58, 0x05,
    0xb8, 0xb3, 0x45, 0x06, 0xd0, 0x2c, 0x1e, 0x8f, 0xca, 0x3f, 0x0f, 0x02,
    0xc1, 0xaf, 0xbd, 0x03, 0x01, 0x13, 0x8a, 0x6b, 0x3a, 0x91, 0x11, 0x41,
    0x4f, 0x67, 0xdc, 0xea, 0x97, 0xf2, 0xcf, 0xce, 0xf0, 0xb4, 0xe6, 0x73,
    0x96, 0xac, 0x74, 0x22, 0xe7, 0xad, 0x35, 0x85, 0xe2, 0xf9, 0x37, 0xe8,
    0x1c, 0x75, 0xdf, 0x6e, 0x47, 0xf1, 0x1a, 0x71, 0x1d, 0x29, 0xc5, 0x89,
    0x6f, 0xb7, 0x62, 0x0e, 0xaa, 0x18, 0xbe, 0x1b, 0xfc, 0x56, 0x3e, 0x4b,
    0xc6, 0xd2, 0x79, 0x20, 0x9a, 0xdb, 0xc0, 0xfe, 0x78, 0xcd, 0x5a, 0xf4,
    0x1f, 0xdd, 0xa8, 0x33, 0x88, 0x07, 0xc7, 0x31, 0xb1, 0x12, 0x10, 0x59,
    0x27, 0x80, 0xec, 0x5f, 0x60, 0x51, 0x7f, 0xa9, 0x19, 0xb5, 0x4a, 0x0d,
    0x2d, 0xe5, 0x7a, 0x9f, 0x93, 0xc9, 0x9c, 0xef, 0xa0, 0xe0, 0x3b, 0x4d,
    0xae, 0x2a, 0xf5, 0xb0, 0xc8, 0xeb, 0xbb, 0x3c, 0x83, 0x53, 0x99, 0x61,
    0x17, 0x2b, 0x04, 0x7e, 0xba, 0x77, 0xd6, 0x26, 0xe1, 0x69, 0x14, 0x63,
    0x55, 0x21, 0x0c, 0x7d,
};

// XT is x_time function that muliplies 'x' by 2 in GF(2^8)
#define XT(x) (((x) << 1) ^ ((((x) >> 7) & 1) * 0x1b))

inline __m128i aesenc_128_reference(__m128i a, __m128i b) {
  uint8_t i, t, u, v[4][4];
  for (i = 0; i < 16; ++i) {
    v[((i / 4) + 4 - (i % 4)) % 4][i % 4] =
        crypto_aes_sbox[((SIMDVec *)&a)->m128_u8[i]];
  }
  for (i = 0; i < 4; ++i) {
    t = v[i][0];
    u = v[i][0] ^ v[i][1] ^ v[i][2] ^ v[i][3];
    v[i][0] ^= u ^ XT(v[i][0] ^ v[i][1]);
    v[i][1] ^= u ^ XT(v[i][1] ^ v[i][2]);
    v[i][2] ^= u ^ XT(v[i][2] ^ v[i][3]);
    v[i][3] ^= u ^ XT(v[i][3] ^ t);
  }

  for (i = 0; i < 16; ++i) {
    ((SIMDVec *)&a)->m128_u8[i] = v[i / 4][i % 4] ^ ((SIMDVec *)&b)->m128_u8[i];
  }

  return a;
}

#define MULTIPLY(x, y)                                                         \
  (((y & 1) * x) ^ ((y >> 1 & 1) * XT(x)) ^ ((y >> 2 & 1) * XT(XT(x))) ^       \
   ((y >> 3 & 1) * XT(XT(XT(x)))) ^ ((y >> 4 & 1) * XT(XT(XT(XT(x))))))

inline __m128i aesdec_128_reference(__m128i a, __m128i b) {
  uint8_t i, e, f, g, h, v[4][4];
  for (i = 0; i < 16; ++i) {
    v[((i / 4) + (i % 4)) % 4][i % 4] =
        crypto_aes_rsbox[((SIMDVec *)&a)->m128_u8[i]];
  }

  for (i = 0; i < 4; ++i) {
    e = v[i][0];
    f = v[i][1];
    g = v[i][2];
    h = v[i][3];

    v[i][0] = MULTIPLY(e, 0x0e) ^ MULTIPLY(f, 0x0b) ^ MULTIPLY(g, 0x0d) ^
              MULTIPLY(h, 0x09);
    v[i][1] = MULTIPLY(e, 0x09) ^ MULTIPLY(f, 0x0e) ^ MULTIPLY(g, 0x0b) ^
              MULTIPLY(h, 0x0d);
    v[i][2] = MULTIPLY(e, 0x0d) ^ MULTIPLY(f, 0x09) ^ MULTIPLY(g, 0x0e) ^
              MULTIPLY(h, 0x0b);
    v[i][3] = MULTIPLY(e, 0x0b) ^ MULTIPLY(f, 0x0d) ^ MULTIPLY(g, 0x09) ^
              MULTIPLY(h, 0x0e);
  }

  for (i = 0; i < 16; ++i) {
    ((SIMDVec *)&a)->m128_u8[i] = v[i / 4][i % 4] ^ ((SIMDVec *)&b)->m128_u8[i];
  }
  return a;
}

inline __m128i aesenclast_128_reference(__m128i s, __m128i rk) {
  uint8_t i, v[4][4];
  for (i = 0; i < 16; ++i)
    v[((i / 4) + 4 - (i % 4)) % 4][i % 4] =
        crypto_aes_sbox[((SIMDVec *)&s)->m128_u8[i]];
  for (i = 0; i < 16; ++i)
    ((SIMDVec *)&s)->m128_u8[i] =
        v[i / 4][i % 4] ^ ((SIMDVec *)&rk)->m128_u8[i];
  return s;
}

// Rotates right (circular right shift) value by "amount" positions
static inline uint32_t rotr(uint32_t value, uint32_t amount) {
  return (value >> amount) | (value << ((32 - amount) & 31));
}

static inline uint64_t MUL(uint32_t a, uint32_t b) {
  return (uint64_t)a * (uint64_t)b;
}

// From BearSSL. Performs a 32-bit->64-bit carryless/polynomial
// long multiply.
//
// This implementation was chosen because it is reasonably fast
// without a lookup table or branching.
//
// This does it by splitting up the bits in a way that they
// would not carry, then combine them together with xor (a
// carryless add).
//
// https://www.bearssl.org/gitweb/?p=BearSSL;a=blob;f=src/hash/ghash_ctmul.c;h=3623202;hb=5f045c7#l164
static uint64_t clmul_32(uint32_t x, uint32_t y) {
  uint32_t x0, x1, x2, x3;
  uint32_t y0, y1, y2, y3;
  uint64_t z0, z1, z2, z3;

  x0 = x & (uint32_t)0x11111111;
  x1 = x & (uint32_t)0x22222222;
  x2 = x & (uint32_t)0x44444444;
  x3 = x & (uint32_t)0x88888888;
  y0 = y & (uint32_t)0x11111111;
  y1 = y & (uint32_t)0x22222222;
  y2 = y & (uint32_t)0x44444444;
  y3 = y & (uint32_t)0x88888888;
  z0 = MUL(x0, y0) ^ MUL(x1, y3) ^ MUL(x2, y2) ^ MUL(x3, y1);
  z1 = MUL(x0, y1) ^ MUL(x1, y0) ^ MUL(x2, y3) ^ MUL(x3, y2);
  z2 = MUL(x0, y2) ^ MUL(x1, y1) ^ MUL(x2, y0) ^ MUL(x3, y3);
  z3 = MUL(x0, y3) ^ MUL(x1, y2) ^ MUL(x2, y1) ^ MUL(x3, y0);
  z0 &= (uint64_t)0x1111111111111111;
  z1 &= (uint64_t)0x2222222222222222;
  z2 &= (uint64_t)0x4444444444444444;
  z3 &= (uint64_t)0x8888888888888888;
  return z0 | z1 | z2 | z3;
}

// Performs a 64x64->128-bit carryless/polynomial long
// multiply, using the above routine to calculate the
// subproducts needed for the full-size multiply.
//
// This uses the Karatsuba algorithm.
//
// Normally, the Karatsuba algorithm isn't beneficial
// until very large numbers due to carry tracking and
// multiplication being relatively cheap.
//
// However, we have no carries and multiplication is
// definitely not cheap, so the Karatsuba algorithm is
// a low cost and easy optimization.
//
// https://en.m.wikipedia.org/wiki/Karatsuba_algorithm
//
// Note that addition and subtraction are both
// performed with xor, since all operations are
// carryless.
//
// The comments represent the actual mathematical
// operations being performed (instead of the bitwise
// operations) and to reflect the linked Wikipedia article.
static std::pair<uint64_t, uint64_t> clmul_64(uint64_t x, uint64_t y) {
  // B = 2
  // m = 32
  // x = (x1 * B^m) + x0
  uint32_t x0 = x & 0xffffffff;
  uint32_t x1 = x >> 32;
  // y = (y1 * B^m) + y0
  uint32_t y0 = y & 0xffffffff;
  uint32_t y1 = y >> 32;

  // z0 = x0 * y0
  uint64_t z0 = clmul_32(x0, y0);
  // z2 = x1 * y1
  uint64_t z2 = clmul_32(x1, y1);
  // z1 = (x0 + x1) * (y0 + y1) - z0 - z2
  uint64_t z1 = clmul_32(x0 ^ x1, y0 ^ y1) ^ z0 ^ z2;

  // xy = z0 + (z1 * B^m) + (z2 * B^2m)
  // note: z1 is split between the low and high halves
  uint64_t xy0 = z0 ^ (z1 << 32);
  uint64_t xy1 = z2 ^ (z1 >> 32);

  return std::make_pair(xy0, xy1);
}

/* MMX */
result_t test_mm_empty(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  return TEST_SUCCESS;
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

/* SSE */
result_t test_mm_add_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const float *_a = impl.test_cases_float_pointer1;
  const float *_b = impl.test_cases_float_pointer2;
  float dx = _a[0] + _b[0];
  float dy = _a[1] + _b[1];
  float dz = _a[2] + _b[2];
  float dw = _a[3] + _b[3];

  __m128 a = load_m128(_a);
  __m128 b = load_m128(_b);
  __m128 c = _mm_add_ps(a, b);
  return validate_float(c, dx, dy, dz, dw);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_add_ss(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const float *_a = impl.test_cases_float_pointer1;
  const float *_b = impl.test_cases_float_pointer1;

  float f0 = _a[0] + _b[0];
  float f1 = _a[1];
  float f2 = _a[2];
  float f3 = _a[3];

  __m128 a = _mm_load_ps(_a);
  __m128 b = _mm_load_ps(_b);
  __m128 c = _mm_add_ss(a, b);

  return validate_float(c, f0, f1, f2, f3);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_and_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const float *_a = impl.test_cases_float_pointer1;
  const float *_b = impl.test_cases_float_pointer2;
  __m128 a = load_m128(_a);
  __m128 b = load_m128(_b);
  __m128 c = _mm_and_ps(a, b);
  // now for the assertion...
  const uint32_t *ia = (const uint32_t *)&a;
  const uint32_t *ib = (const uint32_t *)&b;
  uint32_t r[4];
  r[0] = ia[0] & ib[0];
  r[1] = ia[1] & ib[1];
  r[2] = ia[2] & ib[2];
  r[3] = ia[3] & ib[3];
  __m128i ret = do_mm_set_epi32(r[3], r[2], r[1], r[0]);
  result_t res = VALIDATE_INT32_M128(*(const __m128i *)&c, r);
  if (res) {
    res = VALIDATE_INT32_M128(ret, r);
  }
  return res;
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

// r0 := ~a0 & b0
// r1 := ~a1 & b1
// r2 := ~a2 & b2
// r3 := ~a3 & b3
result_t test_mm_andnot_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const float *_a = impl.test_cases_float_pointer1;
  const float *_b = impl.test_cases_float_pointer2;

  __m128 a = load_m128(_a);
  __m128 b = load_m128(_b);
  __m128 c = _mm_andnot_ps(a, b);
  // now for the assertion...
  const uint32_t *ia = (const uint32_t *)&a;
  const uint32_t *ib = (const uint32_t *)&b;
  uint32_t r[4];
  r[0] = ~ia[0] & ib[0];
  r[1] = ~ia[1] & ib[1];
  r[2] = ~ia[2] & ib[2];
  r[3] = ~ia[3] & ib[3];
  __m128i ret = do_mm_set_epi32(r[3], r[2], r[1], r[0]);
  result_t res = TEST_FAIL;
  res = VALIDATE_INT32_M128(*(const __m128i *)&c, r);
  if (res) {
    res = VALIDATE_INT32_M128(ret, r);
  }
  return res;
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_avg_pu16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const uint16_t *_a = (const uint16_t *)impl.test_cases_int_pointer1;
  const uint16_t *_b = (const uint16_t *)impl.test_cases_int_pointer2;
  uint16_t _c[4];
  _c[0] = (_a[0] + _b[0] + 1) >> 1;
  _c[1] = (_a[1] + _b[1] + 1) >> 1;
  _c[2] = (_a[2] + _b[2] + 1) >> 1;
  _c[3] = (_a[3] + _b[3] + 1) >> 1;

  __m64 a = load_m64(_a);
  __m64 b = load_m64(_b);
  __m64 c = _mm_avg_pu16(a, b);
  return VALIDATE_UINT16_M64(c, _c);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_avg_pu8(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const uint8_t *_a = (const uint8_t *)impl.test_cases_int_pointer1;
  const uint8_t *_b = (const uint8_t *)impl.test_cases_int_pointer2;
  uint8_t d[8];
  d[0] = (_a[0] + _b[0] + 1) >> 1;
  d[1] = (_a[1] + _b[1] + 1) >> 1;
  d[2] = (_a[2] + _b[2] + 1) >> 1;
  d[3] = (_a[3] + _b[3] + 1) >> 1;
  d[4] = (_a[4] + _b[4] + 1) >> 1;
  d[5] = (_a[5] + _b[5] + 1) >> 1;
  d[6] = (_a[6] + _b[6] + 1) >> 1;
  d[7] = (_a[7] + _b[7] + 1) >> 1;

  __m64 a = load_m64(_a);
  __m64 b = load_m64(_b);
  __m64 c = _mm_avg_pu8(a, b);

  return VALIDATE_UINT8_M64(c, d);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_cmpeq_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   const float *_b = impl.test_cases_float_pointer2;
  //   __m128 a = load_m128(_a);
  //   __m128 b = load_m128(_b);
  //
  //   int32_t result[4];
  //   result[0] = _a[0] == _b[0] ? -1 : 0;
  //   result[1] = _a[1] == _b[1] ? -1 : 0;
  //   result[2] = _a[2] == _b[2] ? -1 : 0;
  //   result[3] = _a[3] == _b[3] ? -1 : 0;
  //
  //   __m128 ret = _mm_cmpeq_ps(a, b);
  //   __m128i iret = *(const __m128i *)&ret;
  //   return VALIDATE_INT32_M128(iret, result);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cmpeq_ss(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   const float *_b = impl.test_cases_float_pointer2;
  //   __m128 a = load_m128(_a);
  //   __m128 b = load_m128(_b);
  //
  //   float result[4];
  //   result[0] = _a[0] == _b[0] ? ALL_BIT_1_32 : 0;
  //   result[1] = _a[1];
  //   result[2] = _a[2];
  //   result[3] = _a[3];
  //
  //   __m128 ret = _mm_cmpeq_ss(a, b);
  //   return validate_float(ret, result[0], result[1], result[2], result[3]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cmpge_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   const float *_b = impl.test_cases_float_pointer2;
  //   __m128 a = load_m128(_a);
  //   __m128 b = load_m128(_b);
  //
  //   int32_t result[4];
  //   result[0] = _a[0] >= _b[0] ? -1 : 0;
  //   result[1] = _a[1] >= _b[1] ? -1 : 0;
  //   result[2] = _a[2] >= _b[2] ? -1 : 0;
  //   result[3] = _a[3] >= _b[3] ? -1 : 0;
  //
  //   __m128 ret = _mm_cmpge_ps(a, b);
  //   __m128i iret = *(const __m128i *)&ret;
  //   return VALIDATE_INT32_M128(iret, result);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cmpge_ss(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   const float *_b = impl.test_cases_float_pointer2;
  //   __m128 a = load_m128(_a);
  //   __m128 b = load_m128(_b);
  //
  //   float result[4];
  //   result[0] = _a[0] >= _b[0] ? ALL_BIT_1_32 : 0;
  //   result[1] = _a[1];
  //   result[2] = _a[2];
  //   result[3] = _a[3];
  //
  //   __m128 ret = _mm_cmpge_ss(a, b);
  //   return validate_float(ret, result[0], result[1], result[2], result[3]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cmpgt_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   const float *_b = impl.test_cases_float_pointer2;
  //   __m128 a = load_m128(_a);
  //   __m128 b = load_m128(_b);
  //
  //   int32_t result[4];
  //   result[0] = _a[0] > _b[0] ? -1 : 0;
  //   result[1] = _a[1] > _b[1] ? -1 : 0;
  //   result[2] = _a[2] > _b[2] ? -1 : 0;
  //   result[3] = _a[3] > _b[3] ? -1 : 0;
  //
  //   __m128 ret = _mm_cmpgt_ps(a, b);
  //   __m128i iret = *(const __m128i *)&ret;
  //   return VALIDATE_INT32_M128(iret, result);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cmpgt_ss(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   const float *_b = impl.test_cases_float_pointer2;
  //   __m128 a = load_m128(_a);
  //   __m128 b = load_m128(_b);
  //
  //   float result[4];
  //   result[0] = _a[0] > _b[0] ? ALL_BIT_1_32 : 0;
  //   result[1] = _a[1];
  //   result[2] = _a[2];
  //   result[3] = _a[3];
  //
  //   __m128 ret = _mm_cmpgt_ss(a, b);
  //   return validate_float(ret, result[0], result[1], result[2], result[3]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cmple_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   const float *_b = impl.test_cases_float_pointer2;
  //   __m128 a = load_m128(_a);
  //   __m128 b = load_m128(_b);
  //
  //   int32_t result[4];
  //   result[0] = _a[0] <= _b[0] ? -1 : 0;
  //   result[1] = _a[1] <= _b[1] ? -1 : 0;
  //   result[2] = _a[2] <= _b[2] ? -1 : 0;
  //   result[3] = _a[3] <= _b[3] ? -1 : 0;
  //
  //   __m128 ret = _mm_cmple_ps(a, b);
  //   __m128i iret = *(const __m128i *)&ret;
  //   return VALIDATE_INT32_M128(iret, result);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cmple_ss(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   const float *_b = impl.test_cases_float_pointer2;
  //   __m128 a = load_m128(_a);
  //   __m128 b = load_m128(_b);
  //
  //   float result[4];
  //   result[0] = _a[0] <= _b[0] ? ALL_BIT_1_32 : 0;
  //   result[1] = _a[1];
  //   result[2] = _a[2];
  //   result[3] = _a[3];
  //
  //   __m128 ret = _mm_cmple_ss(a, b);
  //   return validate_float(ret, result[0], result[1], result[2], result[3]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cmplt_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   const float *_b = impl.test_cases_float_pointer2;
  //   __m128 a = load_m128(_a);
  //   __m128 b = load_m128(_b);
  //
  //   int32_t result[4];
  //   result[0] = _a[0] < _b[0] ? -1 : 0;
  //   result[1] = _a[1] < _b[1] ? -1 : 0;
  //   result[2] = _a[2] < _b[2] ? -1 : 0;
  //   result[3] = _a[3] < _b[3] ? -1 : 0;
  //
  //   __m128 ret = _mm_cmplt_ps(a, b);
  //   __m128i iret = *(const __m128i *)&ret;
  //   return VALIDATE_INT32_M128(iret, result);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cmplt_ss(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   const float *_b = impl.test_cases_float_pointer2;
  //
  //   __m128 a = load_m128(_a);
  //   __m128 b = load_m128(_b);
  //
  //   float result[4];
  //   result[0] = _a[0] < _b[0] ? ALL_BIT_1_32 : 0;
  //   result[1] = _a[1];
  //   result[2] = _a[2];
  //   result[3] = _a[3];
  //
  //   __m128 ret = _mm_cmplt_ss(a, b);
  //   return validate_float(ret, result[0], result[1], result[2], result[3]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cmpneq_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const float *_a = impl.test_cases_float_pointer1;
  const float *_b = impl.test_cases_float_pointer2;
  __m128 a = load_m128(_a);
  __m128 b = load_m128(_b);

  int32_t result[4];
  result[0] = _a[0] != _b[0] ? -1 : 0;
  result[1] = _a[1] != _b[1] ? -1 : 0;
  result[2] = _a[2] != _b[2] ? -1 : 0;
  result[3] = _a[3] != _b[3] ? -1 : 0;

  __m128 ret = _mm_cmpneq_ps(a, b);
  __m128i iret = *(const __m128i *)&ret;
  return VALIDATE_INT32_M128(iret, result);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_cmpneq_ss(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   const float *_b = impl.test_cases_float_pointer2;
  //   __m128 a = load_m128(_a);
  //   __m128 b = load_m128(_b);
  //
  //   float result[4];
  //   result[0] = _a[0] != _b[0] ? ALL_BIT_1_32 : 0;
  //   result[1] = _a[1];
  //   result[2] = _a[2];
  //   result[3] = _a[3];
  //
  //   __m128 ret = _mm_cmpneq_ss(a, b);
  //   return validate_float(ret, result[0], result[1], result[2], result[3]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cmpnge_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  // const float *_a = impl.test_cases_float_pointer1;
  // const float *_b = impl.test_cases_float_pointer2;
  // __m128 a = load_m128(_a);
  // __m128 b = load_m128(_b);

  // float _c[4];
  // _c[0] = !(_a[0] >= _b[0]) ? UINT32_MAX : 0;
  // _c[1] = !(_a[1] >= _b[1]) ? UINT32_MAX : 0;
  // _c[2] = !(_a[2] >= _b[2]) ? UINT32_MAX : 0;
  // _c[3] = !(_a[3] >= _b[3]) ? UINT32_MAX : 0;

  // __m128 c = _mm_cmpnge_ps(a, b);

  // return validate_float(c, _c[0], _c[1], _c[2], _c[3]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cmpnge_ss(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   const float *_b = impl.test_cases_float_pointer2;
  //   __m128 a = load_m128(_a);
  //   __m128 b = load_m128(_b);
  //
  //   float result[4];
  //   result[0] = !(_a[0] >= _b[0]) ? ALL_BIT_1_32 : 0;
  //   result[1] = _a[1];
  //   result[2] = _a[2];
  //   result[3] = _a[3];
  //
  //   __m128 ret = _mm_cmpnge_ss(a, b);
  //   return validate_float(ret, result[0], result[1], result[2], result[3]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cmpngt_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   const float *_b = impl.test_cases_float_pointer2;
  //   __m128 a = load_m128(_a);
  //   __m128 b = load_m128(_b);
  //
  //   float result[4];
  //   result[0] = !(_a[0] > _b[0]) ? ALL_BIT_1_32 : 0;
  //   result[1] = !(_a[1] > _b[1]) ? ALL_BIT_1_32 : 0;
  //   result[2] = !(_a[2] > _b[2]) ? ALL_BIT_1_32 : 0;
  //   result[3] = !(_a[3] > _b[3]) ? ALL_BIT_1_32 : 0;
  //
  //   __m128 ret = _mm_cmpngt_ps(a, b);
  //   return validate_float(ret, result[0], result[1], result[2], result[3]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cmpngt_ss(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   const float *_b = impl.test_cases_float_pointer2;
  //   __m128 a = load_m128(_a);
  //   __m128 b = load_m128(_b);
  //
  //   float result[4];
  //   result[0] = !(_a[0] > _b[0]) ? ALL_BIT_1_32 : 0;
  //   result[1] = _a[1];
  //   result[2] = _a[2];
  //   result[3] = _a[3];
  //
  //   __m128 ret = _mm_cmpngt_ss(a, b);
  //   return validate_float(ret, result[0], result[1], result[2], result[3]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cmpnle_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   const float *_b = impl.test_cases_float_pointer2;
  //   __m128 a = load_m128(_a);
  //   __m128 b = load_m128(_b);
  //
  //   float result[4];
  //   result[0] = !(_a[0] <= _b[0]) ? ALL_BIT_1_32 : 0;
  //   result[1] = !(_a[1] <= _b[1]) ? ALL_BIT_1_32 : 0;
  //   result[2] = !(_a[2] <= _b[2]) ? ALL_BIT_1_32 : 0;
  //   result[3] = !(_a[3] <= _b[3]) ? ALL_BIT_1_32 : 0;
  //
  //   __m128 ret = _mm_cmpnle_ps(a, b);
  //   return validate_float(ret, result[0], result[1], result[2], result[3]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cmpnle_ss(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   const float *_b = impl.test_cases_float_pointer2;
  //   __m128 a = load_m128(_a);
  //   __m128 b = load_m128(_b);
  //
  //   float result[4];
  //   result[0] = !(_a[0] <= _b[0]) ? ALL_BIT_1_32 : 0;
  //   result[1] = _a[1];
  //   result[2] = _a[2];
  //   result[3] = _a[3];
  //
  //   __m128 ret = _mm_cmpnle_ss(a, b);
  //   return validate_float(ret, result[0], result[1], result[2], result[3]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cmpnlt_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   const float *_b = impl.test_cases_float_pointer2;
  //   __m128 a = load_m128(_a);
  //   __m128 b = load_m128(_b);
  //
  //   float result[4];
  //   result[0] = !(_a[0] < _b[0]) ? ALL_BIT_1_32 : 0;
  //   result[1] = !(_a[1] < _b[1]) ? ALL_BIT_1_32 : 0;
  //   result[2] = !(_a[2] < _b[2]) ? ALL_BIT_1_32 : 0;
  //   result[3] = !(_a[3] < _b[3]) ? ALL_BIT_1_32 : 0;
  //
  //   __m128 ret = _mm_cmpnlt_ps(a, b);
  //   return validate_float(ret, result[0], result[1], result[2], result[3]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cmpnlt_ss(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   const float *_b = impl.test_cases_float_pointer2;
  //   __m128 a = load_m128(_a);
  //   __m128 b = load_m128(_b);
  //
  //   float result[4];
  //   result[0] = !(_a[0] < _b[0]) ? ALL_BIT_1_32 : 0;
  //   result[1] = _a[1];
  //   result[2] = _a[2];
  //   result[3] = _a[3];
  //
  //   __m128 ret = _mm_cmpnlt_ss(a, b);
  //   return validate_float(ret, result[0], result[1], result[2], result[3]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cmpord_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const float *_a = impl.test_cases_float_pointer1;
  const float *_b = impl.test_cases_float_pointer2;
  __m128 a = load_m128(_a);
  __m128 b = load_m128(_b);

  float _c[4];
  for (uint32_t i = 0; i < 4; i++) {
    _c[i] = cmp_noNaN(_a[i], _b[i]);
  }
  __m128 c = _mm_cmpord_ps(a, b);

  return validate_float(c, _c[0], _c[1], _c[2], _c[3]);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_cmpord_ss(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const float *_a = impl.test_cases_float_pointer1;
  const float *_b = impl.test_cases_float_pointer2;
  __m128 a = load_m128(_a);
  __m128 b = load_m128(_b);
  float _c[4];
  _c[0] = cmp_noNaN(_a[0], _b[0]);
  _c[1] = _a[1];
  _c[2] = _a[2];
  _c[3] = _a[3];
  __m128 c = _mm_cmpord_ss(a, b);

  return validate_float(c, _c[0], _c[1], _c[2], _c[3]);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_cmpunord_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const float *_a = impl.test_cases_float_pointer1;
  const float *_b = impl.test_cases_float_pointer2;
  __m128 a = load_m128(_a);
  __m128 b = load_m128(_b);

  float _c[4];
  for (uint32_t i = 0; i < 4; i++) {
    _c[i] = cmp_hasNaN(_a[i], _b[i]);
  }
  __m128 c = _mm_cmpunord_ps(a, b);

  return validate_float(c, _c[0], _c[1], _c[2], _c[3]);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_cmpunord_ss(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const float *_a = impl.test_cases_float_pointer1;
  const float *_b = impl.test_cases_float_pointer2;
  __m128 a = load_m128(_a);
  __m128 b = load_m128(_b);

  float _c[4];
  _c[0] = cmp_hasNaN(_a[0], _b[0]);
  _c[1] = _a[1];
  _c[2] = _a[2];
  _c[3] = _a[3];
  __m128 c = _mm_cmpunord_ss(a, b);

  return validate_float(c, _c[0], _c[1], _c[2], _c[3]);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_comieq_ss(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
// FIXME:
// The GCC does not implement _mm_comieq_ss correctly.
// See https://gcc.gnu.org/bugzilla/show_bug.cgi?id=98612 for more
// information.
#if defined(__GNUC__) && !defined(__clang__)
  return TEST_UNIMPL;
#else
  const float *_a = impl.test_cases_float_pointer1;
  const float *_b = impl.test_cases_float_pointer2;

  int32_t _c = comieq_ss(_a[0], _b[0]);
  __m128 a = load_m128(_a);
  __m128 b = load_m128(_b);
  int32_t c = _mm_comieq_ss(a, b);

  return _c == c ? TEST_SUCCESS : TEST_FAIL;
#endif
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_comige_ss(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const float *_a = impl.test_cases_float_pointer1;
  const float *_b = impl.test_cases_float_pointer2;
  __m128 a = load_m128(_a);
  __m128 b = load_m128(_b);

  int32_t result = comige_ss(_a[0], _b[0]);
  int32_t ret = _mm_comige_ss(a, b);

  return result == ret ? TEST_SUCCESS : TEST_FAIL;
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_comigt_ss(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const float *_a = impl.test_cases_float_pointer1;
  const float *_b = impl.test_cases_float_pointer2;
  __m128 a = load_m128(_a);
  __m128 b = load_m128(_b);

  int32_t result = comigt_ss(_a[0], _b[0]);
  int32_t ret = _mm_comigt_ss(a, b);

  return result == ret ? TEST_SUCCESS : TEST_FAIL;
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_comile_ss(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
// FIXME:
// The GCC does not implement _mm_comile_ss correctly.
// See https://gcc.gnu.org/bugzilla/show_bug.cgi?id=98612 for more
// information.
#if defined(__GNUC__) && !defined(__clang__)
  return TEST_UNIMPL;
#else
  const float *_a = impl.test_cases_float_pointer1;
  const float *_b = impl.test_cases_float_pointer2;
  __m128 a = load_m128(_a);
  __m128 b = load_m128(_b);

  int32_t result = comile_ss(_a[0], _b[0]);
  int32_t ret = _mm_comile_ss(a, b);

  return result == ret ? TEST_SUCCESS : TEST_FAIL;
#endif
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_comilt_ss(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
// FIXME:
// The GCC does not implement _mm_comilt_ss correctly.
// See https://gcc.gnu.org/bugzilla/show_bug.cgi?id=98612 for more
// information.
#if defined(__GNUC__) && !defined(__clang__)
  return TEST_UNIMPL;
#else
  const float *_a = impl.test_cases_float_pointer1;
  const float *_b = impl.test_cases_float_pointer2;
  __m128 a = load_m128(_a);
  __m128 b = load_m128(_b);

  int32_t result = comilt_ss(_a[0], _b[0]);

  int32_t ret = _mm_comilt_ss(a, b);

  return result == ret ? TEST_SUCCESS : TEST_FAIL;
#endif
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_comineq_ss(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
// FIXME:
// The GCC does not implement _mm_comineq_ss correctly.
// See https://gcc.gnu.org/bugzilla/show_bug.cgi?id=98612 for more
// information.
#if defined(__GNUC__) && !defined(__clang__)
  return TEST_UNIMPL;
#else
  const float *_a = impl.test_cases_float_pointer1;
  const float *_b = impl.test_cases_float_pointer2;
  __m128 a = load_m128(_a);
  __m128 b = load_m128(_b);

  int32_t result = comineq_ss(_a[0], _b[0]);
  int32_t ret = _mm_comineq_ss(a, b);

  return result == ret ? TEST_SUCCESS : TEST_FAIL;
#endif
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_cvt_pi2ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   const int32_t *_b = impl.test_cases_int_pointer2;
  //
  //   float dx = (float)_b[0];
  //   float dy = (float)_b[1];
  //   float dz = _a[2];
  //   float dw = _a[3];
  //
  //   __m128 a = load_m128(_a);
  //   __m64 b = load_m64(_b);
  //   __m128 c = _mm_cvt_pi2ps(a, b);
  //
  //   return validate_float(c, dx, dy, dz, dw);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvt_ps2pi(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   int32_t d[2];
  //
  //   for (int idx = 0; idx < 2; idx++) {
  //     switch (iter & 0x3) {
  //     case 0:
  //       _MM_SET_ROUNDING_MODE(_MM_ROUND_NEAREST);
  //       d[idx] = (int32_t)(bankers_rounding(_a[idx]));
  //       break;
  //     case 1:
  //       _MM_SET_ROUNDING_MODE(_MM_ROUND_DOWN);
  //       d[idx] = (int32_t)(floorf(_a[idx]));
  //       break;
  //     case 2:
  //       _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
  //       d[idx] = (int32_t)(ceilf(_a[idx]));
  //       break;
  //     case 3:
  //       _MM_SET_ROUNDING_MODE(_MM_ROUND_TOWARD_ZERO);
  //       d[idx] = (int32_t)(_a[idx]);
  //       break;
  //     }
  //   }
  //
  //   __m128 a = load_m128(_a);
  //   __m64 ret = _mm_cvt_ps2pi(a);
  //
  //   return VALIDATE_INT32_M64(ret, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvt_si2ss(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   const int32_t b = *impl.test_cases_int_pointer2;
  //
  //   float dx = (float)b;
  //   float dy = _a[1];
  //   float dz = _a[2];
  //   float dw = _a[3];
  //
  //   __m128 a = load_m128(_a);
  //   __m128 c = _mm_cvt_si2ss(a, b);
  //
  //   return validate_float(c, dx, dy, dz, dw);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvt_ss2si(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   int32_t d0;
  //
  //   switch (iter & 0x3) {
  //   case 0:
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_NEAREST);
  //     d0 = (int32_t)(bankers_rounding(_a[0]));
  //     break;
  //   case 1:
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_DOWN);
  //     d0 = (int32_t)(floorf(_a[0]));
  //     break;
  //   case 2:
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
  //     d0 = (int32_t)(ceilf(_a[0]));
  //     break;
  //   case 3:
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_TOWARD_ZERO);
  //     d0 = (int32_t)(_a[0]);
  //     break;
  //   }
  //
  //   __m128 a = load_m128(_a);
  //   int32_t ret = _mm_cvt_ss2si(a);
  //   return ret == d0 ? TEST_SUCCESS : TEST_FAIL;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtpi16_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  //
  //   float dx = (float)_a[0];
  //   float dy = (float)_a[1];
  //   float dz = (float)_a[2];
  //   float dw = (float)_a[3];
  //
  //   __m64 a = load_m64(_a);
  //   __m128 c = _mm_cvtpi16_ps(a);
  //
  //   return validate_float(c, dx, dy, dz, dw);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtpi32_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   const int32_t *_b = (const int32_t *)impl.test_cases_int_pointer2;
  //
  //   float dx = (float)_b[0];
  //   float dy = (float)_b[1];
  //   float dz = _a[2];
  //   float dw = _a[3];
  //
  //   __m128 a = load_m128(_a);
  //   __m64 b = load_m64(_b);
  //   __m128 c = _mm_cvtpi32_ps(a, b);
  //
  //   return validate_float(c, dx, dy, dz, dw);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtpi32x2_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int32_t *_a = (const int32_t *)impl.test_cases_int_pointer1;
  //   const int32_t *_b = (const int32_t *)impl.test_cases_int_pointer2;
  //
  //   float dx = (float)_a[0];
  //   float dy = (float)_a[1];
  //   float dz = (float)_b[0];
  //   float dw = (float)_b[1];
  //
  //   __m64 a = load_m64(_a);
  //   __m64 b = load_m64(_b);
  //   __m128 c = _mm_cvtpi32x2_ps(a, b);
  //
  //   return validate_float(c, dx, dy, dz, dw);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtpi8_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int8_t *_a = (const int8_t *)impl.test_cases_int_pointer1;
  //
  //   float dx = (float)_a[0];
  //   float dy = (float)_a[1];
  //   float dz = (float)_a[2];
  //   float dw = (float)_a[3];
  //
  //   __m64 a = load_m64(_a);
  //   __m128 c = _mm_cvtpi8_ps(a);
  //
  //   return validate_float(c, dx, dy, dz, dw);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtps_pi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   int16_t rnd[4];
  //
  //   for (int i = 0; i < 4; i++) {
  //     if ((float)INT16_MAX <= _a[i] && _a[i] <= (float)INT32_MAX) {
  //       rnd[i] = INT16_MAX;
  //     } else if (INT16_MIN < _a[i] && _a[i] < INT16_MAX) {
  //       switch (iter & 0x3) {
  //       case 0:
  //         _MM_SET_ROUNDING_MODE(_MM_ROUND_NEAREST);
  //         rnd[i] = (int16_t)bankers_rounding(_a[i]);
  //         break;
  //       case 1:
  //         _MM_SET_ROUNDING_MODE(_MM_ROUND_DOWN);
  //         rnd[i] = (int16_t)floorf(_a[i]);
  //         break;
  //       case 2:
  //         _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
  //         rnd[i] = (int16_t)ceilf(_a[i]);
  //         break;
  //       case 3:
  //         _MM_SET_ROUNDING_MODE(_MM_ROUND_TOWARD_ZERO);
  //         rnd[i] = (int16_t)_a[i];
  //         break;
  //       }
  //     } else {
  //       rnd[i] = INT16_MIN;
  //     }
  //   }
  //
  //   __m128 a = load_m128(_a);
  //   __m64 ret = _mm_cvtps_pi16(a);
  //   return VALIDATE_INT16_M64(ret, rnd);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtps_pi32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   int32_t d[2];
  //
  //   switch (iter & 0x3) {
  //   case 0:
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_NEAREST);
  //     d[0] = (int32_t)bankers_rounding(_a[0]);
  //     d[1] = (int32_t)bankers_rounding(_a[1]);
  //     break;
  //   case 1:
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_DOWN);
  //     d[0] = (int32_t)floorf(_a[0]);
  //     d[1] = (int32_t)floorf(_a[1]);
  //     break;
  //   case 2:
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
  //     d[0] = (int32_t)ceilf(_a[0]);
  //     d[1] = (int32_t)ceilf(_a[1]);
  //     break;
  //   case 3:
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_TOWARD_ZERO);
  //     d[0] = (int32_t)_a[0];
  //     d[1] = (int32_t)_a[1];
  //     break;
  //   }
  //
  //   __m128 a = load_m128(_a);
  //   __m64 ret = _mm_cvtps_pi32(a);
  //
  //   return VALIDATE_INT32_M64(ret, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtps_pi8(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   int8_t rnd[8] = {};
  //
  //   for (int i = 0; i < 4; i++) {
  //     if ((float)INT8_MAX <= _a[i] && _a[i] <= (float)INT32_MAX) {
  //       rnd[i] = INT8_MAX;
  //     } else if (INT8_MIN < _a[i] && _a[i] < INT8_MAX) {
  //       switch (iter & 0x3) {
  //       case 0:
  //         _MM_SET_ROUNDING_MODE(_MM_ROUND_NEAREST);
  //         rnd[i] = (int8_t)bankers_rounding(_a[i]);
  //         break;
  //       case 1:
  //         _MM_SET_ROUNDING_MODE(_MM_ROUND_DOWN);
  //         rnd[i] = (int8_t)floorf(_a[i]);
  //         break;
  //       case 2:
  //         _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
  //         rnd[i] = (int8_t)ceilf(_a[i]);
  //         break;
  //       case 3:
  //         _MM_SET_ROUNDING_MODE(_MM_ROUND_TOWARD_ZERO);
  //         rnd[i] = (int8_t)_a[i];
  //         break;
  //       }
  //     } else {
  //       rnd[i] = INT8_MIN;
  //     }
  //   }
  //
  //   __m128 a = load_m128(_a);
  //   __m64 ret = _mm_cvtps_pi8(a);
  //   return VALIDATE_INT8_M64(ret, rnd);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtpu16_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const uint16_t *_a = (const uint16_t *)impl.test_cases_int_pointer1;
  //
  //   float dx = (float)_a[0];
  //   float dy = (float)_a[1];
  //   float dz = (float)_a[2];
  //   float dw = (float)_a[3];
  //
  //   __m64 a = load_m64(_a);
  //   __m128 c = _mm_cvtpu16_ps(a);
  //
  //   return validate_float(c, dx, dy, dz, dw);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtpu8_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const uint8_t *_a = (const uint8_t *)impl.test_cases_int_pointer1;
  //
  //   float dx = (float)_a[0];
  //   float dy = (float)_a[1];
  //   float dz = (float)_a[2];
  //   float dw = (float)_a[3];
  //
  //   __m64 a = load_m64(_a);
  //   __m128 c = _mm_cvtpu8_ps(a);
  //
  //   return validate_float(c, dx, dy, dz, dw);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtsi32_ss(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   const int32_t b = *impl.test_cases_int_pointer2;
  //
  //   float dx = (float)b;
  //   float dy = _a[1];
  //   float dz = _a[2];
  //   float dw = _a[3];
  //
  //   __m128 a = load_m128(_a);
  //   __m128 c = _mm_cvtsi32_ss(a, b);
  //
  //   return validate_float(c, dx, dy, dz, dw);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtsi64_ss(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   const int64_t b = *(int64_t *)impl.test_cases_int_pointer2;
  //
  //   float dx = (float)b;
  //   float dy = _a[1];
  //   float dz = _a[2];
  //   float dw = _a[3];
  //
  //   __m128 a = load_m128(_a);
  //   __m128 c = _mm_cvtsi64_ss(a, b);
  //
  //   return validate_float(c, dx, dy, dz, dw);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtss_f32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //
  //   float f = _a[0];
  //
  //   __m128 a = load_m128(_a);
  //   float c = _mm_cvtss_f32(a);
  //
  //   return f == c ? TEST_SUCCESS : TEST_FAIL;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtss_si32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //
  //   int32_t d0;
  //   switch (iter & 0x3) {
  //   case 0:
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_NEAREST);
  //     d0 = (int32_t)(bankers_rounding(_a[0]));
  //     break;
  //   case 1:
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_DOWN);
  //     d0 = (int32_t)(floorf(_a[0]));
  //     break;
  //   case 2:
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
  //     d0 = (int32_t)(ceilf(_a[0]));
  //     break;
  //   case 3:
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_TOWARD_ZERO);
  //     d0 = (int32_t)(_a[0]);
  //     break;
  //   }
  //
  //   __m128 a = load_m128(_a);
  //   int32_t ret = _mm_cvtss_si32(a);
  //
  //   return ret == d0 ? TEST_SUCCESS : TEST_FAIL;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtss_si64(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //
  //   int64_t d0;
  //   switch (iter & 0x3) {
  //   case 0:
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_NEAREST);
  //     d0 = (int64_t)(bankers_rounding(_a[0]));
  //     break;
  //   case 1:
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_DOWN);
  //     d0 = (int64_t)(floorf(_a[0]));
  //     break;
  //   case 2:
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
  //     d0 = (int64_t)(ceilf(_a[0]));
  //     break;
  //   case 3:
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_TOWARD_ZERO);
  //     d0 = (int64_t)(_a[0]);
  //     break;
  //   }
  //
  //   __m128 a = load_m128(_a);
  //   int64_t ret = _mm_cvtss_si64(a);
  //
  //   return ret == d0 ? TEST_SUCCESS : TEST_FAIL;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtt_ps2pi(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   int32_t d[2];
  //
  //   d[0] = (int32_t)_a[0];
  //   d[1] = (int32_t)_a[1];
  //
  //   __m128 a = load_m128(_a);
  //   __m64 ret = _mm_cvtt_ps2pi(a);
  //
  //   return VALIDATE_INT32_M64(ret, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtt_ss2si(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //
  //   __m128 a = load_m128(_a);
  //   int ret = _mm_cvtt_ss2si(a);
  //
  //   return ret == (int32_t)_a[0] ? TEST_SUCCESS : TEST_FAIL;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvttps_pi32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   int32_t d[2];
  //
  //   d[0] = (int32_t)_a[0];
  //   d[1] = (int32_t)_a[1];
  //
  //   __m128 a = load_m128(_a);
  //   __m64 ret = _mm_cvttps_pi32(a);
  //
  //   return VALIDATE_INT32_M64(ret, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvttss_si32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //
  //   __m128 a = load_m128(_a);
  //   int ret = _mm_cvttss_si32(a);
  //
  //   return ret == (int32_t)_a[0] ? TEST_SUCCESS : TEST_FAIL;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvttss_si64(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //
  //   __m128 a = load_m128(_a);
  //   int64_t ret = _mm_cvttss_si64(a);
  //
  //   return ret == (int64_t)_a[0] ? TEST_SUCCESS : TEST_FAIL;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_div_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   const float *_b = impl.test_cases_float_pointer2;
  //   float f0 = _a[0] / _b[0];
  //   float f1 = _a[1] / _b[1];
  //   float f2 = _a[2] / _b[2];
  //   float f3 = _a[3] / _b[3];
  //
  //   __m128 a = load_m128(_a);
  //   __m128 b = load_m128(_b);
  //   __m128 c = _mm_div_ps(a, b);
  //
  //   return validate_float(c, f0, f1, f2, f3);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_div_ss(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   const float *_b = impl.test_cases_float_pointer2;
  //
  //   float d0 = _a[0] / _b[0];
  //   float d1 = _a[1];
  //   float d2 = _a[2];
  //   float d3 = _a[3];
  //
  //   __m128 a = load_m128(_a);
  //   __m128 b = load_m128(_b);
  //   __m128 c = _mm_div_ss(a, b);
  //
  //   return validate_float(c, d0, d1, d2, d3);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_extract_pi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  // FIXME GCC has bug on "_mm_extract_pi16" intrinsics. We will enable this
  // test when GCC fix this bug.
  // see https://gcc.gnu.org/bugzilla/show_bug.cgi?id=98495 for more
  // information
  // #if defined(__clang__) || defined(_MSC_VER)
  //   uint64_t *_a = (uint64_t *)impl.test_cases_int_pointer1;
  //   const int idx = iter & 0x3;
  //
  //   __m64 a = load_m64(_a);
  //   int c;
  //   switch (idx) {
  //   case 0:
  //     c = _mm_extract_pi16(a, 0);
  //     break;
  //   case 1:
  //     c = _mm_extract_pi16(a, 1);
  //     break;
  //   case 2:
  //     c = _mm_extract_pi16(a, 2);
  //     break;
  //   case 3:
  //     c = _mm_extract_pi16(a, 3);
  //     break;
  //   }
  //
  //   ASSERT_RETURN((uint64_t)c == ((*_a >> (idx * 16)) & 0xFFFF));
  //   ASSERT_RETURN(0 == ((uint64_t)c & 0xFFFF0000));
  //   return TEST_SUCCESS;
  // #else
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
  // #endif
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_malloc(const SSE2RVV_TEST_IMPL &impl, uint32_t iter);
// #ifdef ENABLE_TEST_ALL

result_t test_mm_free(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   /* We verify _mm_malloc first, and there is no need to check _mm_free .
  //   */ return test_mm_malloc(impl, iter);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_get_flush_zero_mode(const SSE2RVV_TEST_IMPL &impl,
                                     // #ifdef ENABLE_TEST_ALL
                                     uint32_t iter) {
  //   int res_flush_zero_on, res_flush_zero_off;
  //   _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
  //   res_flush_zero_on = _MM_GET_FLUSH_ZERO_MODE() == _MM_FLUSH_ZERO_ON;
  //   _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_OFF);
  //   res_flush_zero_off = _MM_GET_FLUSH_ZERO_MODE() == _MM_FLUSH_ZERO_OFF;
  //
  //   return (res_flush_zero_on && res_flush_zero_off) ? TEST_SUCCESS :
  //   TEST_FAIL;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_get_rounding_mode(const SSE2RVV_TEST_IMPL &impl,
                                   // #ifdef ENABLE_TEST_ALL
                                   uint32_t iter) {
  //   int res_toward_zero, res_to_neg_inf, res_to_pos_inf, res_nearest;
  //   _MM_SET_ROUNDING_MODE(_MM_ROUND_TOWARD_ZERO);
  //   res_toward_zero = _MM_GET_ROUNDING_MODE() == _MM_ROUND_TOWARD_ZERO ? 1 :
  //   0; _MM_SET_ROUNDING_MODE(_MM_ROUND_DOWN); res_to_neg_inf =
  //   _MM_GET_ROUNDING_MODE() == _MM_ROUND_DOWN ? 1 : 0;
  //   _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
  //   res_to_pos_inf = _MM_GET_ROUNDING_MODE() == _MM_ROUND_UP ? 1 : 0;
  //   _MM_SET_ROUNDING_MODE(_MM_ROUND_NEAREST);
  //   res_nearest = _MM_GET_ROUNDING_MODE() == _MM_ROUND_NEAREST ? 1 : 0;
  //
  //   if (res_toward_zero && res_to_neg_inf && res_to_pos_inf && res_nearest) {
  //     return TEST_SUCCESS;
  //   } else {
  //     return TEST_FAIL;
  //   }
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_getcsr(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  // store original csr value for post test restoring
  //   unsigned int originalCsr = _mm_getcsr();
  //
  //   unsigned int roundings[] = {_MM_ROUND_TOWARD_ZERO, _MM_ROUND_DOWN,
  //                               _MM_ROUND_UP, _MM_ROUND_NEAREST};
  //   for (size_t i = 0; i < sizeof(roundings) / sizeof(roundings[0]); i++) {
  //     _mm_setcsr(_mm_getcsr() | roundings[i]);
  //     if ((_mm_getcsr() & roundings[i]) != roundings[i]) {
  //       return TEST_FAIL;
  //     }
  //   }
  //
  // restore original csr value for remaining tests
  //   _mm_setcsr(originalCsr);
  //
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_insert_pi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  //   const int16_t insert = (int16_t)impl.test_cases_ints[iter];
  //   __m64 a;
  //   __m64 b;
  //
  // #define TEST_IMPL(IDX)
  //   int16_t d##IDX[4];
  //   for (int i = 0; i < 4; i++) {
  //     d##IDX[i] = _a[i];
  //   }
  //   d##IDX[IDX] = insert;
  //
  //   a = load_m64(_a);
  //   b = _mm_insert_pi16(a, insert, IDX);
  //   CHECK_RESULT(VALIDATE_INT16_M64(b, d##IDX))
  //
  //   IMM_4_ITER
  // #undef TEST_IMPL
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_load_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const float *addr = impl.test_cases_float_pointer1;

  __m128 ret = _mm_load_ps(addr);

  return validate_float(ret, addr[0], addr[1], addr[2], addr[3]);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_load_ps1(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const float *addr = impl.test_cases_float_pointer1;

  __m128 ret = _mm_load_ps1(addr);

  return validate_float(ret, addr[0], addr[0], addr[0], addr[0]);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_load_ss(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const float *addr = impl.test_cases_float_pointer1;

  __m128 ret = _mm_load_ss(addr);

  return validate_float(ret, addr[0], 0, 0, 0);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_load1_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const float *p = impl.test_cases_float_pointer1;
  __m128 a = _mm_load1_ps(p);
  return validate_float(a, p[0], p[0], p[0], p[0]);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_loadh_pi(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *p1 = impl.test_cases_float_pointer1;
  //   const float *p2 = impl.test_cases_float_pointer2;
  //   const __m64 *b = (const __m64 *)p2;
  //   __m128 a = _mm_load_ps(p1);
  //   __m128 c = _mm_loadh_pi(a, b);
  //
  //   return validate_float(c, p1[0], p1[1], p2[0], p2[1]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_loadl_pi(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *p1 = impl.test_cases_float_pointer1;
  //   const float *p2 = impl.test_cases_float_pointer2;
  //   __m128 a = _mm_load_ps(p1);
  //   const __m64 *b = (const __m64 *)p2;
  //   __m128 c = _mm_loadl_pi(a, b);
  //
  //   return validate_float(c, p2[0], p2[1], p1[2], p1[3]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_loadr_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *addr = impl.test_cases_float_pointer1;
  //
  //   __m128 ret = _mm_loadr_ps(addr);
  //
  //   return validate_float(ret, addr[3], addr[2], addr[1], addr[0]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_loadu_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *addr = impl.test_cases_float_pointer1;
  //
  //   __m128 ret = _mm_loadu_ps(addr);
  //
  //   return validate_float(ret, addr[0], addr[1], addr[2], addr[3]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_loadu_si16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  // The GCC version before 11 does not implement intrinsic function
  // _mm_loadu_si16. Check https://gcc.gnu.org/bugzilla/show_bug.cgi?id=95483
  // for more information.
  // #if (defined(__GNUC__) && !defined(__clang__)) && (__GNUC__ <= 10)
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
  // #else
  //   const int16_t *addr = (const int16_t *)impl.test_cases_int_pointer1;
  //
  //   __m128i ret = _mm_loadu_si16((const void *)addr);
  //
  //   return validate_int16(ret, addr[0], 0, 0, 0, 0, 0, 0, 0);
  // #endif
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_loadu_si64(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  // Versions of GCC prior to 9 do not implement intrinsic function
  // _mm_loadu_si64. Check https://gcc.gnu.org/bugzilla/show_bug.cgi?id=78782
  // for more information.
  // #if (defined(__GNUC__) && !defined(__clang__)) && (__GNUC__ < 9)
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
  // #else
  //   const int64_t *addr = (const int64_t *)impl.test_cases_int_pointer1;
  //
  //   __m128i ret = _mm_loadu_si64((const void *)addr);
  //
  //   return validate_int64(ret, addr[0], 0);
  // #endif
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_malloc(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const size_t *a = (const size_t *)impl.test_cases_int_pointer1;
  //   const size_t *b = (const size_t *)impl.test_cases_int_pointer2;
  //   size_t size = *a % (1024 * 16) + 1;
  //   size_t align = 2 << (*b % 5);
  //
  //   void *p = _mm_malloc(size, align);
  //   if (!p)
  //     return TEST_FAIL;
  //   result_t res = (((uintptr_t)p % align) == 0) ? TEST_SUCCESS : TEST_FAIL;
  //   _mm_free(p);
  //   return res;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_maskmove_si64(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const uint8_t *_a = (const uint8_t *)impl.test_cases_int_pointer1;
  //   const uint8_t *_mask = (const uint8_t *)impl.test_cases_int_pointer2;
  //   char mem_addr[16];
  //
  //   const __m64 *a = (const __m64 *)_a;
  //   const __m64 *mask = (const __m64 *)_mask;
  //   _mm_maskmove_si64(*a, *mask, (char *)mem_addr);
  //
  //   for (int i = 0; i < 8; i++) {
  //     if (_mask[i] >> 7) {
  //       ASSERT_RETURN(_a[i] == (uint8_t)mem_addr[i]);
  //     }
  //   }
  //
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_m_maskmovq(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  // return test_mm_maskmove_si64(impl, iter);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_max_pi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  //   const int16_t *_b = (const int16_t *)impl.test_cases_int_pointer2;
  //   int16_t c[4];
  //
  //   c[0] = _a[0] > _b[0] ? _a[0] : _b[0];
  //   c[1] = _a[1] > _b[1] ? _a[1] : _b[1];
  //   c[2] = _a[2] > _b[2] ? _a[2] : _b[2];
  //   c[3] = _a[3] > _b[3] ? _a[3] : _b[3];
  //
  //   __m64 a = load_m64(_a);
  //   __m64 b = load_m64(_b);
  //   __m64 ret = _mm_max_pi16(a, b);
  //   return VALIDATE_INT16_M64(ret, c);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_max_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   const float *_b = impl.test_cases_float_pointer2;
  //   float c[4];
  //
  //   c[0] = _a[0] > _b[0] ? _a[0] : _b[0];
  //   c[1] = _a[1] > _b[1] ? _a[1] : _b[1];
  //   c[2] = _a[2] > _b[2] ? _a[2] : _b[2];
  //   c[3] = _a[3] > _b[3] ? _a[3] : _b[3];
  //
  //   __m128 a = load_m128(_a);
  //   __m128 b = load_m128(_b);
  //   __m128 ret = _mm_max_ps(a, b);
  //   return validate_float(ret, c[0], c[1], c[2], c[3]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_max_pu8(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const uint8_t *_a = (const uint8_t *)impl.test_cases_int_pointer1;
  //   const uint8_t *_b = (const uint8_t *)impl.test_cases_int_pointer2;
  //   uint8_t c[8];
  //
  //   c[0] = _a[0] > _b[0] ? _a[0] : _b[0];
  //   c[1] = _a[1] > _b[1] ? _a[1] : _b[1];
  //   c[2] = _a[2] > _b[2] ? _a[2] : _b[2];
  //   c[3] = _a[3] > _b[3] ? _a[3] : _b[3];
  //   c[4] = _a[4] > _b[4] ? _a[4] : _b[4];
  //   c[5] = _a[5] > _b[5] ? _a[5] : _b[5];
  //   c[6] = _a[6] > _b[6] ? _a[6] : _b[6];
  //   c[7] = _a[7] > _b[7] ? _a[7] : _b[7];
  //
  //   __m64 a = load_m64(_a);
  //   __m64 b = load_m64(_b);
  //   __m64 ret = _mm_max_pu8(a, b);
  //   return VALIDATE_UINT8_M64(ret, c);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_max_ss(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   const float *_b = impl.test_cases_float_pointer1;
  //
  //   float f0 = _a[0] > _b[0] ? _a[0] : _b[0];
  //   float f1 = _a[1];
  //   float f2 = _a[2];
  //   float f3 = _a[3];
  //
  //   __m128 a = _mm_load_ps(_a);
  //   __m128 b = _mm_load_ps(_b);
  //   __m128 c = _mm_max_ss(a, b);
  //
  //   return validate_float(c, f0, f1, f2, f3);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_min_pi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  //   const int16_t *_b = (const int16_t *)impl.test_cases_int_pointer2;
  //   int16_t c[4];
  //
  //   c[0] = _a[0] < _b[0] ? _a[0] : _b[0];
  //   c[1] = _a[1] < _b[1] ? _a[1] : _b[1];
  //   c[2] = _a[2] < _b[2] ? _a[2] : _b[2];
  //   c[3] = _a[3] < _b[3] ? _a[3] : _b[3];
  //
  //   __m64 a = load_m64(_a);
  //   __m64 b = load_m64(_b);
  //   __m64 ret = _mm_min_pi16(a, b);
  //   return VALIDATE_INT16_M64(ret, c);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_min_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   const float *_b = impl.test_cases_float_pointer2;
  //   float c[4];
  //
  //   c[0] = _a[0] < _b[0] ? _a[0] : _b[0];
  //   c[1] = _a[1] < _b[1] ? _a[1] : _b[1];
  //   c[2] = _a[2] < _b[2] ? _a[2] : _b[2];
  //   c[3] = _a[3] < _b[3] ? _a[3] : _b[3];
  //
  //   __m128 a = load_m128(_a);
  //   __m128 b = load_m128(_b);
  //   __m128 ret = _mm_min_ps(a, b);
  //   return validate_float(ret, c[0], c[1], c[2], c[3]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_min_pu8(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const uint8_t *_a = (const uint8_t *)impl.test_cases_int_pointer1;
  //   const uint8_t *_b = (const uint8_t *)impl.test_cases_int_pointer2;
  //   uint8_t c[8];
  //
  //   c[0] = _a[0] < _b[0] ? _a[0] : _b[0];
  //   c[1] = _a[1] < _b[1] ? _a[1] : _b[1];
  //   c[2] = _a[2] < _b[2] ? _a[2] : _b[2];
  //   c[3] = _a[3] < _b[3] ? _a[3] : _b[3];
  //   c[4] = _a[4] < _b[4] ? _a[4] : _b[4];
  //   c[5] = _a[5] < _b[5] ? _a[5] : _b[5];
  //   c[6] = _a[6] < _b[6] ? _a[6] : _b[6];
  //   c[7] = _a[7] < _b[7] ? _a[7] : _b[7];
  //
  //   __m64 a = load_m64(_a);
  //   __m64 b = load_m64(_b);
  //   __m64 ret = _mm_min_pu8(a, b);
  //   return VALIDATE_UINT8_M64(ret, c);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_min_ss(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   const float *_b = impl.test_cases_float_pointer2;
  //   float c;
  //
  //   c = _a[0] < _b[0] ? _a[0] : _b[0];
  //
  //   __m128 a = load_m128(_a);
  //   __m128 b = load_m128(_b);
  //   __m128 ret = _mm_min_ss(a, b);
  //
  //   return validate_float(ret, c, _a[1], _a[2], _a[3]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_move_ss(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   const float *_b = impl.test_cases_float_pointer2;
  //   __m128 a = load_m128(_a);
  //   __m128 b = load_m128(_b);
  //
  //   float result[4];
  //   result[0] = _b[0];
  //   result[1] = _a[1];
  //   result[2] = _a[2];
  //   result[3] = _a[3];
  //
  //   __m128 ret = _mm_move_ss(a, b);
  //   return validate_float(ret, result[0], result[1], result[2], result[3]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_movehl_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   const float *_b = impl.test_cases_float_pointer2;
  //
  //   float f0 = _b[2];
  //   float f1 = _b[3];
  //   float f2 = _a[2];
  //   float f3 = _a[3];
  //
  //   __m128 a = load_m128(_a);
  //   __m128 b = load_m128(_b);
  //   __m128 ret = _mm_movehl_ps(a, b);
  //
  //   return validate_float(ret, f0, f1, f2, f3);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_movelh_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   const float *_b = impl.test_cases_float_pointer2;
  //
  //   float f0 = _a[0];
  //   float f1 = _a[1];
  //   float f2 = _b[0];
  //   float f3 = _b[1];
  //
  //   __m128 a = load_m128(_a);
  //   __m128 b = load_m128(_b);
  //   __m128 ret = _mm_movelh_ps(a, b);
  //
  //   return validate_float(ret, f0, f1, f2, f3);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_movemask_pi8(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const uint8_t *_a = (const uint8_t *)impl.test_cases_int_pointer1;
  //   unsigned int _c = 0;
  //   for (int i = 0; i < 8; i++) {
  //     if (_a[i] & 0x80) {
  //       _c |= (1 << i);
  //     }
  //   }
  //
  //   const __m64 *a = (const __m64 *)_a;
  //   int c = _mm_movemask_pi8(*a);
  //
  //   ASSERT_RETURN((unsigned int)c == _c);
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_movemask_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *p = impl.test_cases_float_pointer1;
  //   int ret = 0;
  //
  //   const uint32_t *ip = (const uint32_t *)p;
  //   if (ip[0] & 0x80000000) {
  //     ret |= 1;
  //   }
  //   if (ip[1] & 0x80000000) {
  //     ret |= 2;
  //   }
  //   if (ip[2] & 0x80000000) {
  //     ret |= 4;
  //   }
  //   if (ip[3] & 0x80000000) {
  //     ret |= 8;
  //   }
  //   __m128 a = load_m128(p);
  //   int val = _mm_movemask_ps(a);
  //   return val == ret ? TEST_SUCCESS : TEST_FAIL;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_mul_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   const float *_b = impl.test_cases_float_pointer2;
  //   float dx = _a[0] * _b[0];
  //   float dy = _a[1] * _b[1];
  //   float dz = _a[2] * _b[2];
  //   float dw = _a[3] * _b[3];
  //
  //   __m128 a = load_m128(_a);
  //   __m128 b = load_m128(_b);
  //   __m128 c = _mm_mul_ps(a, b);
  //   return validate_float(c, dx, dy, dz, dw);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_mul_ss(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   const float *_b = impl.test_cases_float_pointer2;
  //
  //   float dx = _a[0] * _b[0];
  //   float dy = _a[1];
  //   float dz = _a[2];
  //   float dw = _a[3];
  //
  //   __m128 a = load_m128(_a);
  //   __m128 b = load_m128(_b);
  //   __m128 c = _mm_mul_ss(a, b);
  //   return validate_float(c, dx, dy, dz, dw);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_mulhi_pu16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const uint16_t *_a = (const uint16_t *)impl.test_cases_int_pointer1;
  //   const uint16_t *_b = (const uint16_t *)impl.test_cases_int_pointer2;
  //   uint16_t d[4];
  //   for (uint32_t i = 0; i < 4; i++) {
  //     uint32_t m = (uint32_t)_a[i] * (uint32_t)_b[i];
  //     d[i] = (uint16_t)(m >> 16);
  //   }
  //
  //   __m64 a = load_m64(_a);
  //   __m64 b = load_m64(_b);
  //   __m64 c = _mm_mulhi_pu16(a, b);
  //   return VALIDATE_UINT16_M64(c, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_or_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   const float *_b = impl.test_cases_float_pointer2;
  //   __m128 a = load_m128(_a);
  //   __m128 b = load_m128(_b);
  //   __m128 c = _mm_or_ps(a, b);
  // now for the assertion...
  //   const uint32_t *ia = (const uint32_t *)&a;
  //   const uint32_t *ib = (const uint32_t *)&b;
  //   uint32_t r[4];
  //   r[0] = ia[0] | ib[0];
  //   r[1] = ia[1] | ib[1];
  //   r[2] = ia[2] | ib[2];
  //   r[3] = ia[3] | ib[3];
  //   __m128i ret = do_mm_set_epi32(r[3], r[2], r[1], r[0]);
  //   result_t res = VALIDATE_INT32_M128(*(const __m128i *)&c, r);
  //   if (res) {
  //     res = VALIDATE_INT32_M128(ret, r);
  //   }
  //
  //   return res;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_m_pavgb(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  return test_mm_avg_pu8(impl, iter);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_m_pavgw(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  return test_mm_avg_pu16(impl, iter);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_m_pextrw(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  // return test_mm_extract_pi16(impl, iter);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_m_pinsrw(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  // return test_mm_insert_pi16(impl, iter);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_m_pmaxsw(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  // return test_mm_max_pi16(impl, iter);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_m_pmaxub(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  // return test_mm_max_pu8(impl, iter);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_m_pminsw(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  // return test_mm_min_pi16(impl, iter);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_m_pminub(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  // return test_mm_min_pu8(impl, iter);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_m_pmovmskb(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  // return test_mm_movemask_pi8(impl, iter);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_m_pmulhuw(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  // return test_mm_mulhi_pu16(impl, iter);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_prefetch(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   typedef struct {
  //     __m128 a;
  //     float r[4];
  //   } prefetch_test_t;
  //   prefetch_test_t test_vec[8] = {
  //       {
  //           _mm_set_ps(-0.1f, 0.2f, 0.3f, 0.4f),
  //           {0.4f, 0.3f, 0.2f, -0.1f},
  //       },
  //       {
  //           _mm_set_ps(0.5f, 0.6f, -0.7f, -0.8f),
  //           {-0.8f, -0.7f, 0.6f, 0.5f},
  //       },
  //       {
  //           _mm_set_ps(0.9f, 0.10f, -0.11f, 0.12f),
  //           {0.12f, -0.11f, 0.10f, 0.9f},
  //       },
  //       {
  //           _mm_set_ps(-1.1f, -2.1f, -3.1f, -4.1f),
  //           {-4.1f, -3.1f, -2.1f, -1.1f},
  //       },
  //       {
  //           _mm_set_ps(100.0f, -110.0f, 120.0f, -130.0f),
  //           {-130.0f, 120.0f, -110.0f, 100.0f},
  //       },
  //       {
  //           _mm_set_ps(200.5f, 210.5f, -220.5f, 230.5f),
  //           {995.74f, -93.04f, 144.03f, 902.50f},
  //       },
  //       {
  //           _mm_set_ps(10.11f, -11.12f, -12.13f, 13.14f),
  //           {13.14f, -12.13f, -11.12f, 10.11f},
  //       },
  //       {
  //           _mm_set_ps(10.1f, -20.2f, 30.3f, 40.4f),
  //           {40.4f, 30.3f, -20.2f, 10.1f},
  //       },
  //   };
  //
  //   for (size_t i = 0; i < (sizeof(test_vec) / (sizeof(test_vec[0]))); i++) {
  //     _mm_prefetch(((const char *)&test_vec[i].a), _MM_HINT_T0);
  //     _mm_prefetch(((const char *)&test_vec[i].a), _MM_HINT_T1);
  //     _mm_prefetch(((const char *)&test_vec[i].a), _MM_HINT_T2);
  //     _mm_prefetch(((const char *)&test_vec[i].a), _MM_HINT_NTA);
  //   }
  //
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_m_psadbw(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  // const uint8_t *_a = (const uint8_t *)impl.test_cases_int_pointer1;
  // const uint8_t *_b = (const uint8_t *)impl.test_cases_int_pointer2;
  // uint16_t d = 0;
  // for (int i = 0; i < 8; i++) {
  //   d += abs(_a[i] - _b[i]);
  // }

  // __m64 a = load_m64(_a);
  // __m64 b = load_m64(_b);
  // __m64 c = _m_psadbw(a, b);
  // return validate_uint16(c, d, 0, 0, 0);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_m_pshufw(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  // return test_mm_shuffle_pi16(impl, iter);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_rcp_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   float dx = 1.0f / _a[0];
  //   float dy = 1.0f / _a[1];
  //   float dz = 1.0f / _a[2];
  //   float dw = 1.0f / _a[3];
  //
  //   __m128 a = load_m128(_a);
  //   __m128 c = _mm_rcp_ps(a);
  //   return validate_float_error(c, dx, dy, dz, dw, 0.001f);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_rcp_ss(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //
  //   float dx = 1.0f / _a[0];
  //   float dy = _a[1];
  //   float dz = _a[2];
  //   float dw = _a[3];
  //   __m128 a = load_m128(_a);
  //   __m128 c = _mm_rcp_ss(a);
  //   return validate_float_error(c, dx, dy, dz, dw, 0.001f);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_rsqrt_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = (const float *)impl.test_cases_float_pointer1;
  //
  //   float f0 = 1 / sqrt(_a[0]);
  //   float f1 = 1 / sqrt(_a[1]);
  //   float f2 = 1 / sqrt(_a[2]);
  //   float f3 = 1 / sqrt(_a[3]);
  //
  //   __m128 a = load_m128(_a);
  //   __m128 c = _mm_rsqrt_ps(a);
  //
  // Here, we ensure the error rate of "_mm_rsqrt_ps()" is under 0.1% compared
  // to the C implementation.
  //   return validate_float_error(c, f0, f1, f2, f3, 0.001f);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_rsqrt_ss(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = (const float *)impl.test_cases_float_pointer1;
  //
  //   float f0 = 1 / sqrt(_a[0]);
  //   float f1 = _a[1];
  //   float f2 = _a[2];
  //   float f3 = _a[3];
  //
  //   __m128 a = load_m128(_a);
  //   __m128 c = _mm_rsqrt_ss(a);
  //
  // Here, we ensure the error rate of "_mm_rsqrt_ps()" is under 0.1% compared
  // to the C implementation.
  //   return validate_float_error(c, f0, f1, f2, f3, 0.001f);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_sad_pu8(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const uint8_t *_a = (const uint8_t *)impl.test_cases_int_pointer1;
  //   const uint8_t *_b = (const uint8_t *)impl.test_cases_int_pointer2;
  //   uint16_t d = 0;
  //   for (int i = 0; i < 8; i++) {
  //     d += abs(_a[i] - _b[i]);
  //   }
  //
  //   __m64 a = load_m64(_a);
  //   __m64 b = load_m64(_b);
  //   __m64 c = _mm_sad_pu8(a, b);
  //   return validate_uint16(c, d, 0, 0, 0);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_set_flush_zero_mode(const SSE2RVV_TEST_IMPL &impl,
                                     // #ifdef ENABLE_TEST_ALL
                                     uint32_t iter) {
  // TODO:
  // After the behavior of denormal number and flush zero mode is fully
  // investigated, the testing would be added.
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_set_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  float e0 = impl.test_cases_floats[iter];
  float e1 = impl.test_cases_floats[iter + 1];
  float e2 = impl.test_cases_floats[iter + 2];
  float e3 = impl.test_cases_floats[iter + 3];
  __m128 a = _mm_set_ps(e3, e2, e1, e0);
  return validate_float(a, e0, e1, e2, e3);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_set_ps1(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  float a = impl.test_cases_floats[iter];

  __m128 ret = _mm_set_ps1(a);

  return validate_float(ret, a, a, a, a);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_set_rounding_mode(const SSE2RVV_TEST_IMPL &impl,
                                   // #ifdef ENABLE_TEST_ALL
                                   uint32_t iter) {
  //   const float *_a = impl.test_cases_float_pointer1;
  //   result_t res_toward_zero, res_to_neg_inf, res_to_pos_inf, res_nearest;
  //
  //   __m128 a = load_m128(_a);
  //   __m128 b, c;
  //
  //   _MM_SET_ROUNDING_MODE(_MM_ROUND_TOWARD_ZERO);
  //   b = _mm_round_ps(a, _MM_FROUND_CUR_DIRECTION);
  //   c = _mm_round_ps(a, _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC);
  //   res_toward_zero = validate_128bits(c, b);
  //
  //   _MM_SET_ROUNDING_MODE(_MM_ROUND_DOWN);
  //   b = _mm_round_ps(a, _MM_FROUND_CUR_DIRECTION);
  //   c = _mm_round_ps(a, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
  //   res_to_neg_inf = validate_128bits(c, b);
  //
  //   _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
  //   b = _mm_round_ps(a, _MM_FROUND_CUR_DIRECTION);
  //   c = _mm_round_ps(a, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
  //   res_to_pos_inf = validate_128bits(c, b);
  //
  //   _MM_SET_ROUNDING_MODE(_MM_ROUND_NEAREST);
  //   b = _mm_round_ps(a, _MM_FROUND_CUR_DIRECTION);
  //   c = _mm_round_ps(a, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
  //   res_nearest = validate_128bits(c, b);
  //
  //   if (res_toward_zero == TEST_SUCCESS && res_to_neg_inf == TEST_SUCCESS &&
  //       res_to_pos_inf == TEST_SUCCESS && res_nearest == TEST_SUCCESS) {
  //     return TEST_SUCCESS;
  //   } else {
  //     return TEST_FAIL;
  //   }
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_set_ss(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  float a = impl.test_cases_floats[iter];
  __m128 c = _mm_set_ss(a);
  return validate_float(c, a, 0, 0, 0);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_set1_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  float w = impl.test_cases_floats[iter];
  __m128 a = _mm_set1_ps(w);
  return validate_float(a, w, w, w, w);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_setcsr(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   return test_mm_set_rounding_mode(impl, iter);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_setr_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  float x = impl.test_cases_floats[iter];
  float y = impl.test_cases_floats[iter + 1];
  float z = impl.test_cases_floats[iter + 2];
  float w = impl.test_cases_floats[iter + 3];

  __m128 ret = _mm_setr_ps(w, z, y, x);

  return validate_float(ret, w, z, y, x);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_setzero_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   __m128 a = _mm_setzero_ps();
  //   return validate_float(a, 0, 0, 0, 0);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_sfence(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   /* FIXME: Assume that memory barriers always function as intended. */
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_shuffle_pi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  //   __m64 a;
  //   __m64 d;
  //
  // #define TEST_IMPL(IDX)
  //   a = load_m64(_a);
  //   d = _mm_shuffle_pi16(a, IDX);
  //
  //   int16_t _d##IDX[4];
  //   _d##IDX[0] = _a[IDX & 0x3];
  //   _d##IDX[1] = _a[(IDX >> 2) & 0x3];
  //   _d##IDX[2] = _a[(IDX >> 4) & 0x3];
  //   _d##IDX[3] = _a[(IDX >> 6) & 0x3];
  //   if (VALIDATE_INT16_M64(d, _d##IDX) != TEST_SUCCESS) {
  //     return TEST_FAIL;
  //   }
  //
  //   IMM_256_ITER
  // #undef TEST_IMPL
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

// Note, NEON does not have a general purpose shuffled command like SSE.
// When invoking this method, there is special code for a number of the most
// common shuffle permutations
result_t test_mm_shuffle_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   const float *_b = impl.test_cases_float_pointer2;
  //   result_t isValid = TEST_SUCCESS;
  //   __m128 a = load_m128(_a);
  //   __m128 b = load_m128(_b);
  // Test many permutations of the shuffle operation, including all
  // permutations which have an optimized/customized implementation
  //   __m128 ret;
  //   ret = _mm_shuffle_ps(a, b, _MM_SHUFFLE(0, 1, 2, 3));
  //   if (!validate_float(ret, _a[3], _a[2], _b[1], _b[0])) {
  //     isValid = TEST_FAIL;
  //   }
  //   ret = _mm_shuffle_ps(a, b, _MM_SHUFFLE(3, 2, 1, 0));
  //   if (!validate_float(ret, _a[0], _a[1], _b[2], _b[3])) {
  //     isValid = TEST_FAIL;
  //   }
  //   ret = _mm_shuffle_ps(a, b, _MM_SHUFFLE(0, 0, 1, 1));
  //   if (!validate_float(ret, _a[1], _a[1], _b[0], _b[0])) {
  //     isValid = TEST_FAIL;
  //   }
  //   ret = _mm_shuffle_ps(a, b, _MM_SHUFFLE(3, 1, 0, 2));
  //   if (!validate_float(ret, _a[2], _a[0], _b[1], _b[3])) {
  //     isValid = TEST_FAIL;
  //   }
  //   ret = _mm_shuffle_ps(a, b, _MM_SHUFFLE(1, 0, 3, 2));
  //   if (!validate_float(ret, _a[2], _a[3], _b[0], _b[1])) {
  //     isValid = TEST_FAIL;
  //   }
  //   ret = _mm_shuffle_ps(a, b, _MM_SHUFFLE(2, 3, 0, 1));
  //   if (!validate_float(ret, _a[1], _a[0], _b[3], _b[2])) {
  //     isValid = TEST_FAIL;
  //   }
  //   ret = _mm_shuffle_ps(a, b, _MM_SHUFFLE(0, 0, 2, 2));
  //   if (!validate_float(ret, _a[2], _a[2], _b[0], _b[0])) {
  //     isValid = TEST_FAIL;
  //   }
  //   ret = _mm_shuffle_ps(a, b, _MM_SHUFFLE(2, 2, 0, 0));
  //   if (!validate_float(ret, _a[0], _a[0], _b[2], _b[2])) {
  //     isValid = TEST_FAIL;
  //   }
  //   ret = _mm_shuffle_ps(a, b, _MM_SHUFFLE(3, 2, 0, 2));
  //   if (!validate_float(ret, _a[2], _a[0], _b[2], _b[3])) {
  //     isValid = TEST_FAIL;
  //   }
  //   ret = _mm_shuffle_ps(a, b, _MM_SHUFFLE(1, 1, 3, 3));
  //   if (!validate_float(ret, _a[3], _a[3], _b[1], _b[1])) {
  //     isValid = TEST_FAIL;
  //   }
  //   ret = _mm_shuffle_ps(a, b, _MM_SHUFFLE(2, 0, 1, 0));
  //   if (!validate_float(ret, _a[0], _a[1], _b[0], _b[2])) {
  //     isValid = TEST_FAIL;
  //   }
  //   ret = _mm_shuffle_ps(a, b, _MM_SHUFFLE(2, 0, 0, 1));
  //   if (!validate_float(ret, _a[1], _a[0], _b[0], _b[2])) {
  //     isValid = TEST_FAIL;
  //   }
  //   ret = _mm_shuffle_ps(a, b, _MM_SHUFFLE(2, 0, 3, 2));
  //   if (!validate_float(ret, _a[2], _a[3], _b[0], _b[2])) {
  //     isValid = TEST_FAIL;
  //   }
  //
  //   return isValid;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_sqrt_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = (const float *)impl.test_cases_float_pointer1;
  //
  //   float f0 = sqrt(_a[0]);
  //   float f1 = sqrt(_a[1]);
  //   float f2 = sqrt(_a[2]);
  //   float f3 = sqrt(_a[3]);
  //
  //   __m128 a = load_m128(_a);
  //   __m128 c = _mm_sqrt_ps(a);
  //
  //   return validate_float_error(c, f0, f1, f2, f3, 0.000001f);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_sqrt_ss(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = (const float *)impl.test_cases_float_pointer1;
  //
  //   float f0 = sqrt(_a[0]);
  //   float f1 = _a[1];
  //   float f2 = _a[2];
  //   float f3 = _a[3];
  //
  //   __m128 a = load_m128(_a);
  //   __m128 c = _mm_sqrt_ss(a);
  //
  //   return validate_float_error(c, f0, f1, f2, f3, 0.000001f);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_store_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  int32_t *p = impl.test_cases_int_pointer1;
  int32_t x = impl.test_cases_ints[iter];
  int32_t y = impl.test_cases_ints[iter + 1];
  int32_t z = impl.test_cases_ints[iter + 2];
  int32_t w = impl.test_cases_ints[iter + 3];
  __m128i a = _mm_set_epi32(x, y, z, w);
  _mm_store_ps((float *)p, *(const __m128 *)&a);
  ASSERT_RETURN(p[0] == w);
  ASSERT_RETURN(p[1] == z);
  ASSERT_RETURN(p[2] == y);
  ASSERT_RETURN(p[3] == x);
  return TEST_SUCCESS;
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_store_ps1(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  float *p = impl.test_cases_float_pointer1;
  float d[4];

  __m128 a = load_m128(p);
  _mm_store_ps1(d, a);

  ASSERT_RETURN(d[0] == *p);
  ASSERT_RETURN(d[1] == *p);
  ASSERT_RETURN(d[2] == *p);
  ASSERT_RETURN(d[3] == *p);
  return TEST_SUCCESS;
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_store_ss(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  float x = impl.test_cases_floats[iter];
  float p[4];

  __m128 a = _mm_set_ss(x);
  _mm_store_ss(p, a);
  ASSERT_RETURN(p[0] == x);
  return TEST_SUCCESS;
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_store1_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  float *p = impl.test_cases_float_pointer1;
  float d[4];

  __m128 a = load_m128(p);
  _mm_store1_ps(d, a);

  ASSERT_RETURN(d[0] == *p);
  ASSERT_RETURN(d[1] == *p);
  ASSERT_RETURN(d[2] == *p);
  ASSERT_RETURN(d[3] == *p);
  return TEST_SUCCESS;
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_storeh_pi(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *p = impl.test_cases_float_pointer1;
  //   float d[4] = {1.0f, 2.0f, 3.0f, 4.0f};
  //   __m128 a = _mm_load_ps(p);
  //   __m64 *b = (__m64 *)d;
  //
  //   _mm_storeh_pi(b, a);
  //   ASSERT_RETURN(d[0] == p[2]);
  //   ASSERT_RETURN(d[1] == p[3]);
  //   ASSERT_RETURN(d[2] == 3.0f);
  //   ASSERT_RETURN(d[3] == 4.0f);
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_storel_pi(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *p = impl.test_cases_float_pointer1;
  //   float d[4] = {1.0f, 2.0f, 3.0f, 4.0f};
  //   __m128 a = _mm_load_ps(p);
  //   __m64 *b = (__m64 *)d;
  //
  //   _mm_storel_pi(b, a);
  //   ASSERT_RETURN(d[0] == p[0]);
  //   ASSERT_RETURN(d[1] == p[1]);
  //   ASSERT_RETURN(d[2] == 3.0f);
  //   ASSERT_RETURN(d[3] == 4.0f);
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_storer_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   float *p = impl.test_cases_float_pointer1;
  //   float d[4];
  //
  //   __m128 a = load_m128(p);
  //   _mm_storer_ps(d, a);
  //
  //   ASSERT_RETURN(d[0] == p[3]);
  //   ASSERT_RETURN(d[1] == p[2]);
  //   ASSERT_RETURN(d[2] == p[1]);
  //   ASSERT_RETURN(d[3] == p[0]);
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_storeu_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   float *_a = impl.test_cases_float_pointer1;
  //   float f[4];
  //   __m128 a = _mm_load_ps(_a);
  //
  //   _mm_storeu_ps(f, a);
  //   return validate_float(a, f[0], f[1], f[2], f[3]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_storeu_si16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  // The GCC version before 11 does not implement intrinsic function
  // _mm_storeu_si16. Check https://gcc.gnu.org/bugzilla/show_bug.cgi?id=95483
  // for more information.
  // #if (defined(__GNUC__) && !defined(__clang__)) && (__GNUC__ <= 10)
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
  // #else
  //   const int32_t *_a = (const int32_t *)impl.test_cases_int_pointer1;
  //   __m128i b;
  //   __m128i a = load_m128i(_a);
  //   _mm_storeu_si16(&b, a);
  //   int16_t *_b = (int16_t *)&b;
  //   int16_t *_c = (int16_t *)&a;
  //   return validate_int16(b, _c[0], _b[1], _b[2], _b[3], _b[4], _b[5], _b[6],
  //                         _b[7]);
  // #endif
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_storeu_si64(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  // Versions of GCC prior to 9 do not implement intrinsic function
  // _mm_storeu_si64. Check https://gcc.gnu.org/bugzilla/show_bug.cgi?id=87558
  // for more information.
  // #if (defined(__GNUC__) && !defined(__clang__)) && (__GNUC__ < 9)
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
  // #else
  //   const int32_t *_a = (const int32_t *)impl.test_cases_int_pointer1;
  //   __m128i b;
  //   __m128i a = load_m128i(_a);
  //   _mm_storeu_si64(&b, a);
  //   int64_t *_b = (int64_t *)&b;
  //   int64_t *_c = (int64_t *)&a;
  //   return validate_int64(b, _c[0], _b[1]);
  // #endif
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_stream_pi(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int64_t *_a = (const int64_t *)impl.test_cases_int_pointer1;
  //   __m64 a = load_m64(_a);
  //   __m64 p;
  //
  //   _mm_stream_pi(&p, a);
  //   return validate_int64(p, _a[0]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_stream_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   __m128 a = load_m128(_a);
  //   alignas(16) float p[4];
  //
  //   _mm_stream_ps(p, a);
  //   ASSERT_RETURN(p[0] == _a[0]);
  //   ASSERT_RETURN(p[1] == _a[1]);
  //   ASSERT_RETURN(p[2] == _a[2]);
  //   ASSERT_RETURN(p[3] == _a[3]);
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_sub_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const float *_a = impl.test_cases_float_pointer1;
  const float *_b = impl.test_cases_float_pointer2;
  float dx = _a[0] - _b[0];
  float dy = _a[1] - _b[1];
  float dz = _a[2] - _b[2];
  float dw = _a[3] - _b[3];

  __m128 a = load_m128(_a);
  __m128 b = load_m128(_b);
  __m128 c = _mm_sub_ps(a, b);
  return validate_float(c, dx, dy, dz, dw);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_sub_ss(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const float *_a = impl.test_cases_float_pointer1;
  const float *_b = impl.test_cases_float_pointer2;
  float dx = _a[0] - _b[0];
  float dy = _a[1];
  float dz = _a[2];
  float dw = _a[3];

  __m128 a = load_m128(_a);
  __m128 b = load_m128(_b);
  __m128 c = _mm_sub_ss(a, b);
  return validate_float(c, dx, dy, dz, dw);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_ucomieq_ss(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
#if defined(__GNUC__) && !defined(__clang__)
  return TEST_UNIMPL;
#else
  const float *_a = impl.test_cases_float_pointer1;
  const float *_b = impl.test_cases_float_pointer2;

  int32_t _c = comieq_ss(_a[0], _b[0]);
  __m128 a = load_m128(_a);
  __m128 b = load_m128(_b);
  int32_t c = _mm_comieq_ss(a, b);

  return _c == c ? TEST_SUCCESS : TEST_FAIL;
#endif
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_ucomige_ss(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  // _mm_ucomige_ss is equal to _mm_comige_ss
  //   return test_mm_comige_ss(impl, iter);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_ucomigt_ss(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  // _mm_ucomigt_ss is equal to _mm_comigt_ss
  //   return test_mm_comigt_ss(impl, iter);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_ucomile_ss(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  // _mm_ucomile_ss is equal to _mm_comile_ss
  //   return test_mm_comile_ss(impl, iter);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_ucomilt_ss(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  // _mm_ucomilt_ss is equal to _mm_comilt_ss
  //   return test_mm_comilt_ss(impl, iter);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_ucomineq_ss(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  // _mm_ucomineq_ss is equal to _mm_comineq_ss
  //   return test_mm_comineq_ss(impl, iter);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_undefined_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   __m128 a = _mm_undefined_ps();
  //   a = _mm_xor_ps(a, a);
  //   return validate_float(a, 0, 0, 0, 0);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_unpackhi_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   float *_a = impl.test_cases_float_pointer1;
  //   float *_b = impl.test_cases_float_pointer1;
  //
  //   float f0 = _a[2];
  //   float f1 = _b[2];
  //   float f2 = _a[3];
  //   float f3 = _b[3];
  //
  //   __m128 a = _mm_load_ps(_a);
  //   __m128 b = _mm_load_ps(_b);
  //   __m128 c = _mm_unpackhi_ps(a, b);
  //   return validate_float(c, f0, f1, f2, f3);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_unpacklo_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   float *_a = impl.test_cases_float_pointer1;
  //   float *_b = impl.test_cases_float_pointer1;
  //
  //   float f0 = _a[0];
  //   float f1 = _b[0];
  //   float f2 = _a[1];
  //   float f3 = _b[1];
  //
  //   __m128 a = _mm_load_ps(_a);
  //   __m128 b = _mm_load_ps(_b);
  //   __m128 c = _mm_unpacklo_ps(a, b);
  //
  //   return validate_float(c, f0, f1, f2, f3);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_xor_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int32_t *_a = (const int32_t *)impl.test_cases_float_pointer1;
  //   const int32_t *_b = (const int32_t *)impl.test_cases_float_pointer2;
  //
  //   int32_t d0 = _a[0] ^ _b[0];
  //   int32_t d1 = _a[1] ^ _b[1];
  //   int32_t d2 = _a[2] ^ _b[2];
  //   int32_t d3 = _a[3] ^ _b[3];
  //
  //   __m128 a = load_m128(_a);
  //   __m128 b = load_m128(_b);
  //   __m128 c = _mm_xor_ps(a, b);
  //
  //   return validate_float(c, *((float *)&d0), *((float *)&d1), *((float
  //   *)&d2),
  //                         *((float *)&d3));
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

/* SSE2 */
result_t test_mm_add_epi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  const int16_t *_b = (const int16_t *)impl.test_cases_int_pointer2;

  int16_t d[8];
  d[0] = _a[0] + _b[0];
  d[1] = _a[1] + _b[1];
  d[2] = _a[2] + _b[2];
  d[3] = _a[3] + _b[3];
  d[4] = _a[4] + _b[4];
  d[5] = _a[5] + _b[5];
  d[6] = _a[6] + _b[6];
  d[7] = _a[7] + _b[7];

  __m128i a = load_m128i(_a);
  __m128i b = load_m128i(_b);
  __m128i c = _mm_add_epi16(a, b);

  return VALIDATE_INT16_M128(c, d);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_add_epi32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int32_t *_a = impl.test_cases_int_pointer1;
  const int32_t *_b = impl.test_cases_int_pointer2;
  int32_t d[4];
  d[0] = _a[0] + _b[0];
  d[1] = _a[1] + _b[1];
  d[2] = _a[2] + _b[2];
  d[3] = _a[3] + _b[3];

  __m128i a = load_m128i(_a);
  __m128i b = load_m128i(_b);
  __m128i c = _mm_add_epi32(a, b);
  return VALIDATE_INT32_M128(c, d);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_add_epi64(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int64_t *_a = (const int64_t *)impl.test_cases_int_pointer1;
  const int64_t *_b = (const int64_t *)impl.test_cases_int_pointer2;

  int64_t d0 = _a[0] + _b[0];
  int64_t d1 = _a[1] + _b[1];

  __m128i a = load_m128i(_a);
  __m128i b = load_m128i(_b);
  __m128i c = _mm_add_epi64(a, b);

  return validate_int64(c, d0, d1);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_add_epi8(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int8_t *_a = (const int8_t *)impl.test_cases_int_pointer1;
  const int8_t *_b = (const int8_t *)impl.test_cases_int_pointer2;
  int8_t d[16];
  d[0] = _a[0] + _b[0];
  d[1] = _a[1] + _b[1];
  d[2] = _a[2] + _b[2];
  d[3] = _a[3] + _b[3];
  d[4] = _a[4] + _b[4];
  d[5] = _a[5] + _b[5];
  d[6] = _a[6] + _b[6];
  d[7] = _a[7] + _b[7];
  d[8] = _a[8] + _b[8];
  d[9] = _a[9] + _b[9];
  d[10] = _a[10] + _b[10];
  d[11] = _a[11] + _b[11];
  d[12] = _a[12] + _b[12];
  d[13] = _a[13] + _b[13];
  d[14] = _a[14] + _b[14];
  d[15] = _a[15] + _b[15];

  __m128i a = load_m128i(_a);
  __m128i b = load_m128i(_b);
  __m128i c = _mm_add_epi8(a, b);
  return VALIDATE_INT8_M128(c, d);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_add_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const double *_a = (const double *)impl.test_cases_float_pointer1;
  const double *_b = (const double *)impl.test_cases_float_pointer2;
  double d0 = _a[0] + _b[0];
  double d1 = _a[1] + _b[1];

  __m128d a = load_m128d(_a);
  __m128d b = load_m128d(_b);
  __m128d c = _mm_add_pd(a, b);
  return validate_double(c, d0, d1);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_add_sd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const double *_a = (const double *)impl.test_cases_float_pointer1;
  const double *_b = (const double *)impl.test_cases_float_pointer2;
  double d0 = _a[0] + _b[0];
  double d1 = _a[1];

  __m128d a = load_m128d(_a);
  __m128d b = load_m128d(_b);
  __m128d c = _mm_add_sd(a, b);
  return validate_double(c, d0, d1);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_add_si64(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int64_t *_a = (const int64_t *)impl.test_cases_int_pointer1;
  const int64_t *_b = (const int64_t *)impl.test_cases_int_pointer2;

  int64_t d0 = _a[0] + _b[0];

  __m64 a = load_m64(_a);
  __m64 b = load_m64(_b);
  __m64 c = _mm_add_si64(a, b);

  return validate_int64(c, d0);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_adds_epi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  const int16_t *_b = (const int16_t *)impl.test_cases_int_pointer2;
  int32_t d[8];
  d[0] = (int32_t)_a[0] + (int32_t)_b[0];
  if (d[0] > 32767)
    d[0] = 32767;
  if (d[0] < -32768)
    d[0] = -32768;
  d[1] = (int32_t)_a[1] + (int32_t)_b[1];
  if (d[1] > 32767)
    d[1] = 32767;
  if (d[1] < -32768)
    d[1] = -32768;
  d[2] = (int32_t)_a[2] + (int32_t)_b[2];
  if (d[2] > 32767)
    d[2] = 32767;
  if (d[2] < -32768)
    d[2] = -32768;
  d[3] = (int32_t)_a[3] + (int32_t)_b[3];
  if (d[3] > 32767)
    d[3] = 32767;
  if (d[3] < -32768)
    d[3] = -32768;
  d[4] = (int32_t)_a[4] + (int32_t)_b[4];
  if (d[4] > 32767)
    d[4] = 32767;
  if (d[4] < -32768)
    d[4] = -32768;
  d[5] = (int32_t)_a[5] + (int32_t)_b[5];
  if (d[5] > 32767)
    d[5] = 32767;
  if (d[5] < -32768)
    d[5] = -32768;
  d[6] = (int32_t)_a[6] + (int32_t)_b[6];
  if (d[6] > 32767)
    d[6] = 32767;
  if (d[6] < -32768)
    d[6] = -32768;
  d[7] = (int32_t)_a[7] + (int32_t)_b[7];
  if (d[7] > 32767)
    d[7] = 32767;
  if (d[7] < -32768)
    d[7] = -32768;

  __m128i a = load_m128i(_a);
  __m128i b = load_m128i(_b);

  __m128i c = _mm_adds_epi16(a, b);
  return VALIDATE_INT16_M128(c, d);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_adds_epi8(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int8_t *_a = (const int8_t *)impl.test_cases_int_pointer1;
  const int8_t *_b = (const int8_t *)impl.test_cases_int_pointer2;

  int16_t d[16];
  for (int i = 0; i < 16; i++) {
    d[i] = (int16_t)_a[i] + (int16_t)_b[i];
    if (d[i] > 127)
      d[i] = 127;
    if (d[i] < -128)
      d[i] = -128;
  }

  __m128i a = load_m128i(_a);
  __m128i b = load_m128i(_b);
  __m128i c = _mm_adds_epi8(a, b);

  return VALIDATE_INT8_M128(c, d);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_adds_epu16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  uint32_t max = 0xFFFF;
  const uint16_t *_a = (const uint16_t *)impl.test_cases_int_pointer1;
  const uint16_t *_b = (const uint16_t *)impl.test_cases_int_pointer2;

  uint16_t d[8];
  d[0] = (uint32_t)_a[0] + (uint32_t)_b[0] > max ? max : _a[0] + _b[0];
  d[1] = (uint32_t)_a[1] + (uint32_t)_b[1] > max ? max : _a[1] + _b[1];
  d[2] = (uint32_t)_a[2] + (uint32_t)_b[2] > max ? max : _a[2] + _b[2];
  d[3] = (uint32_t)_a[3] + (uint32_t)_b[3] > max ? max : _a[3] + _b[3];
  d[4] = (uint32_t)_a[4] + (uint32_t)_b[4] > max ? max : _a[4] + _b[4];
  d[5] = (uint32_t)_a[5] + (uint32_t)_b[5] > max ? max : _a[5] + _b[5];
  d[6] = (uint32_t)_a[6] + (uint32_t)_b[6] > max ? max : _a[6] + _b[6];
  d[7] = (uint32_t)_a[7] + (uint32_t)_b[7] > max ? max : _a[7] + _b[7];

  __m128i a = load_m128i(_a);
  __m128i b = load_m128i(_b);
  __m128i c = _mm_adds_epu16(a, b);

  return VALIDATE_INT16_M128(c, d);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_adds_epu8(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int8_t *_a = (const int8_t *)impl.test_cases_int_pointer1;
  const int8_t *_b = (const int8_t *)impl.test_cases_int_pointer2;
  uint8_t d[16];
  d[0] = (uint8_t)_a[0] + (uint8_t)_b[0];
  if (d[0] < (uint8_t)_a[0])
    d[0] = 255;
  d[1] = (uint8_t)_a[1] + (uint8_t)_b[1];
  if (d[1] < (uint8_t)_a[1])
    d[1] = 255;
  d[2] = (uint8_t)_a[2] + (uint8_t)_b[2];
  if (d[2] < (uint8_t)_a[2])
    d[2] = 255;
  d[3] = (uint8_t)_a[3] + (uint8_t)_b[3];
  if (d[3] < (uint8_t)_a[3])
    d[3] = 255;
  d[4] = (uint8_t)_a[4] + (uint8_t)_b[4];
  if (d[4] < (uint8_t)_a[4])
    d[4] = 255;
  d[5] = (uint8_t)_a[5] + (uint8_t)_b[5];
  if (d[5] < (uint8_t)_a[5])
    d[5] = 255;
  d[6] = (uint8_t)_a[6] + (uint8_t)_b[6];
  if (d[6] < (uint8_t)_a[6])
    d[6] = 255;
  d[7] = (uint8_t)_a[7] + (uint8_t)_b[7];
  if (d[7] < (uint8_t)_a[7])
    d[7] = 255;
  d[8] = (uint8_t)_a[8] + (uint8_t)_b[8];
  if (d[8] < (uint8_t)_a[8])
    d[8] = 255;
  d[9] = (uint8_t)_a[9] + (uint8_t)_b[9];
  if (d[9] < (uint8_t)_a[9])
    d[9] = 255;
  d[10] = (uint8_t)_a[10] + (uint8_t)_b[10];
  if (d[10] < (uint8_t)_a[10])
    d[10] = 255;
  d[11] = (uint8_t)_a[11] + (uint8_t)_b[11];
  if (d[11] < (uint8_t)_a[11])
    d[11] = 255;
  d[12] = (uint8_t)_a[12] + (uint8_t)_b[12];
  if (d[12] < (uint8_t)_a[12])
    d[12] = 255;
  d[13] = (uint8_t)_a[13] + (uint8_t)_b[13];
  if (d[13] < (uint8_t)_a[13])
    d[13] = 255;
  d[14] = (uint8_t)_a[14] + (uint8_t)_b[14];
  if (d[14] < (uint8_t)_a[14])
    d[14] = 255;
  d[15] = (uint8_t)_a[15] + (uint8_t)_b[15];
  if (d[15] < (uint8_t)_a[15])
    d[15] = 255;

  __m128i a = load_m128i(_a);
  __m128i b = load_m128i(_b);
  __m128i c = _mm_adds_epu8(a, b);
  return VALIDATE_INT8_M128(c, d);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_and_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int64_t *_a = (const int64_t *)impl.test_cases_float_pointer1;
  const int64_t *_b = (const int64_t *)impl.test_cases_float_pointer2;

  int64_t d0 = _a[0] & _b[0];
  int64_t d1 = _a[1] & _b[1];

  __m128d a = load_m128d(_a);
  __m128d b = load_m128d(_b);
  __m128d c = _mm_and_pd(a, b);

  return validate_double(c, *((double *)&d0), *((double *)&d1));
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_and_si128(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int32_t *_a = impl.test_cases_int_pointer1;
  const int32_t *_b = impl.test_cases_int_pointer2;
  __m128i a = load_m128i(_a);
  __m128i b = load_m128i(_b);
  __m128 fc = _mm_and_ps(*(const __m128 *)&a, *(const __m128 *)&b);
  __m128i c = *(const __m128i *)&fc;
  // now for the assertion...
  const uint32_t *ia = (const uint32_t *)&a;
  const uint32_t *ib = (const uint32_t *)&b;
  uint32_t r[4];
  r[0] = ia[0] & ib[0];
  r[1] = ia[1] & ib[1];
  r[2] = ia[2] & ib[2];
  r[3] = ia[3] & ib[3];
  __m128i ret = do_mm_set_epi32(r[3], r[2], r[1], r[0]);
  result_t res = VALIDATE_INT32_M128(c, r);
  if (res) {
    res = VALIDATE_INT32_M128(ret, r);
  }
  return res;
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_andnot_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const double *_a = (const double *)impl.test_cases_float_pointer1;
  const double *_b = (const double *)impl.test_cases_float_pointer2;

  __m128d a = load_m128d(_a);
  __m128d b = load_m128d(_b);
  __m128d c = _mm_andnot_pd(a, b);

  // Take AND operation a complement of 'a' and 'b'. Bitwise operations are
  // not allowed on float/double datatype, so 'a' and 'b' are calculated in
  // uint64_t datatype.
  const uint64_t *ia = (const uint64_t *)&a;
  const uint64_t *ib = (const uint64_t *)&b;
  uint64_t r0 = ~ia[0] & ib[0];
  uint64_t r1 = ~ia[1] & ib[1];
  return validate_uint64(*(const __m128i *)&c, r0, r1);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_andnot_si128(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int32_t *_a = impl.test_cases_int_pointer1;
  const int32_t *_b = impl.test_cases_int_pointer2;
  __m128i a = load_m128i(_a);
  __m128i b = load_m128i(_b);
  __m128 fc = _mm_andnot_ps(*(const __m128 *)&a, *(const __m128 *)&b);
  __m128i c = *(const __m128i *)&fc;
  // now for the assertion...
  const uint32_t *ia = (const uint32_t *)&a;
  const uint32_t *ib = (const uint32_t *)&b;
  uint32_t r[4];
  r[0] = ~ia[0] & ib[0];
  r[1] = ~ia[1] & ib[1];
  r[2] = ~ia[2] & ib[2];
  r[3] = ~ia[3] & ib[3];
  __m128i ret = do_mm_set_epi32(r[3], r[2], r[1], r[0]);
  result_t res = TEST_SUCCESS;
  res = VALIDATE_INT32_M128(c, r);
  if (res) {
    res = VALIDATE_INT32_M128(ret, r);
  }
  return res;
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_avg_epu16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  const int16_t *_b = (const int16_t *)impl.test_cases_int_pointer2;
  uint16_t d[8];
  d[0] = ((uint16_t)_a[0] + (uint16_t)_b[0] + 1) >> 1;
  d[1] = ((uint16_t)_a[1] + (uint16_t)_b[1] + 1) >> 1;
  d[2] = ((uint16_t)_a[2] + (uint16_t)_b[2] + 1) >> 1;
  d[3] = ((uint16_t)_a[3] + (uint16_t)_b[3] + 1) >> 1;
  d[4] = ((uint16_t)_a[4] + (uint16_t)_b[4] + 1) >> 1;
  d[5] = ((uint16_t)_a[5] + (uint16_t)_b[5] + 1) >> 1;
  d[6] = ((uint16_t)_a[6] + (uint16_t)_b[6] + 1) >> 1;
  d[7] = ((uint16_t)_a[7] + (uint16_t)_b[7] + 1) >> 1;
  __m128i a = load_m128i(_a);
  __m128i b = load_m128i(_b);
  __m128i c = _mm_avg_epu16(a, b);
  return VALIDATE_UINT16_M128(c, d);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_avg_epu8(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int8_t *_a = (const int8_t *)impl.test_cases_int_pointer1;
  const int8_t *_b = (const int8_t *)impl.test_cases_int_pointer2;
  uint8_t d[16];
  d[0] = ((uint8_t)_a[0] + (uint8_t)_b[0] + 1) >> 1;
  d[1] = ((uint8_t)_a[1] + (uint8_t)_b[1] + 1) >> 1;
  d[2] = ((uint8_t)_a[2] + (uint8_t)_b[2] + 1) >> 1;
  d[3] = ((uint8_t)_a[3] + (uint8_t)_b[3] + 1) >> 1;
  d[4] = ((uint8_t)_a[4] + (uint8_t)_b[4] + 1) >> 1;
  d[5] = ((uint8_t)_a[5] + (uint8_t)_b[5] + 1) >> 1;
  d[6] = ((uint8_t)_a[6] + (uint8_t)_b[6] + 1) >> 1;
  d[7] = ((uint8_t)_a[7] + (uint8_t)_b[7] + 1) >> 1;
  d[8] = ((uint8_t)_a[8] + (uint8_t)_b[8] + 1) >> 1;
  d[9] = ((uint8_t)_a[9] + (uint8_t)_b[9] + 1) >> 1;
  d[10] = ((uint8_t)_a[10] + (uint8_t)_b[10] + 1) >> 1;
  d[11] = ((uint8_t)_a[11] + (uint8_t)_b[11] + 1) >> 1;
  d[12] = ((uint8_t)_a[12] + (uint8_t)_b[12] + 1) >> 1;
  d[13] = ((uint8_t)_a[13] + (uint8_t)_b[13] + 1) >> 1;
  d[14] = ((uint8_t)_a[14] + (uint8_t)_b[14] + 1) >> 1;
  d[15] = ((uint8_t)_a[15] + (uint8_t)_b[15] + 1) >> 1;
  __m128i a = load_m128i(_a);
  __m128i b = load_m128i(_b);
  __m128i c = _mm_avg_epu8(a, b);
  return VALIDATE_UINT8_M128(c, d);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_bslli_si128(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   return test_mm_slli_si128(impl, iter);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_bsrli_si128(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   return test_mm_srli_si128(impl, iter);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_castpd_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const float *_a = impl.test_cases_float_pointer1;
  const __m128d a = load_m128d(_a);
  const __m128 _c = load_m128(_a);

  __m128 r = _mm_castpd_ps(a);

  return validate_128bits(r, _c);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_castpd_si128(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const float *_a = impl.test_cases_float_pointer1;
  const __m128d a = load_m128d(_a);
  const __m128i *_c = (const __m128i *)_a;

  __m128i r = _mm_castpd_si128(a);

  return validate_128bits(r, *_c);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_castps_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const float *_a = impl.test_cases_float_pointer1;
  const __m128 a = load_m128(_a);
  const __m128d *_c = (const __m128d *)_a;

  __m128d r = _mm_castps_pd(a);

  return validate_128bits(r, *_c);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_castps_si128(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const float *_a = impl.test_cases_float_pointer1;

  const __m128i *_c = (const __m128i *)_a;

  const __m128 a = load_m128(_a);
  __m128i r = _mm_castps_si128(a);

  return validate_128bits(r, *_c);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_castsi128_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int32_t *_a = impl.test_cases_int_pointer1;

  const __m128d *_c = (const __m128d *)_a;

  const __m128i a = load_m128i(_a);
  __m128d r = _mm_castsi128_pd(a);

  return validate_128bits(r, *_c);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_castsi128_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int32_t *_a = impl.test_cases_int_pointer1;

  const __m128 *_c = (const __m128 *)_a;

  const __m128i a = load_m128i(_a);
  __m128 r = _mm_castsi128_ps(a);

  return validate_128bits(r, *_c);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_clflush(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   /* FIXME: Assume that we have portable mechanisms to flush cache. */
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cmpeq_epi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  const int16_t *_b = (const int16_t *)impl.test_cases_int_pointer2;
  int16_t d[8];
  d[0] = (_a[0] == _b[0]) ? ~UINT16_C(0) : 0x0;
  d[1] = (_a[1] == _b[1]) ? ~UINT16_C(0) : 0x0;
  d[2] = (_a[2] == _b[2]) ? ~UINT16_C(0) : 0x0;
  d[3] = (_a[3] == _b[3]) ? ~UINT16_C(0) : 0x0;
  d[4] = (_a[4] == _b[4]) ? ~UINT16_C(0) : 0x0;
  d[5] = (_a[5] == _b[5]) ? ~UINT16_C(0) : 0x0;
  d[6] = (_a[6] == _b[6]) ? ~UINT16_C(0) : 0x0;
  d[7] = (_a[7] == _b[7]) ? ~UINT16_C(0) : 0x0;

  __m128i a = load_m128i(_a);
  __m128i b = load_m128i(_b);
  __m128i c = _mm_cmpeq_epi16(a, b);
  return VALIDATE_INT16_M128(c, d);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_cmpeq_epi32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int32_t *_a = impl.test_cases_int_pointer1;
  const int32_t *_b = impl.test_cases_int_pointer2;

  int32_t d[4];
  d[0] = (_a[0] == _b[0]) ? ~UINT32_C(0) : 0x0;
  d[1] = (_a[1] == _b[1]) ? ~UINT32_C(0) : 0x0;
  d[2] = (_a[2] == _b[2]) ? ~UINT32_C(0) : 0x0;
  d[3] = (_a[3] == _b[3]) ? ~UINT32_C(0) : 0x0;

  __m128i a = load_m128i(_a);
  __m128i b = load_m128i(_b);
  __m128i c = _mm_cmpeq_epi32(a, b);

  return VALIDATE_INT32_M128(c, d);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_cmpeq_epi8(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int8_t *_a = (const int8_t *)impl.test_cases_int_pointer1;
  const int8_t *_b = (const int8_t *)impl.test_cases_int_pointer2;
  int8_t d[16];
  d[0] = (_a[0] == _b[0]) ? ~UINT8_C(0) : 0x00;
  d[1] = (_a[1] == _b[1]) ? ~UINT8_C(0) : 0x00;
  d[2] = (_a[2] == _b[2]) ? ~UINT8_C(0) : 0x00;
  d[3] = (_a[3] == _b[3]) ? ~UINT8_C(0) : 0x00;
  d[4] = (_a[4] == _b[4]) ? ~UINT8_C(0) : 0x00;
  d[5] = (_a[5] == _b[5]) ? ~UINT8_C(0) : 0x00;
  d[6] = (_a[6] == _b[6]) ? ~UINT8_C(0) : 0x00;
  d[7] = (_a[7] == _b[7]) ? ~UINT8_C(0) : 0x00;
  d[8] = (_a[8] == _b[8]) ? ~UINT8_C(0) : 0x00;
  d[9] = (_a[9] == _b[9]) ? ~UINT8_C(0) : 0x00;
  d[10] = (_a[10] == _b[10]) ? ~UINT8_C(0) : 0x00;
  d[11] = (_a[11] == _b[11]) ? ~UINT8_C(0) : 0x00;
  d[12] = (_a[12] == _b[12]) ? ~UINT8_C(0) : 0x00;
  d[13] = (_a[13] == _b[13]) ? ~UINT8_C(0) : 0x00;
  d[14] = (_a[14] == _b[14]) ? ~UINT8_C(0) : 0x00;
  d[15] = (_a[15] == _b[15]) ? ~UINT8_C(0) : 0x00;

  __m128i a = load_m128i(_a);
  __m128i b = load_m128i(_b);
  __m128i c = _mm_cmpeq_epi8(a, b);
  return VALIDATE_INT8_M128(c, d);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_cmpeq_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const double *_a = (const double *)impl.test_cases_float_pointer1;
  const double *_b = (const double *)impl.test_cases_float_pointer2;
  uint64_t d0 = (_a[0] == _b[0]) ? 0xffffffffffffffff : 0;
  uint64_t d1 = (_a[1] == _b[1]) ? 0xffffffffffffffff : 0;

  __m128d a = load_m128d(_a);
  __m128d b = load_m128d(_b);
  __m128d c = _mm_cmpeq_pd(a, b);
  return validate_double(c, *(double *)&d0, *(double *)&d1);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_cmpeq_sd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const double *_a = (const double *)impl.test_cases_float_pointer1;
  const double *_b = (const double *)impl.test_cases_float_pointer2;
  const uint64_t d0 = (_a[0] == _b[0]) ? ~UINT64_C(0) : 0;
  const uint64_t d1 = ((const uint64_t *)_a)[1];

  __m128d a = load_m128d(_a);
  __m128d b = load_m128d(_b);
  __m128d c = _mm_cmpeq_sd(a, b);

  return validate_double(c, *(const double *)&d0, *(const double *)&d1);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_cmpge_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const double *_a = (const double *)impl.test_cases_float_pointer1;
  const double *_b = (const double *)impl.test_cases_float_pointer2;
  uint64_t d0 = (_a[0] >= _b[0]) ? ~UINT64_C(0) : 0;
  uint64_t d1 = (_a[1] >= _b[1]) ? ~UINT64_C(0) : 0;

  __m128d a = load_m128d(_a);
  __m128d b = load_m128d(_b);
  __m128d c = _mm_cmpge_pd(a, b);

  return validate_double(c, *(double *)&d0, *(double *)&d1);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_cmpge_sd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  double *_a = (double *)impl.test_cases_float_pointer1;
  double *_b = (double *)impl.test_cases_float_pointer2;
  uint64_t d0 = (_a[0] >= _b[0]) ? ~UINT64_C(0) : 0;
  uint64_t d1 = ((uint64_t *)_a)[1];

  __m128d a = load_m128d(_a);
  __m128d b = load_m128d(_b);
  __m128d c = _mm_cmpge_sd(a, b);

  return validate_double(c, *(double *)&d0, *(double *)&d1);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_cmpgt_epi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  const int16_t *_b = (const int16_t *)impl.test_cases_int_pointer2;
  uint16_t d[8];
  d[0] = _a[0] > _b[0] ? ~UINT16_C(0) : 0;
  d[1] = _a[1] > _b[1] ? ~UINT16_C(0) : 0;
  d[2] = _a[2] > _b[2] ? ~UINT16_C(0) : 0;
  d[3] = _a[3] > _b[3] ? ~UINT16_C(0) : 0;
  d[4] = _a[4] > _b[4] ? ~UINT16_C(0) : 0;
  d[5] = _a[5] > _b[5] ? ~UINT16_C(0) : 0;
  d[6] = _a[6] > _b[6] ? ~UINT16_C(0) : 0;
  d[7] = _a[7] > _b[7] ? ~UINT16_C(0) : 0;

  __m128i a = load_m128i(_a);
  __m128i b = load_m128i(_b);
  __m128i c = _mm_cmpgt_epi16(a, b);

  return VALIDATE_INT16_M128(c, d);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_cmpgt_epi32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int32_t *_a = impl.test_cases_int_pointer1;
  const int32_t *_b = impl.test_cases_int_pointer2;
  __m128i a = load_m128i(_a);
  __m128i b = load_m128i(_b);

  int32_t result[4];

  result[0] = _a[0] > _b[0] ? -1 : 0;
  result[1] = _a[1] > _b[1] ? -1 : 0;
  result[2] = _a[2] > _b[2] ? -1 : 0;
  result[3] = _a[3] > _b[3] ? -1 : 0;

  __m128i iret = _mm_cmpgt_epi32(a, b);
  return VALIDATE_INT32_M128(iret, result);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_cmpgt_epi8(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int8_t *_a = (const int8_t *)impl.test_cases_int_pointer1;
  const int8_t *_b = (const int8_t *)impl.test_cases_int_pointer2;
  int8_t d[16];
  d[0] = (_a[0] > _b[0]) ? ~UINT8_C(0) : 0x00;
  d[1] = (_a[1] > _b[1]) ? ~UINT8_C(0) : 0x00;
  d[2] = (_a[2] > _b[2]) ? ~UINT8_C(0) : 0x00;
  d[3] = (_a[3] > _b[3]) ? ~UINT8_C(0) : 0x00;
  d[4] = (_a[4] > _b[4]) ? ~UINT8_C(0) : 0x00;
  d[5] = (_a[5] > _b[5]) ? ~UINT8_C(0) : 0x00;
  d[6] = (_a[6] > _b[6]) ? ~UINT8_C(0) : 0x00;
  d[7] = (_a[7] > _b[7]) ? ~UINT8_C(0) : 0x00;
  d[8] = (_a[8] > _b[8]) ? ~UINT8_C(0) : 0x00;
  d[9] = (_a[9] > _b[9]) ? ~UINT8_C(0) : 0x00;
  d[10] = (_a[10] > _b[10]) ? ~UINT8_C(0) : 0x00;
  d[11] = (_a[11] > _b[11]) ? ~UINT8_C(0) : 0x00;
  d[12] = (_a[12] > _b[12]) ? ~UINT8_C(0) : 0x00;
  d[13] = (_a[13] > _b[13]) ? ~UINT8_C(0) : 0x00;
  d[14] = (_a[14] > _b[14]) ? ~UINT8_C(0) : 0x00;
  d[15] = (_a[15] > _b[15]) ? ~UINT8_C(0) : 0x00;

  __m128i a = load_m128i(_a);
  __m128i b = load_m128i(_b);
  __m128i c = _mm_cmpgt_epi8(a, b);
  return VALIDATE_INT8_M128(c, d);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_cmpgt_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const double *_a = (const double *)impl.test_cases_float_pointer1;
  const double *_b = (const double *)impl.test_cases_float_pointer2;
  uint64_t d0 = (_a[0] > _b[0]) ? ~UINT64_C(0) : 0;
  uint64_t d1 = (_a[1] > _b[1]) ? ~UINT64_C(0) : 0;

  __m128d a = load_m128d(_a);
  __m128d b = load_m128d(_b);
  __m128d c = _mm_cmpgt_pd(a, b);

  return validate_double(c, *(double *)&d0, *(double *)&d1);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_cmpgt_sd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  double *_a = (double *)impl.test_cases_float_pointer1;
  double *_b = (double *)impl.test_cases_float_pointer2;
  uint64_t d0 = (_a[0] > _b[0]) ? ~UINT64_C(0) : 0;
  uint64_t d1 = ((uint64_t *)_a)[1];

  __m128d a = load_m128d(_a);
  __m128d b = load_m128d(_b);
  __m128d c = _mm_cmpgt_sd(a, b);

  return validate_double(c, *(double *)&d0, *(double *)&d1);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_cmple_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const double *_a = (const double *)impl.test_cases_float_pointer1;
  const double *_b = (const double *)impl.test_cases_float_pointer2;
  uint64_t d0 = (_a[0] <= _b[0]) ? ~UINT64_C(0) : 0;
  uint64_t d1 = (_a[1] <= _b[1]) ? ~UINT64_C(0) : 0;

  __m128d a = load_m128d(_a);
  __m128d b = load_m128d(_b);
  __m128d c = _mm_cmple_pd(a, b);

  return validate_double(c, *(double *)&d0, *(double *)&d1);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_cmple_sd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  double *_a = (double *)impl.test_cases_float_pointer1;
  double *_b = (double *)impl.test_cases_float_pointer2;
  uint64_t d0 = (_a[0] <= _b[0]) ? ~UINT64_C(0) : 0;
  uint64_t d1 = ((uint64_t *)_a)[1];

  __m128d a = load_m128d(_a);
  __m128d b = load_m128d(_b);
  __m128d c = _mm_cmple_sd(a, b);

  return validate_double(c, *(double *)&d0, *(double *)&d1);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_cmplt_epi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  const int16_t *_b = (const int16_t *)impl.test_cases_int_pointer2;
  uint16_t d[8];
  d[0] = _a[0] < _b[0] ? ~UINT16_C(0) : 0;
  d[1] = _a[1] < _b[1] ? ~UINT16_C(0) : 0;
  d[2] = _a[2] < _b[2] ? ~UINT16_C(0) : 0;
  d[3] = _a[3] < _b[3] ? ~UINT16_C(0) : 0;
  d[4] = _a[4] < _b[4] ? ~UINT16_C(0) : 0;
  d[5] = _a[5] < _b[5] ? ~UINT16_C(0) : 0;
  d[6] = _a[6] < _b[6] ? ~UINT16_C(0) : 0;
  d[7] = _a[7] < _b[7] ? ~UINT16_C(0) : 0;

  __m128i a = load_m128i(_a);
  __m128i b = load_m128i(_b);
  __m128i c = _mm_cmplt_epi16(a, b);

  return VALIDATE_UINT16_M128(c, d);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_cmplt_epi32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int32_t *_a = impl.test_cases_int_pointer1;
  const int32_t *_b = impl.test_cases_int_pointer2;
  __m128i a = load_m128i(_a);
  __m128i b = load_m128i(_b);

  int32_t result[4];
  result[0] = _a[0] < _b[0] ? -1 : 0;
  result[1] = _a[1] < _b[1] ? -1 : 0;
  result[2] = _a[2] < _b[2] ? -1 : 0;
  result[3] = _a[3] < _b[3] ? -1 : 0;

  __m128i iret = _mm_cmplt_epi32(a, b);
  return VALIDATE_INT32_M128(iret, result);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_cmplt_epi8(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int8_t *_a = (const int8_t *)impl.test_cases_int_pointer1;
  const int8_t *_b = (const int8_t *)impl.test_cases_int_pointer2;
  int8_t d[16];
  d[0] = (_a[0] < _b[0]) ? ~UINT8_C(0) : 0x00;
  d[1] = (_a[1] < _b[1]) ? ~UINT8_C(0) : 0x00;
  d[2] = (_a[2] < _b[2]) ? ~UINT8_C(0) : 0x00;
  d[3] = (_a[3] < _b[3]) ? ~UINT8_C(0) : 0x00;
  d[4] = (_a[4] < _b[4]) ? ~UINT8_C(0) : 0x00;
  d[5] = (_a[5] < _b[5]) ? ~UINT8_C(0) : 0x00;
  d[6] = (_a[6] < _b[6]) ? ~UINT8_C(0) : 0x00;
  d[7] = (_a[7] < _b[7]) ? ~UINT8_C(0) : 0x00;
  d[8] = (_a[8] < _b[8]) ? ~UINT8_C(0) : 0x00;
  d[9] = (_a[9] < _b[9]) ? ~UINT8_C(0) : 0x00;
  d[10] = (_a[10] < _b[10]) ? ~UINT8_C(0) : 0x00;
  d[11] = (_a[11] < _b[11]) ? ~UINT8_C(0) : 0x00;
  d[12] = (_a[12] < _b[12]) ? ~UINT8_C(0) : 0x00;
  d[13] = (_a[13] < _b[13]) ? ~UINT8_C(0) : 0x00;
  d[14] = (_a[14] < _b[14]) ? ~UINT8_C(0) : 0x00;
  d[15] = (_a[15] < _b[15]) ? ~UINT8_C(0) : 0x00;

  __m128i a = load_m128i(_a);
  __m128i b = load_m128i(_b);
  __m128i c = _mm_cmplt_epi8(a, b);
  return VALIDATE_INT8_M128(c, d);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_cmplt_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const double *_a = (const double *)impl.test_cases_float_pointer1;
  const double *_b = (const double *)impl.test_cases_float_pointer2;

  int64_t f0 = (_a[0] < _b[0]) ? ~UINT64_C(0) : UINT64_C(0);
  int64_t f1 = (_a[1] < _b[1]) ? ~UINT64_C(0) : UINT64_C(0);

  __m128d a = load_m128d(_a);
  __m128d b = load_m128d(_b);
  __m128d c = _mm_cmplt_pd(a, b);

  return validate_double(c, *(double *)&f0, *(double *)&f1);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_cmplt_sd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  double *_a = (double *)impl.test_cases_float_pointer1;
  double *_b = (double *)impl.test_cases_float_pointer2;
  uint64_t d0 = (_a[0] < _b[0]) ? ~UINT64_C(0) : 0;
  uint64_t d1 = ((uint64_t *)_a)[1];

  __m128d a = load_m128d(_a);
  __m128d b = load_m128d(_b);
  __m128d c = _mm_cmplt_sd(a, b);

  return validate_double(c, *(double *)&d0, *(double *)&d1);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_cmpneq_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const double *_a = (const double *)impl.test_cases_float_pointer1;
  const double *_b = (const double *)impl.test_cases_float_pointer2;

  int64_t f0 = (_a[0] != _b[0]) ? ~UINT64_C(0) : UINT64_C(0);
  int64_t f1 = (_a[1] != _b[1]) ? ~UINT64_C(0) : UINT64_C(0);

  __m128d a = load_m128d(_a);
  __m128d b = load_m128d(_b);
  __m128d c = _mm_cmpneq_pd(a, b);

  return validate_double(c, *(double *)&f0, *(double *)&f1);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_cmpneq_sd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  double *_a = (double *)impl.test_cases_float_pointer1;
  double *_b = (double *)impl.test_cases_float_pointer2;

  int64_t f0 = (_a[0] != _b[0]) ? ~UINT64_C(0) : UINT64_C(0);
  int64_t f1 = ((int64_t *)_a)[1];

  __m128d a = load_m128d(_a);
  __m128d b = load_m128d(_b);
  __m128d c = _mm_cmpneq_sd(a, b);

  return validate_double(c, *(double *)&f0, *(double *)&f1);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_cmpnge_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  // const double *_a = (const double *)impl.test_cases_float_pointer1;
  // const double *_b = (const double *)impl.test_cases_float_pointer2;
  // uint64_t d0 = !(_a[0] >= _b[0]) ? ~UINT64_C(0) : 0;
  // uint64_t d1 = !(_a[1] >= _b[1]) ? ~UINT64_C(0) : 0;

  // __m128d a = load_m128d(_a);
  // __m128d b = load_m128d(_b);
  // __m128d c = _mm_cmpnge_pd(a, b);

  // return validate_double(c, *(double *)&d0, *(double *)&d1);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cmpnge_sd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  // double *_a = (double *)impl.test_cases_float_pointer1;
  // double *_b = (double *)impl.test_cases_float_pointer2;
  // uint64_t d0 = !(_a[0] >= _b[0]) ? ~UINT64_C(0) : 0;
  // uint64_t d1 = ((uint64_t *)_a)[1];

  // __m128d a = load_m128d(_a);
  // __m128d b = load_m128d(_b);
  // __m128d c = _mm_cmpnge_sd(a, b);

  // return validate_double(c, *(double *)&d0, *(double *)&d1);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cmpngt_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  // const double *_a = (const double *)impl.test_cases_float_pointer1;
  // const double *_b = (const double *)impl.test_cases_float_pointer2;
  // uint64_t d0 = !(_a[0] > _b[0]) ? ~UINT64_C(0) : 0;
  // uint64_t d1 = !(_a[1] > _b[1]) ? ~UINT64_C(0) : 0;

  // __m128d a = load_m128d(_a);
  // __m128d b = load_m128d(_b);
  // __m128d c = _mm_cmpngt_pd(a, b);

  // return validate_double(c, *(double *)&d0, *(double *)&d1);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cmpngt_sd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  // double *_a = (double *)impl.test_cases_float_pointer1;
  // double *_b = (double *)impl.test_cases_float_pointer2;
  // uint64_t d0 = !(_a[0] > _b[0]) ? ~UINT64_C(0) : 0;
  // uint64_t d1 = ((uint64_t *)_a)[1];

  // __m128d a = load_m128d(_a);
  // __m128d b = load_m128d(_b);
  // __m128d c = _mm_cmpngt_sd(a, b);

  // return validate_double(c, *(double *)&d0, *(double *)&d1);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cmpnle_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  // const double *_a = (const double *)impl.test_cases_float_pointer1;
  // const double *_b = (const double *)impl.test_cases_float_pointer2;
  // uint64_t d0 = !(_a[0] <= _b[0]) ? ~UINT64_C(0) : 0;
  // uint64_t d1 = !(_a[1] <= _b[1]) ? ~UINT64_C(0) : 0;

  // __m128d a = load_m128d(_a);
  // __m128d b = load_m128d(_b);
  // __m128d c = _mm_cmpnle_pd(a, b);

  // return validate_double(c, *(double *)&d0, *(double *)&d1);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cmpnle_sd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  // double *_a = (double *)impl.test_cases_float_pointer1;
  // double *_b = (double *)impl.test_cases_float_pointer2;
  // uint64_t d0 = !(_a[0] <= _b[0]) ? ~UINT64_C(0) : 0;
  // uint64_t d1 = ((uint64_t *)_a)[1];

  // __m128d a = load_m128d(_a);
  // __m128d b = load_m128d(_b);
  // __m128d c = _mm_cmpnle_sd(a, b);

  // return validate_double(c, *(double *)&d0, *(double *)&d1);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cmpnlt_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  // const double *_a = (const double *)impl.test_cases_float_pointer1;
  // const double *_b = (const double *)impl.test_cases_float_pointer2;
  // uint64_t d0 = !(_a[0] < _b[0]) ? ~UINT64_C(0) : 0;
  // uint64_t d1 = !(_a[1] < _b[1]) ? ~UINT64_C(0) : 0;

  // __m128d a = load_m128d(_a);
  // __m128d b = load_m128d(_b);
  // __m128d c = _mm_cmpnlt_pd(a, b);

  // return validate_double(c, *(double *)&d0, *(double *)&d1);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cmpnlt_sd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  // double *_a = (double *)impl.test_cases_float_pointer1;
  // double *_b = (double *)impl.test_cases_float_pointer2;
  // uint64_t d0 = !(_a[0] < _b[0]) ? ~UINT64_C(0) : 0;
  // uint64_t d1 = ((uint64_t *)_a)[1];

  // __m128d a = load_m128d(_a);
  // __m128d b = load_m128d(_b);
  // __m128d c = _mm_cmpnlt_sd(a, b);

  // return validate_double(c, *(double *)&d0, *(double *)&d1);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cmpord_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const double *_a = (const double *)impl.test_cases_float_pointer1;
  const double *_b = (const double *)impl.test_cases_float_pointer2;
  __m128d a = _mm_load_pd(_a);
  __m128d b = _mm_load_pd(_b);

  double _c[2];
  for (uint32_t i = 0; i < 2; i++) {
    _c[i] = cmp_noNaN(_a[i], _b[i]);
  }
  __m128d c = _mm_cmpord_pd(a, b);

  return validate_double(c, _c[0], _c[1]);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_cmpord_sd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const double *_a = (const double *)impl.test_cases_float_pointer1;
  const double *_b = (const double *)impl.test_cases_float_pointer2;
  __m128d a = _mm_load_pd(_a);
  __m128d b = _mm_load_pd(_b);

  double c0 = cmp_noNaN(_a[0], _b[0]);
  double c1 = _a[1];

  __m128d ret = _mm_cmpord_sd(a, b);
  return validate_double(ret, c0, c1);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_cmpunord_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const double *_a = (const double *)impl.test_cases_float_pointer1;
  const double *_b = (const double *)impl.test_cases_float_pointer2;
  __m128d a = _mm_load_pd(_a);
  __m128d b = _mm_load_pd(_b);

  double _c[2];
  _c[0] = cmp_hasNaN(_a[0], _b[0]);
  _c[1] = cmp_hasNaN(_a[1], _b[1]);
  __m128d c = _mm_cmpunord_pd(a, b);

  return validate_double(c, _c[0], _c[1]);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_cmpunord_sd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  double *_a = (double *)impl.test_cases_float_pointer1;
  double *_b = (double *)impl.test_cases_float_pointer2;
  __m128d a = _mm_load_pd(_a);
  __m128d b = _mm_load_pd(_b);

  double result[2];
  result[0] = cmp_hasNaN(_a[0], _b[0]);
  result[1] = _a[1];

  __m128d ret = _mm_cmpunord_sd(a, b);
  return validate_double(ret, result[0], result[1]);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_comieq_sd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
// FIXME:
// The GCC does not implement _mm_comieq_sd correctly.
// See https://gcc.gnu.org/bugzilla/show_bug.cgi?id=98612 for more
// information.
#if defined(__GNUC__) && !defined(__clang__)
  return TEST_UNIMPL;
#else
  const double *_a = (const double *)impl.test_cases_float_pointer1;
  const double *_b = (const double *)impl.test_cases_float_pointer2;
  int32_t _c = (_a[0] == _b[0]) ? 1 : 0;

  __m128d a = load_m128d(_a);
  __m128d b = load_m128d(_b);
  int32_t c = _mm_comieq_sd(a, b);

  ASSERT_RETURN(c == _c);
  return TEST_SUCCESS;
#endif
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_comige_sd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const double *_a = (const double *)impl.test_cases_float_pointer1;
  const double *_b = (const double *)impl.test_cases_float_pointer2;
  int32_t _c = (_a[0] >= _b[0]) ? 1 : 0;

  __m128d a = load_m128d(_a);
  __m128d b = load_m128d(_b);
  int32_t c = _mm_comige_sd(a, b);

  ASSERT_RETURN(c == _c);
  return TEST_SUCCESS;
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_comigt_sd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const double *_a = (const double *)impl.test_cases_float_pointer1;
  const double *_b = (const double *)impl.test_cases_float_pointer2;
  int32_t _c = (_a[0] > _b[0]) ? 1 : 0;

  __m128d a = load_m128d(_a);
  __m128d b = load_m128d(_b);
  int32_t c = _mm_comigt_sd(a, b);

  ASSERT_RETURN(c == _c);
  return TEST_SUCCESS;
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_comile_sd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
// FIXME:
// The GCC does not implement _mm_comile_sd correctly.
// See https://gcc.gnu.org/bugzilla/show_bug.cgi?id=98612 for more
// information.
#if defined(__GNUC__) && !defined(__clang__)
  return TEST_UNIMPL;
#else
  const double *_a = (const double *)impl.test_cases_float_pointer1;
  const double *_b = (const double *)impl.test_cases_float_pointer2;
  int32_t _c = (_a[0] <= _b[0]) ? 1 : 0;

  __m128d a = load_m128d(_a);
  __m128d b = load_m128d(_b);
  int32_t c = _mm_comile_sd(a, b);

  ASSERT_RETURN(c == _c);
  return TEST_SUCCESS;
#endif
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_comilt_sd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
// FIXME:
// The GCC does not implement _mm_comilt_sd correctly.
// See https://gcc.gnu.org/bugzilla/show_bug.cgi?id=98612 for more
// information.
#if defined(__GNUC__) && !defined(__clang__)
  return TEST_UNIMPL;
#else
  const double *_a = (const double *)impl.test_cases_float_pointer1;
  const double *_b = (const double *)impl.test_cases_float_pointer2;
  int32_t _c = (_a[0] < _b[0]) ? 1 : 0;

  __m128d a = load_m128d(_a);
  __m128d b = load_m128d(_b);
  int32_t c = _mm_comilt_sd(a, b);

  ASSERT_RETURN(c == _c);
  return TEST_SUCCESS;
#endif
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_comineq_sd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
// FIXME:
// The GCC does not implement _mm_comineq_sd correctly.
// See https://gcc.gnu.org/bugzilla/show_bug.cgi?id=98612 for more
// information.
#if defined(__GNUC__) && !defined(__clang__)
  return TEST_UNIMPL;
#else
  const double *_a = (const double *)impl.test_cases_float_pointer1;
  const double *_b = (const double *)impl.test_cases_float_pointer2;
  int32_t _c = (_a[0] != _b[0]) ? 1 : 0;

  __m128d a = load_m128d(_a);
  __m128d b = load_m128d(_b);
  int32_t c = _mm_comineq_sd(a, b);

  ASSERT_RETURN(c == _c);
  return TEST_SUCCESS;
#endif
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_cvtepi32_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int32_t *_a = impl.test_cases_int_pointer1;
  //   __m128i a = load_m128i(_a);
  //   double trun[2] = {(double)_a[0], (double)_a[1]};
  //
  //   __m128d ret = _mm_cvtepi32_pd(a);
  //   return validate_double(ret, trun[0], trun[1]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtepi32_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int32_t *_a = impl.test_cases_int_pointer1;
  //   __m128i a = load_m128i(_a);
  //   float trun[4];
  //   for (uint32_t i = 0; i < 4; i++) {
  //     trun[i] = (float)_a[i];
  //   }
  //
  //   __m128 ret = _mm_cvtepi32_ps(a);
  //   return validate_float(ret, trun[0], trun[1], trun[2], trun[3]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtpd_epi32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const double *_a = (const double *)impl.test_cases_float_pointer1;
  //   int32_t d[2];
  //
  //   switch (iter & 0x3) {
  //   case 0:
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_NEAREST);
  //     d[0] = (int32_t)(bankers_rounding(_a[0]));
  //     d[1] = (int32_t)(bankers_rounding(_a[1]));
  //     break;
  //   case 1:
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_DOWN);
  //     d[0] = (int32_t)(floor(_a[0]));
  //     d[1] = (int32_t)(floor(_a[1]));
  //     break;
  //   case 2:
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
  //     d[0] = (int32_t)(ceil(_a[0]));
  //     d[1] = (int32_t)(ceil(_a[1]));
  //     break;
  //   case 3:
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_TOWARD_ZERO);
  //     d[0] = (int32_t)(_a[0]);
  //     d[1] = (int32_t)(_a[1]);
  //     break;
  //   }
  //
  //   __m128d a = load_m128d(_a);
  //   __m128i ret = _mm_cvtpd_epi32(a);
  //
  //   return validate_int32(ret, d[0], d[1], 0, 0);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtpd_pi32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const double *_a = (const double *)impl.test_cases_float_pointer1;
  //   int32_t d[2];
  //
  //   switch (iter & 0x3) {
  //   case 0:
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_NEAREST);
  //     d[0] = (int32_t)(bankers_rounding(_a[0]));
  //     d[1] = (int32_t)(bankers_rounding(_a[1]));
  //     break;
  //   case 1:
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_DOWN);
  //     d[0] = (int32_t)(floor(_a[0]));
  //     d[1] = (int32_t)(floor(_a[1]));
  //     break;
  //   case 2:
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
  //     d[0] = (int32_t)(ceil(_a[0]));
  //     d[1] = (int32_t)(ceil(_a[1]));
  //     break;
  //   case 3:
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_TOWARD_ZERO);
  //     d[0] = (int32_t)(_a[0]);
  //     d[1] = (int32_t)(_a[1]);
  //     break;
  //   }
  //
  //   __m128d a = load_m128d(_a);
  //   __m64 ret = _mm_cvtpd_pi32(a);
  //
  //   return VALIDATE_INT32_M64(ret, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtpd_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const double *_a = (const double *)impl.test_cases_float_pointer1;
  //   float f0 = (float)_a[0];
  //   float f1 = (float)_a[1];
  //   const __m128d a = load_m128d(_a);
  //
  //   __m128 r = _mm_cvtpd_ps(a);
  //
  //   return validate_float(r, f0, f1, 0, 0);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtpi32_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int32_t *_a = impl.test_cases_int_pointer1;
  //   __m64 a = load_m64(_a);
  //
  //   double trun[2] = {(double)_a[0], (double)_a[1]};
  //
  //   __m128d ret = _mm_cvtpi32_pd(a);
  //
  //   return validate_double(ret, trun[0], trun[1]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtps_epi32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   __m128 a = load_m128(_a);
  //   int32_t d[4];
  //   switch (iter & 0x3) {
  //   case 0:
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_NEAREST);
  //     for (uint32_t i = 0; i < 4; i++) {
  //       d[i] = (int32_t)(bankers_rounding(_a[i]));
  //     }
  //     break;
  //   case 1:
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_DOWN);
  //     for (uint32_t i = 0; i < 4; i++) {
  //       d[i] = (int32_t)(floorf(_a[i]));
  //     }
  //     break;
  //   case 2:
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
  //     for (uint32_t i = 0; i < 4; i++) {
  //       d[i] = (int32_t)(ceilf(_a[i]));
  //     }
  //     break;
  //   case 3:
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_TOWARD_ZERO);
  //     for (uint32_t i = 0; i < 4; i++) {
  //       d[i] = (int32_t)(_a[i]);
  //     }
  //     break;
  //   }
  //
  //   __m128i ret = _mm_cvtps_epi32(a);
  //   return VALIDATE_INT32_M128(ret, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtps_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   double d0 = (double)_a[0];
  //   double d1 = (double)_a[1];
  //   const __m128 a = load_m128(_a);
  //
  //   __m128d r = _mm_cvtps_pd(a);
  //
  //   return validate_double(r, d0, d1);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtsd_f64(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const double *_a = (const double *)impl.test_cases_float_pointer1;
  //
  //   double d = _a[0];
  //
  //   const __m128d *a = (const __m128d *)_a;
  //   double r = _mm_cvtsd_f64(*a);
  //
  //   return r == d ? TEST_SUCCESS : TEST_FAIL;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtsd_si32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const double *_a = (const double *)impl.test_cases_float_pointer1;
  //   int32_t d;
  //
  //   switch (iter & 0x3) {
  //   case 0:
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_NEAREST);
  //     d = (int32_t)(bankers_rounding(_a[0]));
  //     break;
  //   case 1:
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_DOWN);
  //     d = (int32_t)(floor(_a[0]));
  //     break;
  //   case 2:
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
  //     d = (int32_t)(ceil(_a[0]));
  //     break;
  //   case 3:
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_TOWARD_ZERO);
  //     d = (int32_t)(_a[0]);
  //     break;
  //   }
  //
  //   __m128d a = load_m128d(_a);
  //   int32_t ret = _mm_cvtsd_si32(a);
  //
  //   return ret == d ? TEST_SUCCESS : TEST_FAIL;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtsd_si64(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const double *_a = (const double *)impl.test_cases_float_pointer1;
  //   int64_t d;
  //
  //   switch (iter & 0x3) {
  //   case 0:
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_NEAREST);
  //     d = (int64_t)(bankers_rounding(_a[0]));
  //     break;
  //   case 1:
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_DOWN);
  //     d = (int64_t)(floor(_a[0]));
  //     break;
  //   case 2:
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
  //     d = (int64_t)(ceil(_a[0]));
  //     break;
  //   case 3:
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_TOWARD_ZERO);
  //     d = (int64_t)(_a[0]);
  //     break;
  //   }
  //
  //   __m128d a = load_m128d(_a);
  //   int64_t ret = _mm_cvtsd_si64(a);
  //
  //   return ret == d ? TEST_SUCCESS : TEST_FAIL;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtsd_si64x(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   return test_mm_cvtsd_si64(impl, iter);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtsd_ss(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   const double *_b = (const double *)impl.test_cases_float_pointer2;
  //
  //   float f0 = _b[0];
  //   float f1 = _a[1];
  //   float f2 = _a[2];
  //   float f3 = _a[3];
  //
  //   __m128 a = load_m128(_a);
  //   __m128d b = load_m128d(_b);
  //   __m128 c = _mm_cvtsd_ss(a, b);
  //
  //   return validate_float(c, f0, f1, f2, f3);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtsi128_si32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int32_t *_a = (const int32_t *)impl.test_cases_int_pointer1;
  //
  //   int32_t d = _a[0];
  //
  //   __m128i a = load_m128i(_a);
  //   int c = _mm_cvtsi128_si32(a);
  //
  //   return d == c ? TEST_SUCCESS : TEST_FAIL;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtsi128_si64(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int64_t *_a = (const int64_t *)impl.test_cases_int_pointer1;
  //
  //   int64_t d = _a[0];
  //
  //   __m128i a = load_m128i(_a);
  //   int64_t c = _mm_cvtsi128_si64(a);
  //
  //   return d == c ? TEST_SUCCESS : TEST_FAIL;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtsi128_si64x(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   return test_mm_cvtsi128_si64(impl, iter);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtsi32_sd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const double *_a = (const double *)impl.test_cases_float_pointer1;
  //   const int32_t b = (const int32_t)impl.test_cases_ints[iter];
  //
  //   __m128d a = load_m128d(_a);
  //   __m128d c = _mm_cvtsi32_sd(a, b);
  //
  //   return validate_double(c, b, _a[1]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtsi32_si128(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int32_t *_a = (const int32_t *)impl.test_cases_int_pointer1;
  //
  //   int32_t d = _a[0];
  //
  //   __m128i c = _mm_cvtsi32_si128(*_a);
  //
  //   return validate_int32(c, d, 0, 0, 0);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtsi64_sd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const double *_a = (const double *)impl.test_cases_float_pointer1;
  //   const int64_t b = (const int64_t)impl.test_cases_ints[iter];
  //
  //   __m128d a = load_m128d(_a);
  //   __m128d c = _mm_cvtsi64_sd(a, b);
  //
  //   return validate_double(c, b, _a[1]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtsi64_si128(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int64_t *_a = (const int64_t *)impl.test_cases_int_pointer1;
  //
  //   int64_t d = _a[0];
  //
  //   __m128i c = _mm_cvtsi64_si128(*_a);
  //
  //   return validate_int64(c, d, 0);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtsi64x_sd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   return test_mm_cvtsi64_sd(impl, iter);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtsi64x_si128(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   return test_mm_cvtsi64_si128(impl, iter);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtss_sd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const double *_a = (const double *)impl.test_cases_float_pointer1;
  //   const float *_b = impl.test_cases_float_pointer2;
  //   double d0 = double(_b[0]);
  //   double d1 = _a[1];
  //
  //   __m128d a = load_m128d(_a);
  //   __m128 b = load_m128(_b);
  //   __m128d c = _mm_cvtss_sd(a, b);
  //   return validate_double(c, d0, d1);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvttpd_epi32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const double *_a = (const double *)impl.test_cases_float_pointer1;
  //
  //   __m128d a = load_m128d(_a);
  //   int32_t d0 = (int32_t)(_a[0]);
  //   int32_t d1 = (int32_t)(_a[1]);
  //
  //   __m128i ret = _mm_cvttpd_epi32(a);
  //   return validate_int32(ret, d0, d1, 0, 0);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvttpd_pi32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const double *_a = (const double *)impl.test_cases_float_pointer1;
  //
  //   __m128d a = load_m128d(_a);
  //   int32_t d0 = (int32_t)(_a[0]);
  //   int32_t d1 = (int32_t)(_a[1]);
  //
  //   __m64 ret = _mm_cvttpd_pi32(a);
  //   return validate_int32(ret, d0, d1);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvttps_epi32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   __m128 a = load_m128(_a);
  //   int32_t trun[4];
  //   for (uint32_t i = 0; i < 4; i++) {
  //     trun[i] = (int32_t)_a[i];
  //   }
  //
  //   __m128i ret = _mm_cvttps_epi32(a);
  //   return VALIDATE_INT32_M128(ret, trun);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvttsd_si32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const double *_a = (const double *)impl.test_cases_float_pointer1;
  //
  //   __m128d a = _mm_load_sd(_a);
  //   int32_t ret = _mm_cvttsd_si32(a);
  //
  //   return ret == (int32_t)_a[0] ? TEST_SUCCESS : TEST_FAIL;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvttsd_si64(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const double *_a = (const double *)impl.test_cases_float_pointer1;
  //
  //   __m128d a = _mm_load_sd(_a);
  //   int64_t ret = _mm_cvttsd_si64(a);
  //
  //   return ret == (int64_t)_a[0] ? TEST_SUCCESS : TEST_FAIL;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvttsd_si64x(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  // #if defined(__clang__)
  // The intrinsic _mm_cvttsd_si64x() does not exist in Clang
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
  // #else
  //   const double *_a = (const double *)impl.test_cases_float_pointer1;
  //
  //   __m128d a = _mm_load_sd(_a);
  //   int64_t ret = _mm_cvttsd_si64x(a);
  //
  //   return ret == (int64_t)_a[0] ? TEST_SUCCESS : TEST_FAIL;
  // #endif
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_div_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const double *_a = (const double *)impl.test_cases_float_pointer1;
  //   const double *_b = (const double *)impl.test_cases_float_pointer2;
  //   double d0 = 0.0, d1 = 0.0;
  //
  //   if (_b[0] != 0.0)
  //     d0 = _a[0] / _b[0];
  //   if (_b[1] != 0.0)
  //     d1 = _a[1] / _b[1];
  //
  //   __m128d a = load_m128d(_a);
  //   __m128d b = load_m128d(_b);
  //   __m128d c = _mm_div_pd(a, b);
  //   return validate_double(c, d0, d1);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_div_sd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const double *_a = (const double *)impl.test_cases_float_pointer1;
  //   const double *_b = (const double *)impl.test_cases_float_pointer2;
  //
  //   double d0 = _a[0] / _b[0];
  //   double d1 = _a[1];
  //
  //   __m128d a = load_m128d(_a);
  //   __m128d b = load_m128d(_b);
  //
  //   __m128d c = _mm_div_sd(a, b);
  //
  //   return validate_double(c, d0, d1);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_extract_epi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   uint16_t *_a = (uint16_t *)impl.test_cases_int_pointer1;
  //   const int idx = iter & 0x7;
  //   __m128i a = load_m128i(_a);
  //   int c;
  //   switch (idx) {
  //   case 0:
  //     c = _mm_extract_epi16(a, 0);
  //     break;
  //   case 1:
  //     c = _mm_extract_epi16(a, 1);
  //     break;
  //   case 2:
  //     c = _mm_extract_epi16(a, 2);
  //     break;
  //   case 3:
  //     c = _mm_extract_epi16(a, 3);
  //     break;
  //   case 4:
  //     c = _mm_extract_epi16(a, 4);
  //     break;
  //   case 5:
  //     c = _mm_extract_epi16(a, 5);
  //     break;
  //   case 6:
  //     c = _mm_extract_epi16(a, 6);
  //     break;
  //   case 7:
  //     c = _mm_extract_epi16(a, 7);
  //     break;
  //   }
  //
  //   ASSERT_RETURN(c == *(_a + idx));
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_insert_epi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  //   const int16_t insert = (int16_t)*impl.test_cases_int_pointer2;
  //
  // #define TEST_IMPL(IDX)
  //   int16_t d##IDX[8];
  //   for (int i = 0; i < 8; i++) {
  //     d##IDX[i] = _a[i];
  //   }
  //   d##IDX[IDX] = insert;
  //
  //   __m128i a##IDX = load_m128i(_a);
  //   __m128i b##IDX = _mm_insert_epi16(a##IDX, insert, IDX);
  //   CHECK_RESULT(VALIDATE_INT16_M128(b##IDX, d##IDX))
  //
  //   IMM_8_ITER
  // #undef TEST_IMPL
  //
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_lfence(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   /* FIXME: Assume that memory barriers always function as intended. */
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_load_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const double *p = (const double *)impl.test_cases_float_pointer1;
  __m128d a = _mm_load_pd(p);
  return validate_double(a, p[0], p[1]);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_load_pd1(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const double *p = (const double *)impl.test_cases_float_pointer1;
  __m128d a = _mm_load_pd1(p);
  return validate_double(a, p[0], p[0]);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_load_sd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const double *p = (const double *)impl.test_cases_float_pointer1;
  __m128d a = _mm_load_sd(p);
  return validate_double(a, p[0], 0);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_load_si128(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int32_t *addr = impl.test_cases_int_pointer1;

  __m128i ret = _mm_load_si128((const __m128i *)addr);

  return VALIDATE_INT32_M128(ret, addr);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_load1_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const double *addr = (const double *)impl.test_cases_float_pointer1;

  __m128d ret = _mm_load1_pd(addr);

  return validate_double(ret, addr[0], addr[0]);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_loadh_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const double *_a = (const double *)impl.test_cases_float_pointer1;
  //   const double *addr = (const double *)impl.test_cases_float_pointer2;
  //
  //   __m128d a = load_m128d(_a);
  //   __m128d ret = _mm_loadh_pd(a, addr);
  //
  //   return validate_double(ret, _a[0], addr[0]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_loadl_epi64(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int64_t *addr = (const int64_t *)impl.test_cases_int_pointer1;
  //
  //   __m128i ret = _mm_loadl_epi64((const __m128i *)addr);
  //
  //   return validate_int64(ret, addr[0], 0);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_loadl_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const double *_a = (const double *)impl.test_cases_float_pointer1;
  //   const double *addr = (const double *)impl.test_cases_float_pointer2;
  //
  //   __m128d a = load_m128d(_a);
  //   __m128d ret = _mm_loadl_pd(a, addr);
  //
  //   return validate_double(ret, addr[0], _a[1]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_loadr_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const double *addr = (const double *)impl.test_cases_float_pointer1;
  //
  //   __m128d ret = _mm_loadr_pd(addr);
  //
  //   return validate_double(ret, addr[1], addr[0]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_loadu_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const double *p = (const double *)impl.test_cases_float_pointer1;
  //   __m128d a = _mm_loadu_pd(p);
  //   return validate_double(a, p[0], p[1]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_loadu_si128(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int32_t *_a = (const int32_t *)impl.test_cases_int_pointer1;
  //   __m128i c = _mm_loadu_si128((const __m128i *)_a);
  //   return VALIDATE_INT32_M128(c, _a);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_loadu_si32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  // The GCC version before 11 does not implement intrinsic function
  // _mm_loadu_si32. Check https://gcc.gnu.org/bugzilla/show_bug.cgi?id=95483
  // for more information.
  // #if (defined(__GNUC__) && !defined(__clang__)) && (__GNUC__ <= 10)
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
  // #else
  //   const int32_t *addr = (const int32_t *)impl.test_cases_int_pointer1;
  //
  //   __m128i ret = _mm_loadu_si32((const void *)addr);
  //
  //   return validate_int32(ret, addr[0], 0, 0, 0);
  // #endif
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_madd_epi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  //   const int16_t *_b = (const int16_t *)impl.test_cases_int_pointer2;
  //   int32_t d0 = (int32_t)_a[0] * _b[0];
  //   int32_t d1 = (int32_t)_a[1] * _b[1];
  //   int32_t d2 = (int32_t)_a[2] * _b[2];
  //   int32_t d3 = (int32_t)_a[3] * _b[3];
  //   int32_t d4 = (int32_t)_a[4] * _b[4];
  //   int32_t d5 = (int32_t)_a[5] * _b[5];
  //   int32_t d6 = (int32_t)_a[6] * _b[6];
  //   int32_t d7 = (int32_t)_a[7] * _b[7];
  //
  //   int32_t e[4];
  //   e[0] = d0 + d1;
  //   e[1] = d2 + d3;
  //   e[2] = d4 + d5;
  //   e[3] = d6 + d7;
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //   __m128i c = _mm_madd_epi16(a, b);
  //   return VALIDATE_INT32_M128(c, e);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_maskmoveu_si128(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const uint8_t *_a = (const uint8_t *)impl.test_cases_int_pointer1;
  //   const uint8_t *_mask = (const uint8_t *)impl.test_cases_int_pointer2;
  //   char mem_addr[16];
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i mask = load_m128i(_mask);
  //   _mm_maskmoveu_si128(a, mask, mem_addr);
  //
  //   for (int i = 0; i < 16; i++) {
  //     if (_mask[i] >> 7) {
  //       ASSERT_RETURN(_a[i] == (uint8_t)mem_addr[i]);
  //     }
  //   }
  //
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_max_epi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  //   const int16_t *_b = (const int16_t *)impl.test_cases_int_pointer2;
  //   int16_t d[8];
  //   d[0] = _a[0] > _b[0] ? _a[0] : _b[0];
  //   d[1] = _a[1] > _b[1] ? _a[1] : _b[1];
  //   d[2] = _a[2] > _b[2] ? _a[2] : _b[2];
  //   d[3] = _a[3] > _b[3] ? _a[3] : _b[3];
  //   d[4] = _a[4] > _b[4] ? _a[4] : _b[4];
  //   d[5] = _a[5] > _b[5] ? _a[5] : _b[5];
  //   d[6] = _a[6] > _b[6] ? _a[6] : _b[6];
  //   d[7] = _a[7] > _b[7] ? _a[7] : _b[7];
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //
  //   __m128i c = _mm_max_epi16(a, b);
  //   return VALIDATE_INT16_M128(c, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_max_epu8(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int8_t *_a = (const int8_t *)impl.test_cases_int_pointer1;
  //   const int8_t *_b = (const int8_t *)impl.test_cases_int_pointer2;
  //   uint8_t d[16];
  //   d[0] =
  //       ((uint8_t)_a[0] > (uint8_t)_b[0]) ? ((uint8_t)_a[0]) :
  //       ((uint8_t)_b[0]);
  //   d[1] =
  //       ((uint8_t)_a[1] > (uint8_t)_b[1]) ? ((uint8_t)_a[1]) :
  //       ((uint8_t)_b[1]);
  //   d[2] =
  //       ((uint8_t)_a[2] > (uint8_t)_b[2]) ? ((uint8_t)_a[2]) :
  //       ((uint8_t)_b[2]);
  //   d[3] =
  //       ((uint8_t)_a[3] > (uint8_t)_b[3]) ? ((uint8_t)_a[3]) :
  //       ((uint8_t)_b[3]);
  //   d[4] =
  //       ((uint8_t)_a[4] > (uint8_t)_b[4]) ? ((uint8_t)_a[4]) :
  //       ((uint8_t)_b[4]);
  //   d[5] =
  //       ((uint8_t)_a[5] > (uint8_t)_b[5]) ? ((uint8_t)_a[5]) :
  //       ((uint8_t)_b[5]);
  //   d[6] =
  //       ((uint8_t)_a[6] > (uint8_t)_b[6]) ? ((uint8_t)_a[6]) :
  //       ((uint8_t)_b[6]);
  //   d[7] =
  //       ((uint8_t)_a[7] > (uint8_t)_b[7]) ? ((uint8_t)_a[7]) :
  //       ((uint8_t)_b[7]);
  //   d[8] =
  //       ((uint8_t)_a[8] > (uint8_t)_b[8]) ? ((uint8_t)_a[8]) :
  //       ((uint8_t)_b[8]);
  //   d[9] =
  //       ((uint8_t)_a[9] > (uint8_t)_b[9]) ? ((uint8_t)_a[9]) :
  //       ((uint8_t)_b[9]);
  //   d[10] = ((uint8_t)_a[10] > (uint8_t)_b[10]) ? ((uint8_t)_a[10])
  //                                               : ((uint8_t)_b[10]);
  //   d[11] = ((uint8_t)_a[11] > (uint8_t)_b[11]) ? ((uint8_t)_a[11])
  //                                               : ((uint8_t)_b[11]);
  //   d[12] = ((uint8_t)_a[12] > (uint8_t)_b[12]) ? ((uint8_t)_a[12])
  //                                               : ((uint8_t)_b[12]);
  //   d[13] = ((uint8_t)_a[13] > (uint8_t)_b[13]) ? ((uint8_t)_a[13])
  //                                               : ((uint8_t)_b[13]);
  //   d[14] = ((uint8_t)_a[14] > (uint8_t)_b[14]) ? ((uint8_t)_a[14])
  //                                               : ((uint8_t)_b[14]);
  //   d[15] = ((uint8_t)_a[15] > (uint8_t)_b[15]) ? ((uint8_t)_a[15])
  //                                               : ((uint8_t)_b[15]);
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //   __m128i c = _mm_max_epu8(a, b);
  //   return VALIDATE_INT8_M128(c, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_max_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const double *_a = (const double *)impl.test_cases_float_pointer1;
  //   const double *_b = (const double *)impl.test_cases_float_pointer2;
  //
  //   double f0 = _a[0] > _b[0] ? _a[0] : _b[0];
  //   double f1 = _a[1] > _b[1] ? _a[1] : _b[1];
  //
  //   __m128d a = load_m128d(_a);
  //   __m128d b = load_m128d(_b);
  //   __m128d c = _mm_max_pd(a, b);
  //
  //   return validate_double(c, f0, f1);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_max_sd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const double *_a = (const double *)impl.test_cases_float_pointer1;
  //   const double *_b = (const double *)impl.test_cases_float_pointer2;
  //   double d0 = _a[0] > _b[0] ? _a[0] : _b[0];
  //   double d1 = _a[1];
  //
  //   __m128d a = load_m128d(_a);
  //   __m128d b = load_m128d(_b);
  //   __m128d c = _mm_max_sd(a, b);
  //
  //   return validate_double(c, d0, d1);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_mfence(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   /* FIXME: Assume that memory barriers always function as intended. */
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_min_epi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  //   const int16_t *_b = (const int16_t *)impl.test_cases_int_pointer2;
  //   int16_t d[8];
  //   d[0] = _a[0] < _b[0] ? _a[0] : _b[0];
  //   d[1] = _a[1] < _b[1] ? _a[1] : _b[1];
  //   d[2] = _a[2] < _b[2] ? _a[2] : _b[2];
  //   d[3] = _a[3] < _b[3] ? _a[3] : _b[3];
  //   d[4] = _a[4] < _b[4] ? _a[4] : _b[4];
  //   d[5] = _a[5] < _b[5] ? _a[5] : _b[5];
  //   d[6] = _a[6] < _b[6] ? _a[6] : _b[6];
  //   d[7] = _a[7] < _b[7] ? _a[7] : _b[7];
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //   __m128i c = _mm_min_epi16(a, b);
  //   return VALIDATE_INT16_M128(c, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_min_epu8(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int8_t *_a = (const int8_t *)impl.test_cases_int_pointer1;
  //   const int8_t *_b = (const int8_t *)impl.test_cases_int_pointer2;
  //   uint8_t d[16];
  //   d[0] = ((uint8_t)_a[0] < (uint8_t)_b[0]) ? (uint8_t)_a[0] :
  //   (uint8_t)_b[0]; d[1] = ((uint8_t)_a[1] < (uint8_t)_b[1]) ? (uint8_t)_a[1]
  //   : (uint8_t)_b[1]; d[2] = ((uint8_t)_a[2] < (uint8_t)_b[2]) ?
  //   (uint8_t)_a[2] : (uint8_t)_b[2]; d[3] = ((uint8_t)_a[3] < (uint8_t)_b[3])
  //   ? (uint8_t)_a[3] : (uint8_t)_b[3]; d[4] = ((uint8_t)_a[4] <
  //   (uint8_t)_b[4]) ? (uint8_t)_a[4] : (uint8_t)_b[4]; d[5] = ((uint8_t)_a[5]
  //   < (uint8_t)_b[5]) ? (uint8_t)_a[5] : (uint8_t)_b[5]; d[6] =
  //   ((uint8_t)_a[6] < (uint8_t)_b[6]) ? (uint8_t)_a[6] : (uint8_t)_b[6]; d[7]
  //   = ((uint8_t)_a[7] < (uint8_t)_b[7]) ? (uint8_t)_a[7] : (uint8_t)_b[7];
  //   d[8] = ((uint8_t)_a[8] < (uint8_t)_b[8]) ? (uint8_t)_a[8] :
  //   (uint8_t)_b[8]; d[9] = ((uint8_t)_a[9] < (uint8_t)_b[9]) ? (uint8_t)_a[9]
  //   : (uint8_t)_b[9]; d[10] =
  //       ((uint8_t)_a[10] < (uint8_t)_b[10]) ? (uint8_t)_a[10] :
  //       (uint8_t)_b[10];
  //   d[11] =
  //       ((uint8_t)_a[11] < (uint8_t)_b[11]) ? (uint8_t)_a[11] :
  //       (uint8_t)_b[11];
  //   d[12] =
  //       ((uint8_t)_a[12] < (uint8_t)_b[12]) ? (uint8_t)_a[12] :
  //       (uint8_t)_b[12];
  //   d[13] =
  //       ((uint8_t)_a[13] < (uint8_t)_b[13]) ? (uint8_t)_a[13] :
  //       (uint8_t)_b[13];
  //   d[14] =
  //       ((uint8_t)_a[14] < (uint8_t)_b[14]) ? (uint8_t)_a[14] :
  //       (uint8_t)_b[14];
  //   d[15] =
  //       ((uint8_t)_a[15] < (uint8_t)_b[15]) ? (uint8_t)_a[15] :
  //       (uint8_t)_b[15];
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //   __m128i c = _mm_min_epu8(a, b);
  //   return VALIDATE_INT8_M128(c, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_min_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const double *_a = (const double *)impl.test_cases_float_pointer1;
  //   const double *_b = (const double *)impl.test_cases_float_pointer2;
  //   double f0 = _a[0] < _b[0] ? _a[0] : _b[0];
  //   double f1 = _a[1] < _b[1] ? _a[1] : _b[1];
  //
  //   __m128d a = load_m128d(_a);
  //   __m128d b = load_m128d(_b);
  //
  //   __m128d c = _mm_min_pd(a, b);
  //   return validate_double(c, f0, f1);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_min_sd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const double *_a = (const double *)impl.test_cases_float_pointer1;
  //   const double *_b = (const double *)impl.test_cases_float_pointer2;
  //   double d0 = _a[0] < _b[0] ? _a[0] : _b[0];
  //   double d1 = _a[1];
  //
  //   __m128d a = load_m128d(_a);
  //   __m128d b = load_m128d(_b);
  //   __m128d c = _mm_min_sd(a, b);
  //
  //   return validate_double(c, d0, d1);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_move_epi64(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int64_t *_a = (const int64_t *)impl.test_cases_int_pointer1;
  //
  //   int64_t d0 = _a[0];
  //   int64_t d1 = 0;
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i c = _mm_move_epi64(a);
  //
  //   return validate_int64(c, d0, d1);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_move_sd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const double *_a = (const double *)impl.test_cases_float_pointer1;
  //   const double *_b = (const double *)impl.test_cases_float_pointer2;
  //   __m128d a = load_m128d(_a);
  //   __m128d b = load_m128d(_b);
  //
  //   double result[2];
  //   result[0] = _b[0];
  //   result[1] = _a[1];
  //
  //   __m128d ret = _mm_move_sd(a, b);
  //   return validate_double(ret, result[0], result[1]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_movemask_epi8(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int32_t *_a = impl.test_cases_int_pointer1;
  //   __m128i a = load_m128i(_a);
  //
  //   const uint8_t *ip = (const uint8_t *)_a;
  //   int ret = 0;
  //   uint32_t mask = 1;
  //   for (uint32_t i = 0; i < 16; i++) {
  //     if (ip[i] & 0x80) {
  //       ret |= mask;
  //     }
  //     mask = mask << 1;
  //   }
  //   int test = _mm_movemask_epi8(a);
  //   ASSERT_RETURN(test == ret);
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_movemask_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const double *_a = (const double *)impl.test_cases_float_pointer1;
  //   unsigned int _c = 0;
  //   _c |= ((*(const uint64_t *)_a) >> 63) & 0x1;
  //   _c |= (((*(const uint64_t *)(_a + 1)) >> 62) & 0x2);
  //
  //   __m128d a = load_m128d(_a);
  //   int c = _mm_movemask_pd(a);
  //
  //   ASSERT_RETURN((unsigned int)c == _c);
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_movepi64_pi64(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int64_t *_a = (const int64_t *)impl.test_cases_int_pointer1;
  //
  //   int64_t d0 = _a[0];
  //
  //   __m128i a = load_m128i(_a);
  //   __m64 c = _mm_movepi64_pi64(a);
  //
  //   return validate_int64(c, d0);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_movpi64_epi64(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int64_t *_a = (const int64_t *)impl.test_cases_int_pointer1;
  //
  //   int64_t d0 = _a[0];
  //
  //   __m64 a = load_m64(_a);
  //   __m128i c = _mm_movpi64_epi64(a);
  //
  //   return validate_int64(c, d0, 0);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_mul_epu32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const uint32_t *_a = (const uint32_t *)impl.test_cases_int_pointer1;
  //   const uint32_t *_b = (const uint32_t *)impl.test_cases_int_pointer2;
  //   uint64_t dx = (uint64_t)(_a[0]) * (uint64_t)(_b[0]);
  //   uint64_t dy = (uint64_t)(_a[2]) * (uint64_t)(_b[2]);
  //
  //   __m128i a = _mm_loadu_si128((const __m128i *)_a);
  //   __m128i b = _mm_loadu_si128((const __m128i *)_b);
  //   __m128i r = _mm_mul_epu32(a, b);
  //   return validate_uint64(r, dx, dy);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_mul_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const double *_a = (const double *)impl.test_cases_float_pointer1;
  //   const double *_b = (const double *)impl.test_cases_float_pointer2;
  //   double d0 = _a[0] * _b[0];
  //   double d1 = _a[1] * _b[1];
  //
  //   __m128d a = _mm_load_pd(_a);
  //   __m128d b = _mm_load_pd(_b);
  //   __m128d c = _mm_mul_pd(a, b);
  //   return validate_double(c, d0, d1);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_mul_sd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const double *_a = (const double *)impl.test_cases_float_pointer1;
  //   const double *_b = (const double *)impl.test_cases_float_pointer2;
  //   double dx = _a[0] * _b[0];
  //   double dy = _a[1];
  //
  //   __m128d a = load_m128d(_a);
  //   __m128d b = load_m128d(_b);
  //   __m128d c = _mm_mul_sd(a, b);
  //   return validate_double(c, dx, dy);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_mul_su32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const uint32_t *_a = (const uint32_t *)impl.test_cases_int_pointer1;
  //   const uint32_t *_b = (const uint32_t *)impl.test_cases_int_pointer2;
  //
  //   uint64_t u = (uint64_t)(_a[0]) * (uint64_t)(_b[0]);
  //
  //   __m64 a = load_m64(_a);
  //   __m64 b = load_m64(_b);
  //   __m64 r = _mm_mul_su32(a, b);
  //
  //   return validate_uint64(r, u);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_mulhi_epi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  //   const int16_t *_b = (const int16_t *)impl.test_cases_int_pointer2;
  //   int16_t d[8];
  //   for (uint32_t i = 0; i < 8; i++) {
  //     int32_t m = (int32_t)_a[i] * (int32_t)_b[i];
  //     d[i] = (int16_t)(m >> 16);
  //   }
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //   __m128i c = _mm_mulhi_epi16(a, b);
  //   return VALIDATE_INT16_M128(c, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_mulhi_epu16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const uint16_t *_a = (const uint16_t *)impl.test_cases_int_pointer1;
  //   const uint16_t *_b = (const uint16_t *)impl.test_cases_int_pointer2;
  //   uint16_t d[8];
  //   for (uint32_t i = 0; i < 8; i++) {
  //     uint32_t m = (uint32_t)_a[i] * (uint32_t)_b[i];
  //     d[i] = (uint16_t)(m >> 16);
  //   }
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //   __m128i c = _mm_mulhi_epu16(a, b);
  //   return VALIDATE_INT16_M128(c, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_mullo_epi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  //   const int16_t *_b = (const int16_t *)impl.test_cases_int_pointer2;
  //   int16_t d[8];
  //   d[0] = _a[0] * _b[0];
  //   d[1] = _a[1] * _b[1];
  //   d[2] = _a[2] * _b[2];
  //   d[3] = _a[3] * _b[3];
  //   d[4] = _a[4] * _b[4];
  //   d[5] = _a[5] * _b[5];
  //   d[6] = _a[6] * _b[6];
  //   d[7] = _a[7] * _b[7];
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //   __m128i c = _mm_mullo_epi16(a, b);
  //   return VALIDATE_INT16_M128(c, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_or_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int64_t *_a = (const int64_t *)impl.test_cases_float_pointer1;
  //   const int64_t *_b = (const int64_t *)impl.test_cases_float_pointer2;
  //
  //   int64_t d0 = _a[0] | _b[0];
  //   int64_t d1 = _a[1] | _b[1];
  //
  //   __m128d a = load_m128d(_a);
  //   __m128d b = load_m128d(_b);
  //   __m128d c = _mm_or_pd(a, b);
  //
  //   return validate_double(c, *((double *)&d0), *((double *)&d1));
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_or_si128(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int32_t *_a = impl.test_cases_int_pointer1;
  //   const int32_t *_b = impl.test_cases_int_pointer2;
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //   __m128 fc = _mm_or_ps(*(const __m128 *)&a, *(const __m128 *)&b);
  //   __m128i c = *(const __m128i *)&fc;
  // now for the assertion...
  //   const uint32_t *ia = (const uint32_t *)&a;
  //   const uint32_t *ib = (const uint32_t *)&b;
  //   uint32_t r[4];
  //   r[0] = ia[0] | ib[0];
  //   r[1] = ia[1] | ib[1];
  //   r[2] = ia[2] | ib[2];
  //   r[3] = ia[3] | ib[3];
  //   __m128i ret = do_mm_set_epi32(r[3], r[2], r[1], r[0]);
  //   result_t res = VALIDATE_INT32_M128(c, r);
  //   if (res) {
  //     res = VALIDATE_INT32_M128(ret, r);
  //   }
  //   return res;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_packs_epi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   int8_t max = INT8_MAX;
  //   int8_t min = INT8_MIN;
  //   const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  //   const int16_t *_b = (const int16_t *)impl.test_cases_int_pointer2;
  //
  //   int8_t d[16];
  //   for (int i = 0; i < 8; i++) {
  //     if (_a[i] > max)
  //       d[i] = max;
  //     else if (_a[i] < min)
  //       d[i] = min;
  //     else
  //       d[i] = (int8_t)_a[i];
  //   }
  //   for (int i = 0; i < 8; i++) {
  //     if (_b[i] > max)
  //       d[i + 8] = max;
  //     else if (_b[i] < min)
  //       d[i + 8] = min;
  //     else
  //       d[i + 8] = (int8_t)_b[i];
  //   }
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //   __m128i c = _mm_packs_epi16(a, b);
  //
  //   return VALIDATE_INT8_M128(c, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_packs_epi32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   int16_t max = INT16_MAX;
  //   int16_t min = INT16_MIN;
  //   const int32_t *_a = (const int32_t *)impl.test_cases_int_pointer1;
  //   const int32_t *_b = (const int32_t *)impl.test_cases_int_pointer2;
  //
  //   int16_t d[8];
  //   for (int i = 0; i < 4; i++) {
  //     if (_a[i] > max)
  //       d[i] = max;
  //     else if (_a[i] < min)
  //       d[i] = min;
  //     else
  //       d[i] = (int16_t)_a[i];
  //   }
  //   for (int i = 0; i < 4; i++) {
  //     if (_b[i] > max)
  //       d[i + 4] = max;
  //     else if (_b[i] < min)
  //       d[i + 4] = min;
  //     else
  //       d[i + 4] = (int16_t)_b[i];
  //   }
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //   __m128i c = _mm_packs_epi32(a, b);
  //
  //   return VALIDATE_INT16_M128(c, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_packus_epi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   uint8_t max = UINT8_MAX;
  //   uint8_t min = 0;
  //   const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  //   const int16_t *_b = (const int16_t *)impl.test_cases_int_pointer2;
  //
  //   uint8_t d[16];
  //   for (int i = 0; i < 8; i++) {
  //     if (_a[i] > (int16_t)max)
  //       d[i] = max;
  //     else if (_a[i] < (int16_t)min)
  //       d[i] = min;
  //     else
  //       d[i] = (uint8_t)_a[i];
  //   }
  //   for (int i = 0; i < 8; i++) {
  //     if (_b[i] > (int16_t)max)
  //       d[i + 8] = max;
  //     else if (_b[i] < (int16_t)min)
  //       d[i + 8] = min;
  //     else
  //       d[i + 8] = (uint8_t)_b[i];
  //   }
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //   __m128i c = _mm_packus_epi16(a, b);
  //
  //   return VALIDATE_UINT8_M128(c, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_pause(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   _mm_pause();
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_sad_epu8(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const uint8_t *_a = (const uint8_t *)impl.test_cases_int_pointer1;
  //   const uint8_t *_b = (const uint8_t *)impl.test_cases_int_pointer2;
  //   uint16_t d0 = 0;
  //   uint16_t d1 = 0;
  //   for (int i = 0; i < 8; i++) {
  //     d0 += abs(_a[i] - _b[i]);
  //   }
  //   for (int i = 8; i < 16; i++) {
  //     d1 += abs(_a[i] - _b[i]);
  //   }
  //
  //   const __m128i a = load_m128i(_a);
  //   const __m128i b = load_m128i(_b);
  //   __m128i c = _mm_sad_epu8(a, b);
  //   return validate_uint16(c, d0, 0, 0, 0, d1, 0, 0, 0);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_set_epi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  int16_t d[8];
  d[0] = _a[0];
  d[1] = _a[1];
  d[2] = _a[2];
  d[3] = _a[3];
  d[4] = _a[4];
  d[5] = _a[5];
  d[6] = _a[6];
  d[7] = _a[7];

  __m128i c = _mm_set_epi16(d[7], d[6], d[5], d[4], d[3], d[2], d[1], d[0]);
  return VALIDATE_INT16_M128(c, d);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_set_epi32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  int32_t d[4];
  d[3] = impl.test_cases_ints[iter];
  d[2] = impl.test_cases_ints[iter + 1];
  d[1] = impl.test_cases_ints[iter + 2];
  d[0] = impl.test_cases_ints[iter + 3];
  __m128i a = _mm_set_epi32(d[3], d[2], d[1], d[0]);
  return VALIDATE_INT32_M128(a, d);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_set_epi64(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int64_t *_a = (const int64_t *)impl.test_cases_int_pointer1;

  __m128i ret = _mm_set_epi64(load_m64(&_a[1]), load_m64(&_a[0]));

  return validate_int64(ret, _a[0], _a[1]);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_set_epi64x(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int64_t *_a = (const int64_t *)impl.test_cases_int_pointer1;

  __m128i ret = _mm_set_epi64x(_a[1], _a[0]);

  return validate_int64(ret, _a[0], _a[1]);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_set_epi8(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int8_t *_a = (const int8_t *)impl.test_cases_int_pointer1;
  int8_t d[16];
  d[0] = _a[0];
  d[1] = _a[1];
  d[2] = _a[2];
  d[3] = _a[3];
  d[4] = _a[4];
  d[5] = _a[5];
  d[6] = _a[6];
  d[7] = _a[7];
  d[8] = _a[8];
  d[9] = _a[9];
  d[10] = _a[10];
  d[11] = _a[11];
  d[12] = _a[12];
  d[13] = _a[13];
  d[14] = _a[14];
  d[15] = _a[15];

  __m128i c = _mm_set_epi8(d[15], d[14], d[13], d[12], d[11], d[10], d[9], d[8],
                           d[7], d[6], d[5], d[4], d[3], d[2], d[1], d[0]);
  return VALIDATE_INT8_M128(c, d);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_set_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const double *p = (const double *)impl.test_cases_float_pointer1;
  double x = p[0];
  double y = p[1];
  __m128d a = _mm_set_pd(x, y);
  return validate_double(a, y, x);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_set_pd1(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const double _a = impl.test_cases_floats[iter];

  __m128d a = _mm_set_pd1(_a);

  return validate_double(a, _a, _a);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_set_sd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const double *_a = (const double *)impl.test_cases_float_pointer1;

  double f0 = _a[0];
  double f1 = 0.0;

  __m128d a = _mm_set_sd(_a[0]);
  return validate_double(a, f0, f1);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_set1_epi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  int16_t d0 = _a[0];

  __m128i c = _mm_set1_epi16(d0);
  return validate_int16(c, d0, d0, d0, d0, d0, d0, d0, d0);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_set1_epi32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  int32_t x = impl.test_cases_ints[iter];
  __m128i a = _mm_set1_epi32(x);
  return validate_int32(a, x, x, x, x);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_set1_epi64(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int64_t *_a = (const int64_t *)impl.test_cases_int_pointer1;

  __m128i ret = _mm_set1_epi64(load_m64(&_a[0]));

  return validate_int64(ret, _a[0], _a[0]);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_set1_epi64x(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int64_t *_a = (const int64_t *)impl.test_cases_int_pointer1;

  __m128i ret = _mm_set1_epi64x(_a[0]);

  return validate_int64(ret, _a[0], _a[0]);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_set1_epi8(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int8_t *_a = (const int8_t *)impl.test_cases_int_pointer1;
  int8_t d0 = _a[0];
  __m128i c = _mm_set1_epi8(d0);
  return validate_int8(c, d0, d0, d0, d0, d0, d0, d0, d0, d0, d0, d0, d0, d0,
                       d0, d0, d0);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_set1_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const double *_a = (const double *)impl.test_cases_float_pointer1;
  double d0 = _a[0];
  __m128d c = _mm_set1_pd(d0);
  return validate_double(c, d0, d0);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_setr_epi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;

  __m128i c =
      _mm_setr_epi16(_a[0], _a[1], _a[2], _a[3], _a[4], _a[5], _a[6], _a[7]);

  return VALIDATE_INT16_M128(c, _a);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_setr_epi32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int32_t *_a = (const int32_t *)impl.test_cases_int_pointer1;
  __m128i c = _mm_setr_epi32(_a[0], _a[1], _a[2], _a[3]);
  return VALIDATE_INT32_M128(c, _a);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_setr_epi64(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int64_t *_a = (const int64_t *)impl.test_cases_int_pointer1;
  __m128i c = _mm_setr_epi64(load_m64(&_a[0]), load_m64(&_a[1]));
  return validate_int64(c, _a[0], _a[1]);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_setr_epi8(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int8_t *_a = (const int8_t *)impl.test_cases_int_pointer1;

  __m128i c = _mm_setr_epi8(_a[0], _a[1], _a[2], _a[3], _a[4], _a[5], _a[6],
                            _a[7], _a[8], _a[9], _a[10], _a[11], _a[12], _a[13],
                            _a[14], _a[15]);

  return VALIDATE_INT8_M128(c, _a);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_setr_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const double *p = (const double *)impl.test_cases_float_pointer1;

  double x = p[0];
  double y = p[1];

  __m128d a = _mm_setr_pd(x, y);

  return validate_double(a, x, y);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_setzero_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   __m128d a = _mm_setzero_pd();
  //   return validate_double(a, 0, 0);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_setzero_si128(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   __m128i a = _mm_setzero_si128();
  //   return validate_int32(a, 0, 0, 0, 0);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_shuffle_epi32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int32_t *_a = impl.test_cases_int_pointer1;
  //   __m128i a, c;
  //
  // #define TEST_IMPL(IDX)
  //   int32_t d##IDX[4];
  //   d##IDX[0] = _a[((IDX) & 0x3)];
  //   d##IDX[1] = _a[((IDX >> 2) & 0x3)];
  //   d##IDX[2] = _a[((IDX >> 4) & 0x3)];
  //   d##IDX[3] = _a[((IDX >> 6) & 0x3)];
  //
  //   a = load_m128i(_a);
  //   c = _mm_shuffle_epi32(a, IDX);
  // CHECK_RESULT(VALIDATE_INT32_M128(c, d##IDX))
  //
  //   IMM_256_ITER
  // #undef TEST_IMPL
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_shuffle_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const double *_a = (const double *)impl.test_cases_float_pointer1;
  //   const double *_b = (const double *)impl.test_cases_float_pointer2;
  //   __m128d a, b, c;
  //
  // #define TEST_IMPL(IDX)
  //   a = load_m128d(_a);
  //   b = load_m128d(_b);
  //   c = _mm_shuffle_pd(a, b, IDX);
  //
  //   double d0##IDX = _a[IDX & 0x1];
  //   double d1##IDX = _b[(IDX & 0x2) >> 1];
  // CHECK_RESULT(validate_double(c, d0##IDX, d1##IDX))
  //
  //   IMM_4_ITER
  // #undef TEST_IMPL
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_shufflehi_epi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  //   __m128i a, c;
  //
  // #define TEST_IMPL(IDX)
  //   int16_t d##IDX[8];
  //   d##IDX[0] = _a[0];
  //   d##IDX[1] = _a[1];
  //   d##IDX[2] = _a[2];
  //   d##IDX[3] = _a[3];
  //   d##IDX[4] = ((const int64_t *)_a)[1] >> ((IDX & 0x3) * 16);
  //   d##IDX[5] = ((const int64_t *)_a)[1] >> (((IDX >> 2) & 0x3) * 16);
  //   d##IDX[6] = ((const int64_t *)_a)[1] >> (((IDX >> 4) & 0x3) * 16);
  //   d##IDX[7] = ((const int64_t *)_a)[1] >> (((IDX >> 6) & 0x3) * 16);
  //
  //   a = load_m128i(_a);
  //   c = _mm_shufflehi_epi16(a, IDX);
  //
  //   CHECK_RESULT(VALIDATE_INT16_M128(c, d##IDX))
  //
  //   IMM_256_ITER
  // #undef TEST_IMPL
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_shufflelo_epi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  //   __m128i a, c;
  //
  // #define TEST_IMPL(IDX)
  //   int16_t d##IDX[8];
  //   d##IDX[0] = ((const int64_t *)_a)[0] >> ((IDX & 0x3) * 16);
  //   d##IDX[1] = ((const int64_t *)_a)[0] >> (((IDX >> 2) & 0x3) * 16);
  //   d##IDX[2] = ((const int64_t *)_a)[0] >> (((IDX >> 4) & 0x3) * 16);
  //   d##IDX[3] = ((const int64_t *)_a)[0] >> (((IDX >> 6) & 0x3) * 16);
  //   d##IDX[4] = _a[4];
  //   d##IDX[5] = _a[5];
  //   d##IDX[6] = _a[6];
  //   d##IDX[7] = _a[7];
  //
  //   a = load_m128i(_a);
  //   c = _mm_shufflelo_epi16(a, IDX);
  //
  //   CHECK_RESULT(VALIDATE_INT16_M128(c, d##IDX))
  //
  //   IMM_256_ITER
  // #undef TEST_IMPL
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_sll_epi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  //   __m128i a, b, c;
  //
  // #define TEST_IMPL(IDX)
  //   uint16_t d##IDX[8];
  //   d##IDX[0] = (IDX > 15) ? 0 : _a[0] << IDX;
  //   d##IDX[1] = (IDX > 15) ? 0 : _a[1] << IDX;
  //   d##IDX[2] = (IDX > 15) ? 0 : _a[2] << IDX;
  //   d##IDX[3] = (IDX > 15) ? 0 : _a[3] << IDX;
  //   d##IDX[4] = (IDX > 15) ? 0 : _a[4] << IDX;
  //   d##IDX[5] = (IDX > 15) ? 0 : _a[5] << IDX;
  //   d##IDX[6] = (IDX > 15) ? 0 : _a[6] << IDX;
  //   d##IDX[7] = (IDX > 15) ? 0 : _a[7] << IDX;
  //
  //   a = load_m128i(_a);
  //   b = _mm_set1_epi64x(IDX);
  //   c = _mm_sll_epi16(a, b);
  // CHECK_RESULT(VALIDATE_INT16_M128(c, d##IDX))
  //
  //   IMM_64_ITER
  // #undef TEST_IMPL
  //
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_sll_epi32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int32_t *_a = (const int32_t *)impl.test_cases_int_pointer1;
  //   __m128i a, b, c;
  //
  // #define TEST_IMPL(IDX)
  //   uint32_t d##IDX[4];
  //   d##IDX[0] = (IDX > 31) ? 0 : _a[0] << IDX;
  //   d##IDX[1] = (IDX > 31) ? 0 : _a[1] << IDX;
  //   d##IDX[2] = (IDX > 31) ? 0 : _a[2] << IDX;
  //   d##IDX[3] = (IDX > 31) ? 0 : _a[3] << IDX;
  //
  //   a = load_m128i(_a);
  //   b = _mm_set1_epi64x(IDX);
  //   c = _mm_sll_epi32(a, b);
  // CHECK_RESULT(VALIDATE_INT32_M128(c, d##IDX))
  //
  //   IMM_64_ITER
  // #undef TEST_IMPL
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_sll_epi64(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int64_t *_a = (const int64_t *)impl.test_cases_int_pointer1;
  //   __m128i a, b, c;
  //
  // #define TEST_IMPL(IDX)
  //   uint64_t d0##IDX = (IDX & ~63) ? 0 : _a[0] << IDX;
  //   uint64_t d1##IDX = (IDX & ~63) ? 0 : _a[1] << IDX;
  //
  //   a = load_m128i(_a);
  //   b = _mm_set1_epi64x(IDX);
  //   c = _mm_sll_epi64(a, b);
  //
  //   CHECK_RESULT(validate_int64(c, d0##IDX, d1##IDX))
  //
  //   IMM_64_ITER
  // #undef TEST_IMPL
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_slli_epi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  //   __m128i a, c;
  //
  // #define TEST_IMPL(IDX)
  //   int16_t d##IDX[8];
  //   d##IDX[0] = (IDX > 15) ? 0 : _a[0] << IDX;
  //   d##IDX[1] = (IDX > 15) ? 0 : _a[1] << IDX;
  //   d##IDX[2] = (IDX > 15) ? 0 : _a[2] << IDX;
  //   d##IDX[3] = (IDX > 15) ? 0 : _a[3] << IDX;
  //   d##IDX[4] = (IDX > 15) ? 0 : _a[4] << IDX;
  //   d##IDX[5] = (IDX > 15) ? 0 : _a[5] << IDX;
  //   d##IDX[6] = (IDX > 15) ? 0 : _a[6] << IDX;
  //   d##IDX[7] = (IDX > 15) ? 0 : _a[7] << IDX;
  //
  //   a = load_m128i(_a);
  //   c = _mm_slli_epi16(a, IDX);
  // CHECK_RESULT(VALIDATE_INT16_M128(c, d##IDX))
  //
  //   IMM_64_ITER
  // #undef TEST_IMPL
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_slli_epi32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int32_t *_a = (const int32_t *)impl.test_cases_int_pointer1;
  // #if defined(__clang__)
  // Clang compiler does not allow the second argument of _mm_slli_epi32() to
  // be greater than 31.
  //   const int count = (int)(iter % 33 - 1); // range: -1 ~ 31
  // #else
  //   const int count = (int)(iter % 34 - 1); // range: -1 ~ 32
  // #endif
  //
  //   int32_t d[4];
  //   d[0] = (count & ~31) ? 0 : _a[0] << count;
  //   d[1] = (count & ~31) ? 0 : _a[1] << count;
  //   d[2] = (count & ~31) ? 0 : _a[2] << count;
  //   d[3] = (count & ~31) ? 0 : _a[3] << count;
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i c = _mm_slli_epi32(a, count);
  //   return VALIDATE_INT32_M128(c, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_slli_epi64(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int64_t *_a = (const int64_t *)impl.test_cases_int_pointer1;
  // #if defined(__clang__)
  // Clang compiler does not allow the second argument of "_mm_slli_epi64()"
  // to be greater than 63.
  //   const int count = (int)(iter % 65 - 1); // range: -1 ~ 63
  // #else
  //   const int count = (int)(iter % 66 - 1); // range: -1 ~ 64
  // #endif
  //   int64_t d0 = (count & ~63) ? 0 : _a[0] << count;
  //   int64_t d1 = (count & ~63) ? 0 : _a[1] << count;
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i c = _mm_slli_epi64(a, count);
  //   return validate_int64(c, d0, d1);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_slli_si128(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int32_t *_a = impl.test_cases_int_pointer1;
  //
  //   int8_t d[16];
  //   int count = (iter % 5) << 2;
  //   for (int i = 0; i < 16; i++) {
  //     if (i < count)
  //       d[i] = 0;
  //     else
  //       d[i] = ((const int8_t *)_a)[i - count];
  //   }
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i ret;
  //   switch (iter % 5) {
  //   case 0:
  //     ret = _mm_slli_si128(a, 0);
  //     break;
  //   case 1:
  //     ret = _mm_slli_si128(a, 4);
  //     break;
  //   case 2:
  //     ret = _mm_slli_si128(a, 8);
  //     break;
  //   case 3:
  //     ret = _mm_slli_si128(a, 12);
  //     break;
  //   case 4:
  //     ret = _mm_slli_si128(a, 16);
  //     break;
  //   }
  //
  //   return VALIDATE_INT8_M128(ret, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_sqrt_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const double *_a = (const double *)impl.test_cases_float_pointer1;
  //
  //   double f0 = sqrt(_a[0]);
  //   double f1 = sqrt(_a[1]);
  //
  //   __m128d a = load_m128d(_a);
  //   __m128d c = _mm_sqrt_pd(a);
  //
  //   return validate_double_error(c, f0, f1, 1.0e-15);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_sqrt_sd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const double *_a = (const double *)impl.test_cases_float_pointer1;
  //   const double *_b = (const double *)impl.test_cases_float_pointer2;
  //
  //   double f0 = sqrt(_b[0]);
  //   double f1 = _a[1];
  //
  //   __m128d a = load_m128d(_a);
  //   __m128d b = load_m128d(_b);
  //   __m128d c = _mm_sqrt_sd(a, b);
  //
  //   return validate_double_error(c, f0, f1, 1.0e-15);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_sra_epi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  //   const int64_t count = (int64_t)(iter % 18 - 1); // range: -1 ~ 16
  //
  //   int16_t d[8];
  //   d[0] = (count & ~15) ? (_a[0] < 0 ? ~UINT16_C(0) : 0) : (_a[0] >> count);
  //   d[1] = (count & ~15) ? (_a[1] < 0 ? ~UINT16_C(0) : 0) : (_a[1] >> count);
  //   d[2] = (count & ~15) ? (_a[2] < 0 ? ~UINT16_C(0) : 0) : (_a[2] >> count);
  //   d[3] = (count & ~15) ? (_a[3] < 0 ? ~UINT16_C(0) : 0) : (_a[3] >> count);
  //   d[4] = (count & ~15) ? (_a[4] < 0 ? ~UINT16_C(0) : 0) : (_a[4] >> count);
  //   d[5] = (count & ~15) ? (_a[5] < 0 ? ~UINT16_C(0) : 0) : (_a[5] >> count);
  //   d[6] = (count & ~15) ? (_a[6] < 0 ? ~UINT16_C(0) : 0) : (_a[6] >> count);
  //   d[7] = (count & ~15) ? (_a[7] < 0 ? ~UINT16_C(0) : 0) : (_a[7] >> count);
  //
  //   __m128i a = _mm_load_si128((const __m128i *)_a);
  //   __m128i b = _mm_set1_epi64x(count);
  //   __m128i c = _mm_sra_epi16(a, b);
  //
  //   return VALIDATE_INT16_M128(c, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_sra_epi32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int32_t *_a = (const int32_t *)impl.test_cases_int_pointer1;
  //   const int64_t count = (int64_t)(iter % 34 - 1); // range: -1 ~ 32
  //
  //   int32_t d[4];
  //   d[0] = (count & ~31) ? (_a[0] < 0 ? ~UINT32_C(0) : 0) : _a[0] >> count;
  //   d[1] = (count & ~31) ? (_a[1] < 0 ? ~UINT32_C(0) : 0) : _a[1] >> count;
  //   d[2] = (count & ~31) ? (_a[2] < 0 ? ~UINT32_C(0) : 0) : _a[2] >> count;
  //   d[3] = (count & ~31) ? (_a[3] < 0 ? ~UINT32_C(0) : 0) : _a[3] >> count;
  //
  //   __m128i a = _mm_load_si128((const __m128i *)_a);
  //   __m128i b = _mm_set1_epi64x(count);
  //   __m128i c = _mm_sra_epi32(a, b);
  //
  //   return VALIDATE_INT32_M128(c, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_srai_epi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  //   const int32_t b = (int32_t)(iter % 18 - 1); // range: -1 ~ 16
  //   int16_t d[8];
  //   int count = (b & ~15) ? 15 : b;
  //
  //   for (int i = 0; i < 8; i++) {
  //     d[i] = _a[i] >> count;
  //   }
  //
  //   __m128i a = _mm_load_si128((const __m128i *)_a);
  //   __m128i c = _mm_srai_epi16(a, b);
  //
  //   return VALIDATE_INT16_M128(c, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_srai_epi32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int32_t *_a = (const int32_t *)impl.test_cases_int_pointer1;
  //   const int32_t b = (int32_t)(iter % 34 - 1); // range: -1 ~ 32
  //
  //   int32_t d[4];
  //   int count = (b & ~31) ? 31 : b;
  //   for (int i = 0; i < 4; i++) {
  //     d[i] = _a[i] >> count;
  //   }
  //
  //   __m128i a = _mm_load_si128((const __m128i *)_a);
  //   __m128i c = _mm_srai_epi32(a, b);
  //
  //   return VALIDATE_INT32_M128(c, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_srl_epi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  //   const int64_t count = (int64_t)(iter % 18 - 1); // range: -1 ~ 16
  //
  //   uint16_t d[8];
  //   d[0] = (count & ~15) ? 0 : (uint16_t)(_a[0]) >> count;
  //   d[1] = (count & ~15) ? 0 : (uint16_t)(_a[1]) >> count;
  //   d[2] = (count & ~15) ? 0 : (uint16_t)(_a[2]) >> count;
  //   d[3] = (count & ~15) ? 0 : (uint16_t)(_a[3]) >> count;
  //   d[4] = (count & ~15) ? 0 : (uint16_t)(_a[4]) >> count;
  //   d[5] = (count & ~15) ? 0 : (uint16_t)(_a[5]) >> count;
  //   d[6] = (count & ~15) ? 0 : (uint16_t)(_a[6]) >> count;
  //   d[7] = (count & ~15) ? 0 : (uint16_t)(_a[7]) >> count;
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i b = _mm_set1_epi64x(count);
  //   __m128i c = _mm_srl_epi16(a, b);
  //
  //   return VALIDATE_INT16_M128(c, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_srl_epi32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int32_t *_a = (const int32_t *)impl.test_cases_int_pointer1;
  //   const int64_t count = (int64_t)(iter % 34 - 1); // range: -1 ~ 32
  //
  //   uint32_t d[4];
  //   d[0] = (count & ~31) ? 0 : (uint32_t)(_a[0]) >> count;
  //   d[1] = (count & ~31) ? 0 : (uint32_t)(_a[1]) >> count;
  //   d[2] = (count & ~31) ? 0 : (uint32_t)(_a[2]) >> count;
  //   d[3] = (count & ~31) ? 0 : (uint32_t)(_a[3]) >> count;
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i b = _mm_set1_epi64x(count);
  //   __m128i c = _mm_srl_epi32(a, b);
  //
  //   return VALIDATE_INT32_M128(c, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_srl_epi64(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int64_t *_a = (const int64_t *)impl.test_cases_int_pointer1;
  //   const int64_t count = (int64_t)(iter % 66 - 1); // range: -1 ~ 64
  //
  //   uint64_t d0 = (count & ~63) ? 0 : (uint64_t)(_a[0]) >> count;
  //   uint64_t d1 = (count & ~63) ? 0 : (uint64_t)(_a[1]) >> count;
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i b = _mm_set1_epi64x(count);
  //   __m128i c = _mm_srl_epi64(a, b);
  //
  //   return validate_int64(c, d0, d1);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_srli_epi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  //   const int count = (int)(iter % 18 - 1); // range: -1 ~ 16
  //
  //   int16_t d[8];
  //   d[0] = count & (~15) ? 0 : (uint16_t)(_a[0]) >> count;
  //   d[1] = count & (~15) ? 0 : (uint16_t)(_a[1]) >> count;
  //   d[2] = count & (~15) ? 0 : (uint16_t)(_a[2]) >> count;
  //   d[3] = count & (~15) ? 0 : (uint16_t)(_a[3]) >> count;
  //   d[4] = count & (~15) ? 0 : (uint16_t)(_a[4]) >> count;
  //   d[5] = count & (~15) ? 0 : (uint16_t)(_a[5]) >> count;
  //   d[6] = count & (~15) ? 0 : (uint16_t)(_a[6]) >> count;
  //   d[7] = count & (~15) ? 0 : (uint16_t)(_a[7]) >> count;
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i c = _mm_srli_epi16(a, count);
  //
  //   return VALIDATE_INT16_M128(c, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_srli_epi32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int32_t *_a = (const int32_t *)impl.test_cases_int_pointer1;
  //   const int count = (int)(iter % 34 - 1); // range: -1 ~ 32
  //
  //   int32_t d[4];
  //   d[0] = count & (~31) ? 0 : (uint32_t)(_a[0]) >> count;
  //   d[1] = count & (~31) ? 0 : (uint32_t)(_a[1]) >> count;
  //   d[2] = count & (~31) ? 0 : (uint32_t)(_a[2]) >> count;
  //   d[3] = count & (~31) ? 0 : (uint32_t)(_a[3]) >> count;
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i c = _mm_srli_epi32(a, count);
  //
  //   return VALIDATE_INT32_M128(c, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_srli_epi64(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int64_t *_a = (const int64_t *)impl.test_cases_int_pointer1;
  //   const int count = (int)(iter % 66 - 1); // range: -1 ~ 64
  //
  //   int64_t d0 = count & (~63) ? 0 : (uint64_t)(_a[0]) >> count;
  //   int64_t d1 = count & (~63) ? 0 : (uint64_t)(_a[1]) >> count;
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i c = _mm_srli_epi64(a, count);
  //
  //   return validate_int64(c, d0, d1);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_srli_si128(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int8_t *_a = (const int8_t *)impl.test_cases_int_pointer1;
  //   const int count = (iter % 5) << 2;
  //
  //   int8_t d[16];
  //   for (int i = 0; i < 16; i++) {
  //     if (i >= (16 - count))
  //       d[i] = 0;
  //     else
  //       d[i] = _a[i + count];
  //   }
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i ret;
  //   switch (iter % 5) {
  //   case 0:
  //     ret = _mm_srli_si128(a, 0);
  //     break;
  //   case 1:
  //     ret = _mm_srli_si128(a, 4);
  //     break;
  //   case 2:
  //     ret = _mm_srli_si128(a, 8);
  //     break;
  //   case 3:
  //     ret = _mm_srli_si128(a, 12);
  //     break;
  //   case 4:
  //     ret = _mm_srli_si128(a, 16);
  //     break;
  //   }
  //
  //   return VALIDATE_INT8_M128(ret, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_store_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  double *p = (double *)impl.test_cases_float_pointer1;
  double x = impl.test_cases_floats[iter + 4];
  double y = impl.test_cases_floats[iter + 6];

  __m128d a = _mm_set_pd(x, y);
  _mm_store_pd(p, a);
  ASSERT_RETURN(p[0] == y);
  ASSERT_RETURN(p[1] == x);
  return TEST_SUCCESS;
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_store_pd1(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  double *p = (double *)impl.test_cases_float_pointer1;
  double _a[2] = {(double)impl.test_cases_floats[iter],
                  (double)impl.test_cases_floats[iter + 1]};

  __m128d a = load_m128d(_a);
  _mm_store_pd1(p, a);
  ASSERT_RETURN(p[0] == impl.test_cases_floats[iter]);
  ASSERT_RETURN(p[1] == impl.test_cases_floats[iter]);
  return TEST_SUCCESS;
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_store_sd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  double *p = (double *)impl.test_cases_float_pointer1;
  double _a[2] = {(double)impl.test_cases_floats[iter],
                  (double)impl.test_cases_floats[iter + 1]};

  __m128d a = load_m128d(_a);
  _mm_store_sd(p, a);
  ASSERT_RETURN(p[0] == impl.test_cases_floats[iter]);
  return TEST_SUCCESS;
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_store_si128(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int32_t *_a = (const int32_t *)impl.test_cases_int_pointer1;
  alignas(16) int32_t p[4];

  __m128i a = load_m128i(_a);
  _mm_store_si128((__m128i *)p, a);

  return VALIDATE_INT32_M128(a, p);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_store1_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  return test_mm_store_pd1(impl, iter);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_storeh_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   double *p = (double *)impl.test_cases_float_pointer1;
  //   double mem;
  //
  //   __m128d a = load_m128d(p);
  //   _mm_storeh_pd(&mem, a);
  //
  //   ASSERT_RETURN(mem == p[1]);
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_storel_epi64(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   int64_t *p = (int64_t *)impl.test_cases_int_pointer1;
  //   __m128i mem;
  //
  //   __m128i a = load_m128i(p);
  //   _mm_storel_epi64(&mem, a);
  //
  //   ASSERT_RETURN(((SIMDVec *)&mem)->m128_u64[0] == (uint64_t)p[0]);
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_storel_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   double *p = (double *)impl.test_cases_float_pointer1;
  //   double mem;
  //
  //   __m128d a = load_m128d(p);
  //   _mm_storel_pd(&mem, a);
  //
  //   ASSERT_RETURN(mem == p[0]);
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_storer_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   double *p = (double *)impl.test_cases_float_pointer1;
  //   double mem[2];
  //
  //   __m128d a = load_m128d(p);
  //   _mm_storer_pd(mem, a);
  //
  //   __m128d res = load_m128d(mem);
  //   return validate_double(res, p[1], p[0]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_storeu_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   double *p = (double *)impl.test_cases_float_pointer1;
  //   double x = impl.test_cases_floats[iter + 4];
  //   double y = impl.test_cases_floats[iter + 6];
  //
  //   __m128d a = _mm_set_pd(x, y);
  //   _mm_storeu_pd(p, a);
  //   ASSERT_RETURN(p[0] == y);
  //   ASSERT_RETURN(p[1] == x);
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_storeu_si128(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int32_t *_a = (const int32_t *)impl.test_cases_int_pointer1;
  //   __m128i b;
  //   __m128i a = load_m128i(_a);
  //   _mm_storeu_si128(&b, a);
  //   int32_t *_b = (int32_t *)&b;
  //   return VALIDATE_INT32_M128(a, _b);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_storeu_si32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  // The GCC version before 11 does not implement intrinsic function
  // _mm_storeu_si32. Check https://gcc.gnu.org/bugzilla/show_bug.cgi?id=95483
  // for more information.
  // #if (defined(__GNUC__) && !defined(__clang__)) && (__GNUC__ <= 10)
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
  // #else
  //   const int32_t *_a = (const int32_t *)impl.test_cases_int_pointer1;
  //   __m128i b;
  //   __m128i a = load_m128i(_a);
  //   _mm_storeu_si32(&b, a);
  //   int32_t *_b = (int32_t *)&b;
  //   return validate_int32(b, _a[0], _b[1], _b[2], _b[3]);
  // #endif
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_stream_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const double *_a = (const double *)impl.test_cases_float_pointer1;
  //   double p[2];
  //
  //   __m128d a = load_m128d(_a);
  //   _mm_stream_pd(p, a);
  //
  //   return validate_double(a, p[0], p[1]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_stream_si128(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int32_t *_a = (const int32_t *)impl.test_cases_int_pointer1;
  //   alignas(16) int32_t p[4];
  //
  //   __m128i a = load_m128i(_a);
  //   _mm_stream_si128((__m128i *)p, a);
  //
  //   return VALIDATE_INT32_M128(a, p);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_stream_si32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int32_t a = (const int32_t)impl.test_cases_ints[iter];
  //   int32_t p;
  //
  //   _mm_stream_si32(&p, a);
  //
  //   ASSERT_RETURN(a == p)
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_stream_si64(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int64_t a = (const int64_t)impl.test_cases_ints[iter];
  //   __int64 p[1];
  //   _mm_stream_si64(p, a);
  //   ASSERT_RETURN(p[0] == a);
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_sub_epi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  const int16_t *_b = (const int16_t *)impl.test_cases_int_pointer2;
  int16_t d[8];
  d[0] = _a[0] - _b[0];
  d[1] = _a[1] - _b[1];
  d[2] = _a[2] - _b[2];
  d[3] = _a[3] - _b[3];
  d[4] = _a[4] - _b[4];
  d[5] = _a[5] - _b[5];
  d[6] = _a[6] - _b[6];
  d[7] = _a[7] - _b[7];

  __m128i a = load_m128i(_a);
  __m128i b = load_m128i(_b);
  __m128i c = _mm_sub_epi16(a, b);
  return VALIDATE_INT16_M128(c, d);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_sub_epi32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int32_t *_a = impl.test_cases_int_pointer1;
  const int32_t *_b = impl.test_cases_int_pointer2;
  int32_t d[4];
  d[0] = _a[0] - _b[0];
  d[1] = _a[1] - _b[1];
  d[2] = _a[2] - _b[2];
  d[3] = _a[3] - _b[3];

  __m128i a = load_m128i(_a);
  __m128i b = load_m128i(_b);
  __m128i c = _mm_sub_epi32(a, b);
  return VALIDATE_INT32_M128(c, d);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_sub_epi64(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int64_t *_a = (int64_t *)impl.test_cases_int_pointer1;
  const int64_t *_b = (int64_t *)impl.test_cases_int_pointer2;
  int64_t d0 = _a[0] - _b[0];
  int64_t d1 = _a[1] - _b[1];

  __m128i a = load_m128i(_a);
  __m128i b = load_m128i(_b);
  __m128i c = _mm_sub_epi64(a, b);
  return validate_int64(c, d0, d1);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_sub_epi8(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int8_t *_a = (const int8_t *)impl.test_cases_int_pointer1;
  const int8_t *_b = (const int8_t *)impl.test_cases_int_pointer2;
  int8_t d[16];
  d[0] = _a[0] - _b[0];
  d[1] = _a[1] - _b[1];
  d[2] = _a[2] - _b[2];
  d[3] = _a[3] - _b[3];
  d[4] = _a[4] - _b[4];
  d[5] = _a[5] - _b[5];
  d[6] = _a[6] - _b[6];
  d[7] = _a[7] - _b[7];
  d[8] = _a[8] - _b[8];
  d[9] = _a[9] - _b[9];
  d[10] = _a[10] - _b[10];
  d[11] = _a[11] - _b[11];
  d[12] = _a[12] - _b[12];
  d[13] = _a[13] - _b[13];
  d[14] = _a[14] - _b[14];
  d[15] = _a[15] - _b[15];

  __m128i a = load_m128i(_a);
  __m128i b = load_m128i(_b);
  __m128i c = _mm_sub_epi8(a, b);
  return VALIDATE_INT8_M128(c, d);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_sub_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const double *_a = (const double *)impl.test_cases_float_pointer1;
  const double *_b = (const double *)impl.test_cases_float_pointer2;
  double d0 = _a[0] - _b[0];
  double d1 = _a[1] - _b[1];

  __m128d a = load_m128d(_a);
  __m128d b = load_m128d(_b);
  __m128d c = _mm_sub_pd(a, b);
  return validate_double(c, d0, d1);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_sub_sd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const double *_a = (const double *)impl.test_cases_float_pointer1;
  const double *_b = (const double *)impl.test_cases_float_pointer2;
  double d0 = _a[0] - _b[0];
  double d1 = _a[1];

  __m128d a = load_m128d(_a);
  __m128d b = load_m128d(_b);
  __m128d c = _mm_sub_sd(a, b);
  return validate_double(c, d0, d1);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_sub_si64(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int64_t *_a = (const int64_t *)impl.test_cases_int_pointer1;
  const int64_t *_b = (const int64_t *)impl.test_cases_int_pointer2;

  int64_t d = _a[0] - _b[0];

  __m64 a = load_m64(_a);
  __m64 b = load_m64(_b);
  __m64 c = _mm_sub_si64(a, b);

  return validate_int64(c, d);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_subs_epi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   int32_t max = 32767;
  //   int32_t min = -32768;
  //   const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  //   const int16_t *_b = (const int16_t *)impl.test_cases_int_pointer2;
  //
  //   int16_t d[8];
  //   for (int i = 0; i < 8; i++) {
  //     int32_t res = (int32_t)_a[i] - (int32_t)_b[i];
  //     if (res > max)
  //       d[i] = max;
  //     else if (res < min)
  //       d[i] = min;
  //     else
  //       d[i] = (int16_t)res;
  //   }
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //   __m128i c = _mm_subs_epi16(a, b);
  //
  //   return VALIDATE_INT16_M128(c, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_subs_epi8(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   int16_t max = 127;
  //   int16_t min = -128;
  //   const int8_t *_a = (const int8_t *)impl.test_cases_int_pointer1;
  //   const int8_t *_b = (const int8_t *)impl.test_cases_int_pointer2;
  //
  //   int8_t d[16];
  //   for (int i = 0; i < 16; i++) {
  //     int16_t res = (int16_t)_a[i] - (int16_t)_b[i];
  //     if (res > max)
  //       d[i] = max;
  //     else if (res < min)
  //       d[i] = min;
  //     else
  //       d[i] = (int8_t)res;
  //   }
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //   __m128i c = _mm_subs_epi8(a, b);
  //
  //   return VALIDATE_INT8_M128(c, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_subs_epu16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  //   const int16_t *_b = (const int16_t *)impl.test_cases_int_pointer2;
  //   uint16_t d[8];
  //   d[0] = (uint16_t)_a[0] - (uint16_t)_b[0];
  //   if (d[0] > (uint16_t)_a[0])
  //     d[0] = 0;
  //   d[1] = (uint16_t)_a[1] - (uint16_t)_b[1];
  //   if (d[1] > (uint16_t)_a[1])
  //     d[1] = 0;
  //   d[2] = (uint16_t)_a[2] - (uint16_t)_b[2];
  //   if (d[2] > (uint16_t)_a[2])
  //     d[2] = 0;
  //   d[3] = (uint16_t)_a[3] - (uint16_t)_b[3];
  //   if (d[3] > (uint16_t)_a[3])
  //     d[3] = 0;
  //   d[4] = (uint16_t)_a[4] - (uint16_t)_b[4];
  //   if (d[4] > (uint16_t)_a[4])
  //     d[4] = 0;
  //   d[5] = (uint16_t)_a[5] - (uint16_t)_b[5];
  //   if (d[5] > (uint16_t)_a[5])
  //     d[5] = 0;
  //   d[6] = (uint16_t)_a[6] - (uint16_t)_b[6];
  //   if (d[6] > (uint16_t)_a[6])
  //     d[6] = 0;
  //   d[7] = (uint16_t)_a[7] - (uint16_t)_b[7];
  //   if (d[7] > (uint16_t)_a[7])
  //     d[7] = 0;
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //
  //   __m128i c = _mm_subs_epu16(a, b);
  //   return VALIDATE_INT16_M128(c, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_subs_epu8(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int8_t *_a = (const int8_t *)impl.test_cases_int_pointer1;
  //   const int8_t *_b = (const int8_t *)impl.test_cases_int_pointer2;
  //   uint8_t d[16];
  //   d[0] = (uint8_t)_a[0] - (uint8_t)_b[0];
  //   if (d[0] > (uint8_t)_a[0])
  //     d[0] = 0;
  //   d[1] = (uint8_t)_a[1] - (uint8_t)_b[1];
  //   if (d[1] > (uint8_t)_a[1])
  //     d[1] = 0;
  //   d[2] = (uint8_t)_a[2] - (uint8_t)_b[2];
  //   if (d[2] > (uint8_t)_a[2])
  //     d[2] = 0;
  //   d[3] = (uint8_t)_a[3] - (uint8_t)_b[3];
  //   if (d[3] > (uint8_t)_a[3])
  //     d[3] = 0;
  //   d[4] = (uint8_t)_a[4] - (uint8_t)_b[4];
  //   if (d[4] > (uint8_t)_a[4])
  //     d[4] = 0;
  //   d[5] = (uint8_t)_a[5] - (uint8_t)_b[5];
  //   if (d[5] > (uint8_t)_a[5])
  //     d[5] = 0;
  //   d[6] = (uint8_t)_a[6] - (uint8_t)_b[6];
  //   if (d[6] > (uint8_t)_a[6])
  //     d[6] = 0;
  //   d[7] = (uint8_t)_a[7] - (uint8_t)_b[7];
  //   if (d[7] > (uint8_t)_a[7])
  //     d[7] = 0;
  //   d[8] = (uint8_t)_a[8] - (uint8_t)_b[8];
  //   if (d[8] > (uint8_t)_a[8])
  //     d[8] = 0;
  //   d[9] = (uint8_t)_a[9] - (uint8_t)_b[9];
  //   if (d[9] > (uint8_t)_a[9])
  //     d[9] = 0;
  //   d[10] = (uint8_t)_a[10] - (uint8_t)_b[10];
  //   if (d[10] > (uint8_t)_a[10])
  //     d[10] = 0;
  //   d[11] = (uint8_t)_a[11] - (uint8_t)_b[11];
  //   if (d[11] > (uint8_t)_a[11])
  //     d[11] = 0;
  //   d[12] = (uint8_t)_a[12] - (uint8_t)_b[12];
  //   if (d[12] > (uint8_t)_a[12])
  //     d[12] = 0;
  //   d[13] = (uint8_t)_a[13] - (uint8_t)_b[13];
  //   if (d[13] > (uint8_t)_a[13])
  //     d[13] = 0;
  //   d[14] = (uint8_t)_a[14] - (uint8_t)_b[14];
  //   if (d[14] > (uint8_t)_a[14])
  //     d[14] = 0;
  //   d[15] = (uint8_t)_a[15] - (uint8_t)_b[15];
  //   if (d[15] > (uint8_t)_a[15])
  //     d[15] = 0;
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //   __m128i c = _mm_subs_epu8(a, b);
  //   return VALIDATE_INT8_M128(c, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_ucomieq_sd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
#if defined(__GNUC__) && !defined(__clang__)
  return TEST_UNIMPL;
#else
  const double *_a = (const double *)impl.test_cases_float_pointer1;
  const double *_b = (const double *)impl.test_cases_float_pointer2;
  int32_t _c = (_a[0] == _b[0]) ? 1 : 0;

  __m128d a = load_m128d(_a);
  __m128d b = load_m128d(_b);
  int32_t c = _mm_ucomieq_sd(a, b);

  ASSERT_RETURN(c == _c);
  return TEST_SUCCESS;
#endif
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_ucomige_sd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const double *_a = (const double *)impl.test_cases_float_pointer1;
  const double *_b = (const double *)impl.test_cases_float_pointer2;
  int32_t _c = (_a[0] >= _b[0]) ? 1 : 0;

  __m128d a = load_m128d(_a);
  __m128d b = load_m128d(_b);
  int32_t c = _mm_ucomige_sd(a, b);

  ASSERT_RETURN(c == _c);
  return TEST_SUCCESS;
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_ucomigt_sd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const double *_a = (const double *)impl.test_cases_float_pointer1;
  const double *_b = (const double *)impl.test_cases_float_pointer2;
  int32_t _c = (_a[0] > _b[0]) ? 1 : 0;

  __m128d a = load_m128d(_a);
  __m128d b = load_m128d(_b);
  int32_t c = _mm_ucomigt_sd(a, b);

  ASSERT_RETURN(c == _c);
  return TEST_SUCCESS;
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_ucomile_sd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
#if defined(__GNUC__) && !defined(__clang__)
  return TEST_UNIMPL;
#else
  const double *_a = (const double *)impl.test_cases_float_pointer1;
  const double *_b = (const double *)impl.test_cases_float_pointer2;
  int32_t _c = (!isnan(_a[0]) && !isnan(_b[0]) && (_a[0] <= _b[0])) ? 1 : 0;

  __m128d a = load_m128d(_a);
  __m128d b = load_m128d(_b);
  int32_t c = _mm_ucomile_sd(a, b);

  ASSERT_RETURN(c == _c);
  return TEST_SUCCESS;
#endif
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_ucomilt_sd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
#if defined(__GNUC__) && !defined(__clang__)
  return TEST_UNIMPL;
#else
  const double *_a = (const double *)impl.test_cases_float_pointer1;
  const double *_b = (const double *)impl.test_cases_float_pointer2;
  int32_t _c = (_a[0] < _b[0]) ? 1 : 0;

  __m128d a = load_m128d(_a);
  __m128d b = load_m128d(_b);
  int32_t c = _mm_ucomilt_sd(a, b);

  ASSERT_RETURN(c == _c);
  return TEST_SUCCESS;
#endif
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_ucomineq_sd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
#if defined(__GNUC__) && !defined(__clang__)
  return TEST_UNIMPL;
#else
  const double *_a = (const double *)impl.test_cases_float_pointer1;
  const double *_b = (const double *)impl.test_cases_float_pointer2;
  int32_t _c = (_a[0] != _b[0]) ? 1 : 0;

  __m128d a = load_m128d(_a);
  __m128d b = load_m128d(_b);
  int32_t c = _mm_ucomineq_sd(a, b);

  ASSERT_RETURN(c == _c);
  return TEST_SUCCESS;
#endif
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_undefined_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   __m128d a = _mm_undefined_pd();
  //   a = _mm_xor_pd(a, a);
  //   return validate_double(a, 0, 0);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_undefined_si128(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   __m128i a = _mm_undefined_si128();
  //   a = _mm_xor_si128(a, a);
  //   return validate_int64(a, 0, 0);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_unpackhi_epi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  //   const int16_t *_b = (const int16_t *)impl.test_cases_int_pointer2;
  //
  //   int16_t d[8];
  //   d[0] = _a[4];
  //   d[1] = _b[4];
  //   d[2] = _a[5];
  //   d[3] = _b[5];
  //   d[4] = _a[6];
  //   d[5] = _b[6];
  //   d[6] = _a[7];
  //   d[7] = _b[7];
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //   __m128i ret = _mm_unpackhi_epi16(a, b);
  //
  //   return VALIDATE_INT16_M128(ret, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_unpackhi_epi32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int32_t *_a = (const int32_t *)impl.test_cases_int_pointer1;
  //   const int32_t *_b = (const int32_t *)impl.test_cases_int_pointer2;
  //
  //   int32_t d[4];
  //   d[0] = _a[2];
  //   d[1] = _b[2];
  //   d[2] = _a[3];
  //   d[3] = _b[3];
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //   __m128i ret = _mm_unpackhi_epi32(a, b);
  //
  //   return VALIDATE_INT32_M128(ret, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_unpackhi_epi64(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int64_t *_a = (const int64_t *)impl.test_cases_int_pointer1;
  //   const int64_t *_b = (const int64_t *)impl.test_cases_int_pointer2;
  //
  //   int64_t i0 = _a[1];
  //   int64_t i1 = _b[1];
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //   __m128i ret = _mm_unpackhi_epi64(a, b);
  //
  //   return validate_int64(ret, i0, i1);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_unpackhi_epi8(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int8_t *_a = (const int8_t *)impl.test_cases_int_pointer1;
  //   const int8_t *_b = (const int8_t *)impl.test_cases_int_pointer2;
  //
  //   int8_t d[16];
  //   d[0] = _a[8];
  //   d[1] = _b[8];
  //   d[2] = _a[9];
  //   d[3] = _b[9];
  //   d[4] = _a[10];
  //   d[5] = _b[10];
  //   d[6] = _a[11];
  //   d[7] = _b[11];
  //   d[8] = _a[12];
  //   d[9] = _b[12];
  //   d[10] = _a[13];
  //   d[11] = _b[13];
  //   d[12] = _a[14];
  //   d[13] = _b[14];
  //   d[14] = _a[15];
  //   d[15] = _b[15];
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //   __m128i ret = _mm_unpackhi_epi8(a, b);
  //
  //   return VALIDATE_INT8_M128(ret, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_unpackhi_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const double *_a = (const double *)impl.test_cases_float_pointer1;
  //   const double *_b = (const double *)impl.test_cases_float_pointer2;
  //
  //   __m128d a = load_m128d(_a);
  //   __m128d b = load_m128d(_b);
  //   __m128d ret = _mm_unpackhi_pd(a, b);
  //
  //   return validate_double(ret, _a[1], _b[1]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_unpacklo_epi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  //   const int16_t *_b = (const int16_t *)impl.test_cases_int_pointer2;
  //
  //   int16_t d[8];
  //   d[0] = _a[0];
  //   d[1] = _b[0];
  //   d[2] = _a[1];
  //   d[3] = _b[1];
  //   d[4] = _a[2];
  //   d[5] = _b[2];
  //   d[6] = _a[3];
  //   d[7] = _b[3];
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //   __m128i ret = _mm_unpacklo_epi16(a, b);
  //
  //   return VALIDATE_INT16_M128(ret, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_unpacklo_epi32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int32_t *_a = (const int32_t *)impl.test_cases_int_pointer1;
  //   const int32_t *_b = (const int32_t *)impl.test_cases_int_pointer2;
  //
  //   int32_t d[4];
  //   d[0] = _a[0];
  //   d[1] = _b[0];
  //   d[2] = _a[1];
  //   d[3] = _b[1];
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //   __m128i ret = _mm_unpacklo_epi32(a, b);
  //
  //   return VALIDATE_INT32_M128(ret, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_unpacklo_epi64(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int64_t *_a = (const int64_t *)impl.test_cases_int_pointer1;
  //   const int64_t *_b = (const int64_t *)impl.test_cases_int_pointer2;
  //
  //   int64_t i0 = _a[0];
  //   int64_t i1 = _b[0];
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //   __m128i ret = _mm_unpacklo_epi64(a, b);
  //
  //   return validate_int64(ret, i0, i1);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_unpacklo_epi8(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int8_t *_a = (const int8_t *)impl.test_cases_int_pointer1;
  //   const int8_t *_b = (const int8_t *)impl.test_cases_int_pointer2;
  //
  //   int8_t d[16];
  //   d[0] = _a[0];
  //   d[1] = _b[0];
  //   d[2] = _a[1];
  //   d[3] = _b[1];
  //   d[4] = _a[2];
  //   d[5] = _b[2];
  //   d[6] = _a[3];
  //   d[7] = _b[3];
  //   d[8] = _a[4];
  //   d[9] = _b[4];
  //   d[10] = _a[5];
  //   d[11] = _b[5];
  //   d[12] = _a[6];
  //   d[13] = _b[6];
  //   d[14] = _a[7];
  //   d[15] = _b[7];
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //   __m128i ret = _mm_unpacklo_epi8(a, b);
  //
  //   return VALIDATE_INT8_M128(ret, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_unpacklo_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const double *_a = (const double *)impl.test_cases_float_pointer1;
  //   const double *_b = (const double *)impl.test_cases_float_pointer2;
  //
  //   __m128d a = load_m128d(_a);
  //   __m128d b = load_m128d(_b);
  //   __m128d ret = _mm_unpacklo_pd(a, b);
  //
  //   return validate_double(ret, _a[0], _b[0]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_xor_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int64_t *_a = (const int64_t *)impl.test_cases_float_pointer1;
  //   const int64_t *_b = (const int64_t *)impl.test_cases_float_pointer2;
  //
  //   int64_t d0 = _a[0] ^ _b[0];
  //   int64_t d1 = _a[1] ^ _b[1];
  //
  //   __m128d a = load_m128d(_a);
  //   __m128d b = load_m128d(_b);
  //   __m128d c = _mm_xor_pd(a, b);
  //
  //   return validate_double(c, *((double *)&d0), *((double *)&d1));
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_xor_si128(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int64_t *_a = (const int64_t *)impl.test_cases_int_pointer1;
  //   const int64_t *_b = (const int64_t *)impl.test_cases_int_pointer2;
  //
  //   int64_t d0 = _a[0] ^ _b[0];
  //   int64_t d1 = _a[1] ^ _b[1];
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //   __m128i c = _mm_xor_si128(a, b);
  //
  //   return validate_int64(c, d0, d1);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

/* SSE3 */
result_t test_mm_addsub_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const double *_a = (const double *)impl.test_cases_float_pointer1;
  const double *_b = (const double *)impl.test_cases_float_pointer2;

  double d0 = _a[0] - _b[0];
  double d1 = _a[1] + _b[1];

  __m128d a = load_m128d(_a);
  __m128d b = load_m128d(_b);
  __m128d c = _mm_addsub_pd(a, b);

  return validate_double(c, d0, d1);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_addsub_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  // _MM_SET_ROUNDING_MODE(_MM_ROUND_NEAREST);
  const float *_a = impl.test_cases_float_pointer1;
  const float *_b = impl.test_cases_float_pointer2;

  float f0 = _a[0] - _b[0];
  float f1 = _a[1] + _b[1];
  float f2 = _a[2] - _b[2];
  float f3 = _a[3] + _b[3];

  __m128 a = load_m128(_a);
  __m128 b = load_m128(_b);
  __m128 c = _mm_addsub_ps(a, b);

  return validate_float(c, f0, f1, f2, f3);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_hadd_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const double *_a = (const double *)impl.test_cases_float_pointer1;
  //   const double *_b = (const double *)impl.test_cases_float_pointer2;
  //
  //   double f0 = _a[0] + _a[1];
  //   double f1 = _b[0] + _b[1];
  //
  //   __m128d a = load_m128d(_a);
  //   __m128d b = load_m128d(_b);
  //   __m128d c = _mm_hadd_pd(a, b);
  //
  //   return validate_double(c, f0, f1);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_hadd_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   _MM_SET_ROUNDING_MODE(_MM_ROUND_NEAREST);
  //   const float *_a = impl.test_cases_float_pointer1;
  //   const float *_b = impl.test_cases_float_pointer2;
  //
  //   float f0 = _a[0] + _a[1];
  //   float f1 = _a[2] + _a[3];
  //   float f2 = _b[0] + _b[1];
  //   float f3 = _b[2] + _b[3];
  //
  //   __m128 a = load_m128(_a);
  //   __m128 b = load_m128(_b);
  //   __m128 c = _mm_hadd_ps(a, b);
  //
  //   return validate_float(c, f0, f1, f2, f3);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_hsub_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const double *_a = (const double *)impl.test_cases_float_pointer1;
  //   const double *_b = (const double *)impl.test_cases_float_pointer2;
  //
  //   double f0 = _a[0] - _a[1];
  //   double f1 = _b[0] - _b[1];
  //
  //   __m128d a = load_m128d(_a);
  //   __m128d b = load_m128d(_b);
  //   __m128d c = _mm_hsub_pd(a, b);
  //
  //   return validate_double(c, f0, f1);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_hsub_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   _MM_SET_ROUNDING_MODE(_MM_ROUND_NEAREST);
  //   const float *_a = impl.test_cases_float_pointer1;
  //   const float *_b = impl.test_cases_float_pointer2;
  //
  //   float f0 = _a[0] - _a[1];
  //   float f1 = _a[2] - _a[3];
  //   float f2 = _b[0] - _b[1];
  //   float f3 = _b[2] - _b[3];
  //
  //   __m128 a = load_m128(_a);
  //   __m128 b = load_m128(_b);
  //   __m128 c = _mm_hsub_ps(a, b);
  //
  //   return validate_float(c, f0, f1, f2, f3);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_lddqu_si128(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   return test_mm_loadu_si128(impl, iter);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_loaddup_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const double *addr = (const double *)impl.test_cases_float_pointer1;
  //
  //   __m128d ret = _mm_loaddup_pd(addr);
  //
  //   return validate_double(ret, addr[0], addr[0]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_movedup_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const double *p = (const double *)impl.test_cases_float_pointer1;
  //   __m128d a = load_m128d(p);
  //   __m128d b = _mm_movedup_pd(a);
  //
  //   return validate_double(b, p[0], p[0]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_movehdup_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *p = impl.test_cases_float_pointer1;
  //   __m128 a = load_m128(p);
  //   return validate_float(_mm_movehdup_ps(a), p[1], p[1], p[3], p[3]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_moveldup_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *p = impl.test_cases_float_pointer1;
  //   __m128 a = load_m128(p);
  //   return validate_float(_mm_moveldup_ps(a), p[0], p[0], p[2], p[2]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

/* SSSE3 */
result_t test_mm_abs_epi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  __m128i a = load_m128i(_a);
  __m128i c = _mm_abs_epi16(a);

  uint32_t d[8];
  d[0] = (_a[0] < 0) ? -_a[0] : _a[0];
  d[1] = (_a[1] < 0) ? -_a[1] : _a[1];
  d[2] = (_a[2] < 0) ? -_a[2] : _a[2];
  d[3] = (_a[3] < 0) ? -_a[3] : _a[3];
  d[4] = (_a[4] < 0) ? -_a[4] : _a[4];
  d[5] = (_a[5] < 0) ? -_a[5] : _a[5];
  d[6] = (_a[6] < 0) ? -_a[6] : _a[6];
  d[7] = (_a[7] < 0) ? -_a[7] : _a[7];

  return VALIDATE_UINT16_M128(c, d);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_abs_epi32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int32_t *_a = (const int32_t *)impl.test_cases_int_pointer1;
  __m128i a = load_m128i(_a);
  __m128i c = _mm_abs_epi32(a);

  uint32_t d[4];
  d[0] = (_a[0] < 0) ? -_a[0] : _a[0];
  d[1] = (_a[1] < 0) ? -_a[1] : _a[1];
  d[2] = (_a[2] < 0) ? -_a[2] : _a[2];
  d[3] = (_a[3] < 0) ? -_a[3] : _a[3];

  return VALIDATE_UINT32_M128(c, d);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_abs_epi8(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int8_t *_a = (const int8_t *)impl.test_cases_int_pointer1;
  __m128i a = load_m128i(_a);
  __m128i c = _mm_abs_epi8(a);

  uint32_t d[16];
  for (int i = 0; i < 16; i++) {
    d[i] = (_a[i] < 0) ? -_a[i] : _a[i];
  }

  return VALIDATE_UINT8_M128(c, d);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_abs_pi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  __m64 a = load_m64(_a);
  __m64 c = _mm_abs_pi16(a);

  uint32_t d[4];
  d[0] = (_a[0] < 0) ? -_a[0] : _a[0];
  d[1] = (_a[1] < 0) ? -_a[1] : _a[1];
  d[2] = (_a[2] < 0) ? -_a[2] : _a[2];
  d[3] = (_a[3] < 0) ? -_a[3] : _a[3];

  return VALIDATE_UINT16_M64(c, d);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_abs_pi32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int32_t *_a = (const int32_t *)impl.test_cases_int_pointer1;
  __m64 a = load_m64(_a);
  __m64 c = _mm_abs_pi32(a);

  uint32_t d[2];
  d[0] = (_a[0] < 0) ? -_a[0] : _a[0];
  d[1] = (_a[1] < 0) ? -_a[1] : _a[1];

  return VALIDATE_UINT32_M64(c, d);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_abs_pi8(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int8_t *_a = (const int8_t *)impl.test_cases_int_pointer1;
  __m64 a = load_m64(_a);
  __m64 c = _mm_abs_pi8(a);

  uint32_t d[8];
  d[0] = (_a[0] < 0) ? -_a[0] : _a[0];
  d[1] = (_a[1] < 0) ? -_a[1] : _a[1];
  d[2] = (_a[2] < 0) ? -_a[2] : _a[2];
  d[3] = (_a[3] < 0) ? -_a[3] : _a[3];
  d[4] = (_a[4] < 0) ? -_a[4] : _a[4];
  d[5] = (_a[5] < 0) ? -_a[5] : _a[5];
  d[6] = (_a[6] < 0) ? -_a[6] : _a[6];
  d[7] = (_a[7] < 0) ? -_a[7] : _a[7];

  return VALIDATE_UINT8_M64(c, d);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_alignr_epi8(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  // #if defined(__clang__)
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
  // #else
  //   const uint8_t *_a = (const uint8_t *)impl.test_cases_int_pointer1;
  //   const uint8_t *_b = (const uint8_t *)impl.test_cases_int_pointer2;
  //   unsigned int shift = (iter % 5) << 3;
  //   uint8_t d[32];
  //
  //   if (shift >= 32) {
  //     memset((void *)d, 0, sizeof(d));
  //   } else {
  //     memcpy((void *)d, (const void *)_b, 16);
  //     memcpy((void *)(d + 16), (const void *)_a, 16);
  //     // shifting
  //     for (size_t x = 0; x < sizeof(d); x++) {
  //       if (x + shift >= sizeof(d))
  //         d[x] = 0;
  //       else
  //         d[x] = d[x + shift];
  //     }
  //   }
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //   __m128i ret;
  //   switch (iter % 5) {
  //   case 0:
  //     ret = _mm_alignr_epi8(a, b, 0);
  //     break;
  //   case 1:
  //     ret = _mm_alignr_epi8(a, b, 8);
  //     break;
  //   case 2:
  //     ret = _mm_alignr_epi8(a, b, 16);
  //     break;
  //   case 3:
  //     ret = _mm_alignr_epi8(a, b, 24);
  //     break;
  //   case 4:
  //     ret = _mm_alignr_epi8(a, b, 32);
  //     break;
  //   }
  //
  //   return VALIDATE_UINT8_M128(ret, d);
  // #endif
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_alignr_pi8(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  // #if defined(__clang__)
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
  // #else
  //   const uint8_t *_a = (const uint8_t *)impl.test_cases_int_pointer1;
  //   const uint8_t *_b = (const uint8_t *)impl.test_cases_int_pointer2;
  //   unsigned int shift = (iter % 3) << 3;
  //   uint8_t d[16];
  //
  //   if (shift >= 16) {
  //     memset((void *)d, 0, sizeof(d));
  //   } else {
  //     memcpy((void *)d, (const void *)_b, 8);
  //     memcpy((void *)(d + 8), (const void *)_a, 8);
  //     // shifting
  //     for (size_t x = 0; x < sizeof(d); x++) {
  //       if (x + shift >= sizeof(d))
  //         d[x] = 0;
  //       else
  //         d[x] = d[x + shift];
  //     }
  //   }
  //
  //   __m64 a = load_m64(_a);
  //   __m64 b = load_m64(_b);
  //   __m64 ret;
  //   switch (iter % 3) {
  //   case 0:
  //     ret = _mm_alignr_pi8(a, b, 0);
  //     break;
  //   case 1:
  //     ret = _mm_alignr_pi8(a, b, 8);
  //     break;
  //   case 2:
  //     ret = _mm_alignr_pi8(a, b, 16);
  //     break;
  //   }
  //
  //   return VALIDATE_UINT8_M64(ret, d);
  // #endif
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_hadd_epi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  //   const int16_t *_b = (const int16_t *)impl.test_cases_int_pointer2;
  //   int16_t d[8];
  //   d[0] = _a[0] + _a[1];
  //   d[1] = _a[2] + _a[3];
  //   d[2] = _a[4] + _a[5];
  //   d[3] = _a[6] + _a[7];
  //   d[4] = _b[0] + _b[1];
  //   d[5] = _b[2] + _b[3];
  //   d[6] = _b[4] + _b[5];
  //   d[7] = _b[6] + _b[7];
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //   __m128i ret = _mm_hadd_epi16(a, b);
  //   return VALIDATE_INT16_M128(ret, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_hadd_epi32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int32_t *_a = (const int32_t *)impl.test_cases_int_pointer1;
  //   const int32_t *_b = (const int32_t *)impl.test_cases_int_pointer2;
  //   int32_t d[4];
  //   d[0] = _a[0] + _a[1];
  //   d[1] = _a[2] + _a[3];
  //   d[2] = _b[0] + _b[1];
  //   d[3] = _b[2] + _b[3];
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //   __m128i ret = _mm_hadd_epi32(a, b);
  //   return VALIDATE_INT32_M128(ret, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_hadd_pi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  //   const int16_t *_b = (const int16_t *)impl.test_cases_int_pointer2;
  //   int16_t d[4];
  //   d[0] = _a[0] + _a[1];
  //   d[1] = _a[2] + _a[3];
  //   d[2] = _b[0] + _b[1];
  //   d[3] = _b[2] + _b[3];
  //   __m64 a = load_m64(_a);
  //   __m64 b = load_m64(_b);
  //   __m64 ret = _mm_hadd_pi16(a, b);
  //   return VALIDATE_INT16_M64(ret, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_hadd_pi32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int32_t *_a = (const int32_t *)impl.test_cases_int_pointer1;
  //   const int32_t *_b = (const int32_t *)impl.test_cases_int_pointer2;
  //   int32_t d[2];
  //   d[0] = _a[0] + _a[1];
  //   d[1] = _b[0] + _b[1];
  //   __m64 a = load_m64(_a);
  //   __m64 b = load_m64(_b);
  //   __m64 ret = _mm_hadd_pi32(a, b);
  //   return VALIDATE_INT32_M64(ret, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_hadds_epi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  //   const int16_t *_b = (const int16_t *)impl.test_cases_int_pointer1;
  //
  //   int16_t d16[8];
  //   int32_t d32[8];
  //   d32[0] = (int32_t)_a[0] + (int32_t)_a[1];
  //   d32[1] = (int32_t)_a[2] + (int32_t)_a[3];
  //   d32[2] = (int32_t)_a[4] + (int32_t)_a[5];
  //   d32[3] = (int32_t)_a[6] + (int32_t)_a[7];
  //   d32[4] = (int32_t)_b[0] + (int32_t)_b[1];
  //   d32[5] = (int32_t)_b[2] + (int32_t)_b[3];
  //   d32[6] = (int32_t)_b[4] + (int32_t)_b[5];
  //   d32[7] = (int32_t)_b[6] + (int32_t)_b[7];
  //   for (int i = 0; i < 8; i++) {
  //     if (d32[i] > (int32_t)INT16_MAX)
  //       d16[i] = INT16_MAX;
  //     else if (d32[i] < (int32_t)INT16_MIN)
  //       d16[i] = INT16_MIN;
  //     else
  //       d16[i] = (int16_t)d32[i];
  //   }
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //   __m128i c = _mm_hadds_epi16(a, b);
  //
  //   return VALIDATE_INT16_M128(c, d16);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_hadds_pi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  //   const int16_t *_b = (const int16_t *)impl.test_cases_int_pointer1;
  //
  //   int16_t d16[8];
  //   int32_t d32[8];
  //   d32[0] = (int32_t)_a[0] + (int32_t)_a[1];
  //   d32[1] = (int32_t)_a[2] + (int32_t)_a[3];
  //   d32[2] = (int32_t)_b[0] + (int32_t)_b[1];
  //   d32[3] = (int32_t)_b[2] + (int32_t)_b[3];
  //   for (int i = 0; i < 8; i++) {
  //     if (d32[i] > (int32_t)INT16_MAX)
  //       d16[i] = INT16_MAX;
  //     else if (d32[i] < (int32_t)INT16_MIN)
  //       d16[i] = INT16_MIN;
  //     else
  //       d16[i] = (int16_t)d32[i];
  //   }
  //
  //   __m64 a = load_m64(_a);
  //   __m64 b = load_m64(_b);
  //   __m64 c = _mm_hadds_pi16(a, b);
  //
  //   return VALIDATE_INT16_M64(c, d16);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_hsub_epi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  //   const int16_t *_b = (const int16_t *)impl.test_cases_int_pointer1;
  //
  //   int16_t d[8];
  //   d[0] = _a[0] - _a[1];
  //   d[1] = _a[2] - _a[3];
  //   d[2] = _a[4] - _a[5];
  //   d[3] = _a[6] - _a[7];
  //   d[4] = _b[0] - _b[1];
  //   d[5] = _b[2] - _b[3];
  //   d[6] = _b[4] - _b[5];
  //   d[7] = _b[6] - _b[7];
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //   __m128i c = _mm_hsub_epi16(a, b);
  //
  //   return VALIDATE_INT16_M128(c, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_hsub_epi32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int32_t *_a = impl.test_cases_int_pointer1;
  //   const int32_t *_b = impl.test_cases_int_pointer1;
  //
  //   int32_t d[4];
  //   d[0] = _a[0] - _a[1];
  //   d[1] = _a[2] - _a[3];
  //   d[2] = _b[0] - _b[1];
  //   d[3] = _b[2] - _b[3];
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //   __m128i c = _mm_hsub_epi32(a, b);
  //
  //   return VALIDATE_INT32_M128(c, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_hsub_pi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  //   const int16_t *_b = (const int16_t *)impl.test_cases_int_pointer2;
  //
  //   int16_t d[4];
  //   d[0] = _a[0] - _a[1];
  //   d[1] = _a[2] - _a[3];
  //   d[2] = _b[0] - _b[1];
  //   d[3] = _b[2] - _b[3];
  //   __m64 a = load_m64(_a);
  //   __m64 b = load_m64(_b);
  //   __m64 c = _mm_hsub_pi16(a, b);
  //
  //   return VALIDATE_INT16_M64(c, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_hsub_pi32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int32_t *_a = impl.test_cases_int_pointer1;
  //   const int32_t *_b = impl.test_cases_int_pointer2;
  //
  //   int32_t d[2];
  //   d[0] = _a[0] - _a[1];
  //   d[1] = _b[0] - _b[1];
  //
  //   __m64 a = load_m64(_a);
  //   __m64 b = load_m64(_b);
  //   __m64 c = _mm_hsub_pi32(a, b);
  //
  //   return VALIDATE_INT32_M64(c, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_hsubs_epi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  //   const int16_t *_b = (const int16_t *)impl.test_cases_int_pointer1;
  //
  //   int16_t d16[8];
  //   int32_t d32[8];
  //   d32[0] = (int32_t)_a[0] - (int32_t)_a[1];
  //   d32[1] = (int32_t)_a[2] - (int32_t)_a[3];
  //   d32[2] = (int32_t)_a[4] - (int32_t)_a[5];
  //   d32[3] = (int32_t)_a[6] - (int32_t)_a[7];
  //   d32[4] = (int32_t)_b[0] - (int32_t)_b[1];
  //   d32[5] = (int32_t)_b[2] - (int32_t)_b[3];
  //   d32[6] = (int32_t)_b[4] - (int32_t)_b[5];
  //   d32[7] = (int32_t)_b[6] - (int32_t)_b[7];
  //   for (int i = 0; i < 8; i++) {
  //     if (d32[i] > (int32_t)INT16_MAX)
  //       d16[i] = INT16_MAX;
  //     else if (d32[i] < (int32_t)INT16_MIN)
  //       d16[i] = INT16_MIN;
  //     else
  //       d16[i] = (int16_t)d32[i];
  //   }
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //   __m128i c = _mm_hsubs_epi16(a, b);
  //
  //   return VALIDATE_INT16_M128(c, d16);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_hsubs_pi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  //   const int16_t *_b = (const int16_t *)impl.test_cases_int_pointer1;
  //
  //   int32_t _d[4];
  //   _d[0] = (int32_t)_a[0] - (int32_t)_a[1];
  //   _d[1] = (int32_t)_a[2] - (int32_t)_a[3];
  //   _d[2] = (int32_t)_b[0] - (int32_t)_b[1];
  //   _d[3] = (int32_t)_b[2] - (int32_t)_b[3];
  //
  //   for (int i = 0; i < 4; i++) {
  //     if (_d[i] > (int32_t)INT16_MAX) {
  //       _d[i] = INT16_MAX;
  //     } else if (_d[i] < (int32_t)INT16_MIN) {
  //       _d[i] = INT16_MIN;
  //     }
  //   }
  //
  //   __m64 a = load_m64(_a);
  //   __m64 b = load_m64(_b);
  //   __m64 c = _mm_hsubs_pi16(a, b);
  //
  //   return VALIDATE_INT16_M64(c, _d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_maddubs_epi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const uint8_t *_a = (const uint8_t *)impl.test_cases_int_pointer1;
  //   const int8_t *_b = (const int8_t *)impl.test_cases_int_pointer2;
  //   int32_t d0 = (int32_t)(_a[0] * _b[0]);
  //   int32_t d1 = (int32_t)(_a[1] * _b[1]);
  //   int32_t d2 = (int32_t)(_a[2] * _b[2]);
  //   int32_t d3 = (int32_t)(_a[3] * _b[3]);
  //   int32_t d4 = (int32_t)(_a[4] * _b[4]);
  //   int32_t d5 = (int32_t)(_a[5] * _b[5]);
  //   int32_t d6 = (int32_t)(_a[6] * _b[6]);
  //   int32_t d7 = (int32_t)(_a[7] * _b[7]);
  //   int32_t d8 = (int32_t)(_a[8] * _b[8]);
  //   int32_t d9 = (int32_t)(_a[9] * _b[9]);
  //   int32_t d10 = (int32_t)(_a[10] * _b[10]);
  //   int32_t d11 = (int32_t)(_a[11] * _b[11]);
  //   int32_t d12 = (int32_t)(_a[12] * _b[12]);
  //   int32_t d13 = (int32_t)(_a[13] * _b[13]);
  //   int32_t d14 = (int32_t)(_a[14] * _b[14]);
  //   int32_t d15 = (int32_t)(_a[15] * _b[15]);
  //
  //   int16_t e[8];
  //   e[0] = saturate_16(d0 + d1);
  //   e[1] = saturate_16(d2 + d3);
  //   e[2] = saturate_16(d4 + d5);
  //   e[3] = saturate_16(d6 + d7);
  //   e[4] = saturate_16(d8 + d9);
  //   e[5] = saturate_16(d10 + d11);
  //   e[6] = saturate_16(d12 + d13);
  //   e[7] = saturate_16(d14 + d15);
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //   __m128i c = _mm_maddubs_epi16(a, b);
  //   return VALIDATE_INT16_M128(c, e);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_maddubs_pi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const uint8_t *_a = (const uint8_t *)impl.test_cases_int_pointer1;
  //   const int8_t *_b = (const int8_t *)impl.test_cases_int_pointer2;
  //   int16_t d0 = (int16_t)(_a[0] * _b[0]);
  //   int16_t d1 = (int16_t)(_a[1] * _b[1]);
  //   int16_t d2 = (int16_t)(_a[2] * _b[2]);
  //   int16_t d3 = (int16_t)(_a[3] * _b[3]);
  //   int16_t d4 = (int16_t)(_a[4] * _b[4]);
  //   int16_t d5 = (int16_t)(_a[5] * _b[5]);
  //   int16_t d6 = (int16_t)(_a[6] * _b[6]);
  //   int16_t d7 = (int16_t)(_a[7] * _b[7]);
  //
  //   int16_t e[4];
  //   e[0] = saturate_16(d0 + d1);
  //   e[1] = saturate_16(d2 + d3);
  //   e[2] = saturate_16(d4 + d5);
  //   e[3] = saturate_16(d6 + d7);
  //
  //   __m64 a = load_m64(_a);
  //   __m64 b = load_m64(_b);
  //   __m64 c = _mm_maddubs_pi16(a, b);
  //
  //   return VALIDATE_INT16_M64(c, e);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_mulhrs_epi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  //   const int16_t *_b = (const int16_t *)impl.test_cases_int_pointer2;
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //   int32_t _c[8];
  //   for (int i = 0; i < 8; i++) {
  //     _c[i] = (((((int32_t)_a[i] * (int32_t)_b[i]) >> 14) + 1) & 0x1FFFE) >>
  //     1;
  //   }
  //   __m128i c = _mm_mulhrs_epi16(a, b);
  //
  //   return VALIDATE_INT16_M128(c, _c);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_mulhrs_pi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  //   const int16_t *_b = (const int16_t *)impl.test_cases_int_pointer2;
  //
  //   __m64 a = load_m64(_a);
  //   __m64 b = load_m64(_b);
  //   int32_t _c[4];
  //   for (int i = 0; i < 4; i++) {
  //     _c[i] = (((((int32_t)_a[i] * (int32_t)_b[i]) >> 14) + 1) & 0x1FFFE) >>
  //     1;
  //   }
  //   __m64 c = _mm_mulhrs_pi16(a, b);
  //
  //   return VALIDATE_INT16_M64(c, _c);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_shuffle_epi8(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int8_t *_a = (const int8_t *)impl.test_cases_int_pointer1;
  //   const int8_t *_b = (const int8_t *)impl.test_cases_int_pointer2;
  //   int8_t dst[16];
  //
  //   for (int i = 0; i < 16; i++) {
  //     if (_b[i] & 0x80) {
  //       dst[i] = 0;
  //     } else {
  //       dst[i] = _a[_b[i] & 0x0F];
  //     }
  //   }
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //   __m128i ret = _mm_shuffle_epi8(a, b);
  //
  //   return VALIDATE_INT8_M128(ret, dst);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_shuffle_pi8(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int8_t *_a = (const int8_t *)impl.test_cases_int_pointer1;
  //   const int8_t *_b = (const int8_t *)impl.test_cases_int_pointer2;
  //   int8_t dst[8];
  //
  //   for (int i = 0; i < 8; i++) {
  //     if (_b[i] & 0x80) {
  //       dst[i] = 0;
  //     } else {
  //       dst[i] = _a[_b[i] & 0x07];
  //     }
  //   }
  //
  //   __m64 a = load_m64(_a);
  //   __m64 b = load_m64(_b);
  //   __m64 ret = _mm_shuffle_pi8(a, b);
  //
  //   return VALIDATE_INT8_M64(ret, dst);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_sign_epi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  const int16_t *_b = (const int16_t *)impl.test_cases_int_pointer2;

  int16_t d[8];
  for (int i = 0; i < 8; i++) {
    if (_b[i] < 0) {
      d[i] = -_a[i];
    } else if (_b[i] == 0) {
      d[i] = 0;
    } else {
      d[i] = _a[i];
    }
  }

  __m128i a = load_m128i(_a);
  __m128i b = load_m128i(_b);
  __m128i c = _mm_sign_epi16(a, b);

  return VALIDATE_INT16_M128(c, d);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_sign_epi32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int32_t *_a = (const int32_t *)impl.test_cases_int_pointer1;
  const int32_t *_b = (const int32_t *)impl.test_cases_int_pointer2;

  int32_t d[4];
  for (int i = 0; i < 4; i++) {
    if (_b[i] < 0) {
      d[i] = -_a[i];
    } else if (_b[i] == 0) {
      d[i] = 0;
    } else {
      d[i] = _a[i];
    }
  }

  __m128i a = load_m128i(_a);
  __m128i b = load_m128i(_b);
  __m128i c = _mm_sign_epi32(a, b);

  return VALIDATE_INT32_M128(c, d);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_sign_epi8(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int8_t *_a = (const int8_t *)impl.test_cases_int_pointer1;
  const int8_t *_b = (const int8_t *)impl.test_cases_int_pointer2;

  int8_t d[16];
  for (int i = 0; i < 16; i++) {
    if (_b[i] < 0) {
      d[i] = -_a[i];
    } else if (_b[i] == 0) {
      d[i] = 0;
    } else {
      d[i] = _a[i];
    }
  }

  __m128i a = load_m128i(_a);
  __m128i b = load_m128i(_b);
  __m128i c = _mm_sign_epi8(a, b);

  return VALIDATE_INT8_M128(c, d);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_sign_pi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  const int16_t *_b = (const int16_t *)impl.test_cases_int_pointer2;

  int16_t d[4];
  for (int i = 0; i < 4; i++) {
    if (_b[i] < 0) {
      d[i] = -_a[i];
    } else if (_b[i] == 0) {
      d[i] = 0;
    } else {
      d[i] = _a[i];
    }
  }

  __m64 a = load_m64(_a);
  __m64 b = load_m64(_b);
  __m64 c = _mm_sign_pi16(a, b);

  return VALIDATE_INT16_M64(c, d);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_sign_pi32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int32_t *_a = (const int32_t *)impl.test_cases_int_pointer1;
  const int32_t *_b = (const int32_t *)impl.test_cases_int_pointer2;

  int32_t d[2];
  for (int i = 0; i < 2; i++) {
    if (_b[i] < 0) {
      d[i] = -_a[i];
    } else if (_b[i] == 0) {
      d[i] = 0;
    } else {
      d[i] = _a[i];
    }
  }

  __m64 a = load_m64(_a);
  __m64 b = load_m64(_b);
  __m64 c = _mm_sign_pi32(a, b);

  return VALIDATE_INT32_M64(c, d);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_sign_pi8(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int8_t *_a = (const int8_t *)impl.test_cases_int_pointer1;
  const int8_t *_b = (const int8_t *)impl.test_cases_int_pointer2;

  int8_t d[8];
  for (int i = 0; i < 8; i++) {
    if (_b[i] < 0) {
      d[i] = -_a[i];
    } else if (_b[i] == 0) {
      d[i] = 0;
    } else {
      d[i] = _a[i];
    }
  }

  __m64 a = load_m64(_a);
  __m64 b = load_m64(_b);
  __m64 c = _mm_sign_pi8(a, b);

  return VALIDATE_INT8_M64(c, d);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

/* SSE4.1 */
result_t test_mm_blend_epi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  const int16_t *_b = (const int16_t *)impl.test_cases_int_pointer2;
  int16_t _c[8];
  __m128i a, b, c;

#define TEST_IMPL(IDX)                                                         \
  for (int j = 0; j < 8; j++) {                                                \
    if ((IDX >> j) & 0x1) {                                                    \
      _c[j] = _b[j];                                                           \
    } else {                                                                   \
      _c[j] = _a[j];                                                           \
    }                                                                          \
  }                                                                            \
  a = load_m128i(_a);                                                          \
  b = load_m128i(_b);                                                          \
  c = _mm_blend_epi16(a, b, IDX);                                              \
  CHECK_RESULT(VALIDATE_INT16_M128(c, _c));

  IMM_256_ITER
#undef TEST_IMPL
  return TEST_SUCCESS;
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_blend_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const double *_a = (const double *)impl.test_cases_float_pointer1;
  const double *_b = (const double *)impl.test_cases_float_pointer2;
  __m128d a, b, c;

#define TEST_IMPL(IDX)                                                         \
  double _c##IDX[2];                                                           \
  for (int j = 0; j < 2; j++) {                                                \
    if ((IDX >> j) & 0x1) {                                                    \
      _c##IDX[j] = _b[j];                                                      \
    } else {                                                                   \
      _c##IDX[j] = _a[j];                                                      \
    }                                                                          \
  }                                                                            \
                                                                               \
  a = load_m128d(_a);                                                          \
  b = load_m128d(_b);                                                          \
  c = _mm_blend_pd(a, b, IDX);                                                 \
  CHECK_RESULT(validate_double(c, _c##IDX[0], _c##IDX[1]))

  IMM_4_ITER
#undef TEST_IMPL
  return TEST_SUCCESS;
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_blend_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const float *_a = impl.test_cases_float_pointer1;
  const float *_b = impl.test_cases_float_pointer2;
  __m128 a = load_m128(_a);
  __m128 b = load_m128(_b);
  __m128 c;

// gcc and clang can't compile call to _mm_blend_ps with 3rd argument as
// integer type due 4 bit size limitation.
#define TEST_IMPL(IDX)                                                         \
  float _c##IDX[4];                                                            \
  for (int i = 0; i < 4; i++) {                                                \
    if (IDX & (1 << i)) {                                                      \
      _c##IDX[i] = _b[i];                                                      \
    } else {                                                                   \
      _c##IDX[i] = _a[i];                                                      \
    }                                                                          \
  }                                                                            \
                                                                               \
  c = _mm_blend_ps(a, b, IDX);                                                 \
  CHECK_RESULT(                                                                \
      validate_float(c, _c##IDX[0], _c##IDX[1], _c##IDX[2], _c##IDX[3]))

  IMM_4_ITER
#undef TEST_IMPL
  return TEST_SUCCESS;
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_blendv_epi8(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int8_t *_a = (const int8_t *)impl.test_cases_int_pointer1;
  const int8_t *_b = (const int8_t *)impl.test_cases_int_pointer2;
  const int8_t _mask[16] = {(const int8_t)impl.test_cases_ints[iter],
                            (const int8_t)impl.test_cases_ints[iter + 1],
                            (const int8_t)impl.test_cases_ints[iter + 2],
                            (const int8_t)impl.test_cases_ints[iter + 3],
                            (const int8_t)impl.test_cases_ints[iter + 4],
                            (const int8_t)impl.test_cases_ints[iter + 5],
                            (const int8_t)impl.test_cases_ints[iter + 6],
                            (const int8_t)impl.test_cases_ints[iter + 7]};

  int8_t _c[16];
  for (int i = 0; i < 16; i++) {
    if (_mask[i] >> 7) {
      _c[i] = _b[i];
    } else {
      _c[i] = _a[i];
    }
  }

  __m128i a = load_m128i(_a);
  __m128i b = load_m128i(_b);
  __m128i mask = load_m128i(_mask);
  __m128i c = _mm_blendv_epi8(a, b, mask);

  return VALIDATE_INT8_M128(c, _c);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_blendv_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const double *_a = (const double *)impl.test_cases_float_pointer1;
  const double *_b = (const double *)impl.test_cases_float_pointer2;
  const double _mask[] = {(double)impl.test_cases_floats[iter],
                          (double)impl.test_cases_floats[iter + 1]};

  double _c[2];
  for (int i = 0; i < 2; i++) {
    // signed shift right would return a result which is either all 1's
    // from
    // negative numbers or all 0's from positive numbers
    if ((*(const int64_t *)(_mask + i)) >> 63) {
      _c[i] = _b[i];
    } else {
      _c[i] = _a[i];
    }
  }

  __m128d a = load_m128d(_a);
  __m128d b = load_m128d(_b);
  __m128d mask = load_m128d(_mask);

  __m128d c = _mm_blendv_pd(a, b, mask);

  return validate_double(c, _c[0], _c[1]);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_blendv_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const float *_a = impl.test_cases_float_pointer1;
  const float *_b = impl.test_cases_float_pointer2;
  const float _mask[] = {
      impl.test_cases_floats[iter], impl.test_cases_floats[iter + 1],
      impl.test_cases_floats[iter + 2], impl.test_cases_floats[iter + 3]};

  float _c[4];
  for (int i = 0; i < 4; i++) {
    // signed shift right would return a result which is either all 1's
    // from
    // negative numbers or all 0's from positive numbers
    if ((*(const int32_t *)(_mask + i)) >> 31) {
      _c[i] = _b[i];
    } else {
      _c[i] = _a[i];
    }
  }

  __m128 a = load_m128(_a);
  __m128 b = load_m128(_b);
  __m128 mask = load_m128(_mask);

  __m128 c = _mm_blendv_ps(a, b, mask);

  return validate_float(c, _c[0], _c[1], _c[2], _c[3]);
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_ceil_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const double *_a = (const double *)impl.test_cases_float_pointer1;
  //
  //   double dx = ceil(_a[0]);
  //   double dy = ceil(_a[1]);
  //
  //   __m128d a = load_m128d(_a);
  //   __m128d ret = _mm_ceil_pd(a);
  //
  //   return validate_double(ret, dx, dy);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_ceil_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   float dx = ceilf(_a[0]);
  //   float dy = ceilf(_a[1]);
  //   float dz = ceilf(_a[2]);
  //   float dw = ceilf(_a[3]);
  //
  //   __m128 a = _mm_load_ps(_a);
  //   __m128 c = _mm_ceil_ps(a);
  //   return validate_float(c, dx, dy, dz, dw);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_ceil_sd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const double *_a = (const double *)impl.test_cases_float_pointer1;
  //   const double *_b = (const double *)impl.test_cases_float_pointer2;
  //
  //   double dx = ceil(_b[0]);
  //   double dy = _a[1];
  //
  //   __m128d a = load_m128d(_a);
  //   __m128d b = load_m128d(_b);
  //   __m128d ret = _mm_ceil_sd(a, b);
  //
  //   return validate_double(ret, dx, dy);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_ceil_ss(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   const float *_b = impl.test_cases_float_pointer1;
  //
  //   float f0 = ceilf(_b[0]);
  //
  //   __m128 a = load_m128(_a);
  //   __m128 b = load_m128(_b);
  //   __m128 c = _mm_ceil_ss(a, b);
  //
  //   return validate_float(c, f0, _a[1], _a[2], _a[3]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cmpeq_epi64(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int64_t *_a = (const int64_t *)impl.test_cases_int_pointer1;
  //   const int64_t *_b = (const int64_t *)impl.test_cases_int_pointer2;
  //   int64_t d0 = (_a[0] == _b[0]) ? 0xffffffffffffffff : 0x0;
  //   int64_t d1 = (_a[1] == _b[1]) ? 0xffffffffffffffff : 0x0;
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //   __m128i c = _mm_cmpeq_epi64(a, b);
  //   return validate_int64(c, d0, d1);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtepi16_epi32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  //
  //   int32_t d[4];
  //   d[0] = (int32_t)_a[0];
  //   d[1] = (int32_t)_a[1];
  //   d[2] = (int32_t)_a[2];
  //   d[3] = (int32_t)_a[3];
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i ret = _mm_cvtepi16_epi32(a);
  //
  //   return VALIDATE_INT32_M128(ret, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtepi16_epi64(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  //
  //   int64_t i0 = (int64_t)_a[0];
  //   int64_t i1 = (int64_t)_a[1];
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i ret = _mm_cvtepi16_epi64(a);
  //
  //   return validate_int64(ret, i0, i1);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtepi32_epi64(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int32_t *_a = (const int32_t *)impl.test_cases_int_pointer1;
  //
  //   int64_t i0 = (int64_t)_a[0];
  //   int64_t i1 = (int64_t)_a[1];
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i ret = _mm_cvtepi32_epi64(a);
  //
  //   return validate_int64(ret, i0, i1);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtepi8_epi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int8_t *_a = (const int8_t *)impl.test_cases_int_pointer1;
  //
  //   int16_t d[8];
  //   d[0] = (int16_t)_a[0];
  //   d[1] = (int16_t)_a[1];
  //   d[2] = (int16_t)_a[2];
  //   d[3] = (int16_t)_a[3];
  //   d[4] = (int16_t)_a[4];
  //   d[5] = (int16_t)_a[5];
  //   d[6] = (int16_t)_a[6];
  //   d[7] = (int16_t)_a[7];
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i ret = _mm_cvtepi8_epi16(a);
  //
  //   return VALIDATE_INT16_M128(ret, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtepi8_epi32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int8_t *_a = (const int8_t *)impl.test_cases_int_pointer1;
  //
  //   int32_t d[4];
  //   d[0] = (int32_t)_a[0];
  //   d[1] = (int32_t)_a[1];
  //   d[2] = (int32_t)_a[2];
  //   d[3] = (int32_t)_a[3];
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i ret = _mm_cvtepi8_epi32(a);
  //
  //   return VALIDATE_INT32_M128(ret, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtepi8_epi64(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int8_t *_a = (const int8_t *)impl.test_cases_int_pointer1;
  //
  //   int64_t i0 = (int64_t)_a[0];
  //   int64_t i1 = (int64_t)_a[1];
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i ret = _mm_cvtepi8_epi64(a);
  //
  //   return validate_int64(ret, i0, i1);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtepu16_epi32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const uint16_t *_a = (const uint16_t *)impl.test_cases_int_pointer1;
  //
  //   int32_t d[4];
  //   d[0] = (int32_t)_a[0];
  //   d[1] = (int32_t)_a[1];
  //   d[2] = (int32_t)_a[2];
  //   d[3] = (int32_t)_a[3];
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i ret = _mm_cvtepu16_epi32(a);
  //
  //   return VALIDATE_INT32_M128(ret, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtepu16_epi64(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const uint16_t *_a = (const uint16_t *)impl.test_cases_int_pointer1;
  //
  //   int64_t i0 = (int64_t)_a[0];
  //   int64_t i1 = (int64_t)_a[1];
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i ret = _mm_cvtepu16_epi64(a);
  //
  //   return validate_int64(ret, i0, i1);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtepu32_epi64(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const uint32_t *_a = (const uint32_t *)impl.test_cases_int_pointer1;
  //
  //   int64_t i0 = (int64_t)_a[0];
  //   int64_t i1 = (int64_t)_a[1];
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i ret = _mm_cvtepu32_epi64(a);
  //
  //   return validate_int64(ret, i0, i1);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtepu8_epi16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const uint8_t *_a = (const uint8_t *)impl.test_cases_int_pointer1;
  //
  //   int16_t d[8];
  //   d[0] = (int16_t)_a[0];
  //   d[1] = (int16_t)_a[1];
  //   d[2] = (int16_t)_a[2];
  //   d[3] = (int16_t)_a[3];
  //   d[4] = (int16_t)_a[4];
  //   d[5] = (int16_t)_a[5];
  //   d[6] = (int16_t)_a[6];
  //   d[7] = (int16_t)_a[7];
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i ret = _mm_cvtepu8_epi16(a);
  //
  //   return VALIDATE_INT16_M128(ret, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtepu8_epi32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const uint8_t *_a = (const uint8_t *)impl.test_cases_int_pointer1;
  //
  //   int32_t d[4];
  //   d[0] = (int32_t)_a[0];
  //   d[1] = (int32_t)_a[1];
  //   d[2] = (int32_t)_a[2];
  //   d[3] = (int32_t)_a[3];
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i ret = _mm_cvtepu8_epi32(a);
  //
  //   return VALIDATE_INT32_M128(ret, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cvtepu8_epi64(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const uint8_t *_a = (const uint8_t *)impl.test_cases_int_pointer1;
  //
  //   int64_t i0 = (int64_t)_a[0];
  //   int64_t i1 = (int64_t)_a[1];
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i ret = _mm_cvtepu8_epi64(a);
  //
  //   return validate_int64(ret, i0, i1);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

#define MM_DP_PD_TEST_CASE_WITH(imm8)                                          \
  do {                                                                         \
    const double *_a = (const double *)impl.test_cases_float_pointer1;         \
    const double *_b = (const double *)impl.test_cases_float_pointer2;         \
    const int imm = imm8;                                                      \
    double d[2];                                                               \
    double sum = 0;                                                            \
    for (size_t i = 0; i < 2; i++)                                             \
      sum += ((imm) & (1 << (i + 4))) ? _a[i] * _b[i] : 0;                     \
    for (size_t i = 0; i < 2; i++)                                             \
      d[i] = (imm & (1 << i)) ? sum : 0;                                       \
    __m128d a = load_m128d(_a);                                                \
    __m128d b = load_m128d(_b);                                                \
    __m128d ret = _mm_dp_pd(a, b, imm);                                        \
    if (validate_double(ret, d[0], d[1]) != TEST_SUCCESS)                      \
      return TEST_FAIL;                                                        \
  } while (0)

#define GENERATE_MM_DP_PD_TEST_CASES                                           \
  MM_DP_PD_TEST_CASE_WITH(0xF0);                                               \
  MM_DP_PD_TEST_CASE_WITH(0xF1);                                               \
  MM_DP_PD_TEST_CASE_WITH(0xF2);                                               \
  MM_DP_PD_TEST_CASE_WITH(0xFF);                                               \
  MM_DP_PD_TEST_CASE_WITH(0x10);                                               \
  MM_DP_PD_TEST_CASE_WITH(0x11);                                               \
  MM_DP_PD_TEST_CASE_WITH(0x12);                                               \
  MM_DP_PD_TEST_CASE_WITH(0x13);                                               \
  MM_DP_PD_TEST_CASE_WITH(0x00);                                               \
  MM_DP_PD_TEST_CASE_WITH(0x01);                                               \
  MM_DP_PD_TEST_CASE_WITH(0x02);                                               \
  MM_DP_PD_TEST_CASE_WITH(0x03);                                               \
  MM_DP_PD_TEST_CASE_WITH(0x20);                                               \
  MM_DP_PD_TEST_CASE_WITH(0x21);                                               \
  MM_DP_PD_TEST_CASE_WITH(0x22);                                               \
  MM_DP_PD_TEST_CASE_WITH(0x23);

result_t test_mm_dp_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   GENERATE_MM_DP_PD_TEST_CASES
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

#define MM_DP_PS_TEST_CASE_WITH(IMM)                                           \
  do {                                                                         \
    const float *_a = impl.test_cases_float_pointer1;                          \
    const float *_b = impl.test_cases_float_pointer2;                          \
    const int imm = IMM;                                                       \
    __m128 a = load_m128(_a);                                                  \
    __m128 b = load_m128(_b);                                                  \
    __m128 out = _mm_dp_ps(a, b, imm);                                         \
    float r[4]; /* the reference */                                            \
    float sum = 0;                                                             \
    for (size_t i = 0; i < 4; i++)                                             \
      sum += ((imm) & (1 << (i + 4))) ? _a[i] * _b[i] : 0;                     \
    for (size_t i = 0; i < 4; i++)                                             \
      r[i] = (imm & (1 << i)) ? sum : 0;                                       \
    /* the epsilon has to be large enough, otherwise test suite fails. */      \
    if (validate_float_epsilon(out, r[0], r[1], r[2], r[3], 2050.0f) !=        \
        TEST_SUCCESS)                                                          \
      return TEST_FAIL;                                                        \
  } while (0)

#define GENERATE_MM_DP_PS_TEST_CASES                                           \
  MM_DP_PS_TEST_CASE_WITH(0xFF);                                               \
  MM_DP_PS_TEST_CASE_WITH(0x7F);                                               \
  MM_DP_PS_TEST_CASE_WITH(0x9F);                                               \
  MM_DP_PS_TEST_CASE_WITH(0x2F);                                               \
  MM_DP_PS_TEST_CASE_WITH(0x0F);                                               \
  MM_DP_PS_TEST_CASE_WITH(0x23);                                               \
  MM_DP_PS_TEST_CASE_WITH(0xB5);

result_t test_mm_dp_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   GENERATE_MM_DP_PS_TEST_CASES
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_extract_epi32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   int32_t *_a = (int32_t *)impl.test_cases_int_pointer1;
  //   __m128i a = load_m128i(_a);
  //   int c;
  //
  // #define TEST_IMPL(IDX)
  //   c = _mm_extract_epi32(a, IDX);
  // ASSERT_RETURN(c == *(_a + IDX));
  //
  //   IMM_4_ITER
  // #undef TEST_IMPL
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_extract_epi64(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   int64_t *_a = (int64_t *)impl.test_cases_int_pointer1;
  //   __m128i a = load_m128i(_a);
  //   __int64 c;
  //
  // #define TEST_IMPL(IDX)
  //   c = _mm_extract_epi64(a, IDX);
  //  ASSERT_RETURN(c == *(_a + IDX));
  //
  //   IMM_2_ITER
  // #undef TEST_IMPL
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_extract_epi8(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   uint8_t *_a = (uint8_t *)impl.test_cases_int_pointer1;
  //   __m128i a = load_m128i(_a);
  //   int c;
  //
  // #define TEST_IMPL(IDX)
  //   c = _mm_extract_epi8(a, IDX);
  //  ASSERT_RETURN(c == *(_a + IDX));
  //
  //   IMM_8_ITER
  // #undef TEST_IMPL
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_extract_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = (const float *)impl.test_cases_float_pointer1;
  //
  //   __m128 a = _mm_load_ps(_a);
  //   int32_t c;
  //
  // #define TEST_IMPL(IDX)
  //   c = _mm_extract_ps(a, IDX);
  //  ASSERT_RETURN(c == *(const int32_t *)(_a + IDX));
  //
  //   IMM_4_ITER
  // #undef TEST_IMPL
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_floor_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const double *_a = (const double *)impl.test_cases_float_pointer1;
  //
  //   double dx = floor(_a[0]);
  //   double dy = floor(_a[1]);
  //
  //   __m128d a = load_m128d(_a);
  //   __m128d ret = _mm_floor_pd(a);
  //
  //   return validate_double(ret, dx, dy);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_floor_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   float dx = floorf(_a[0]);
  //   float dy = floorf(_a[1]);
  //   float dz = floorf(_a[2]);
  //   float dw = floorf(_a[3]);
  //
  //   __m128 a = load_m128(_a);
  //   __m128 c = _mm_floor_ps(a);
  //   return validate_float(c, dx, dy, dz, dw);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_floor_sd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const double *_a = (const double *)impl.test_cases_float_pointer1;
  //   const double *_b = (const double *)impl.test_cases_float_pointer2;
  //
  //   double dx = floor(_b[0]);
  //   double dy = _a[1];
  //
  //   __m128d a = load_m128d(_a);
  //   __m128d b = load_m128d(_b);
  //   __m128d ret = _mm_floor_sd(a, b);
  //
  //   return validate_double(ret, dx, dy);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_floor_ss(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   const float *_b = impl.test_cases_float_pointer1;
  //
  //   float f0 = floorf(_b[0]);
  //
  //   __m128 a = load_m128(_a);
  //   __m128 b = load_m128(_b);
  //   __m128 c = _mm_floor_ss(a, b);
  //
  //   return validate_float(c, f0, _a[1], _a[2], _a[3]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_insert_epi32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int32_t *_a = (const int32_t *)impl.test_cases_int_pointer1;
  //   const int32_t insert = (int32_t)*impl.test_cases_int_pointer2;
  //   __m128i a, b;
  //
  // #define TEST_IMPL(IDX)
  //   int32_t d##IDX[4];
  //   for (int i = 0; i < 4; i++) {
  //     d##IDX[i] = _a[i];
  //   }
  //   d##IDX[IDX] = insert;
  //
  //   a = load_m128i(_a);
  //   b = _mm_insert_epi32(a, (int)insert, IDX);
  //   CHECK_RESULT(VALIDATE_INT32_M128(b, d##IDX));
  //
  //   IMM_4_ITER
  // #undef TEST_IMPL
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_insert_epi64(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int64_t *_a = (const int64_t *)impl.test_cases_int_pointer1;
  //   int64_t insert = (int64_t)*impl.test_cases_int_pointer2;
  //
  //   __m128i a, b;
  //   int64_t d[2];
  // #define TEST_IMPL(IDX)
  //   d[0] = _a[0];
  //   d[1] = _a[1];
  //   d[IDX] = insert;
  //   a = load_m128i(_a);
  //   b = _mm_insert_epi64(a, insert, IDX);
  // CHECK_RESULT(validate_int64(b, d[0], d[1]));
  //
  //   IMM_2_ITER
  // #undef TEST_IMPL
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_insert_epi8(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int8_t *_a = (const int8_t *)impl.test_cases_int_pointer1;
  //   const int8_t insert = (int8_t)*impl.test_cases_int_pointer2;
  //   __m128i a, b;
  //   int8_t d[16];
  //
  // #define TEST_IMPL(IDX)
  //   for (int i = 0; i < 16; i++) {
  //     d[i] = _a[i];
  //   }
  //   d[IDX] = insert;
  //   a = load_m128i(_a);
  //   b = _mm_insert_epi8(a, insert, IDX);
  // CHECK_RESULT(VALIDATE_INT8_M128(b, d));
  //
  //   IMM_16_ITER
  // #undef TEST_IMPL
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_insert_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   const float *_b = impl.test_cases_float_pointer2;
  //
  //   __m128 a, b, c;
  // #define TEST_IMPL(IDX)
  //   float d##IDX[4] = {_a[0], _a[1], _a[2], _a[3]};
  //   d##IDX[(IDX >> 4) & 0x3] = _b[(IDX >> 6) & 0x3];
  //
  //   for (int j = 0; j < 4; j++) {
  //     if (IDX & (1 << j)) {
  //       d##IDX[j] = 0;
  //     }
  //   }
  //
  //   a = _mm_load_ps(_a);
  //   b = _mm_load_ps(_b);
  //   c = _mm_insert_ps(a, b, IDX);
  // CHECK_RESULT(validate_float(c, d##IDX[0], d##IDX[1], d##IDX[2],
  // d##IDX[3]));
  //
  //   IMM_256_ITER
  // #undef TEST_IMPL
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_max_epi32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int32_t *_a = (const int32_t *)impl.test_cases_int_pointer1;
  //   const int32_t *_b = (const int32_t *)impl.test_cases_int_pointer2;
  //
  //   int32_t d[4];
  //   d[0] = _a[0] > _b[0] ? _a[0] : _b[0];
  //   d[1] = _a[1] > _b[1] ? _a[1] : _b[1];
  //   d[2] = _a[2] > _b[2] ? _a[2] : _b[2];
  //   d[3] = _a[3] > _b[3] ? _a[3] : _b[3];
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //   __m128i c = _mm_max_epi32(a, b);
  //
  //   return VALIDATE_INT32_M128(c, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_max_epi8(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int8_t *_a = (const int8_t *)impl.test_cases_int_pointer1;
  //   const int8_t *_b = (const int8_t *)impl.test_cases_int_pointer2;
  //   int8_t d[16];
  //   d[0] = _a[0] > _b[0] ? _a[0] : _b[0];
  //   d[1] = _a[1] > _b[1] ? _a[1] : _b[1];
  //   d[2] = _a[2] > _b[2] ? _a[2] : _b[2];
  //   d[3] = _a[3] > _b[3] ? _a[3] : _b[3];
  //   d[4] = _a[4] > _b[4] ? _a[4] : _b[4];
  //   d[5] = _a[5] > _b[5] ? _a[5] : _b[5];
  //   d[6] = _a[6] > _b[6] ? _a[6] : _b[6];
  //   d[7] = _a[7] > _b[7] ? _a[7] : _b[7];
  //   d[8] = _a[8] > _b[8] ? _a[8] : _b[8];
  //   d[9] = _a[9] > _b[9] ? _a[9] : _b[9];
  //   d[10] = _a[10] > _b[10] ? _a[10] : _b[10];
  //   d[11] = _a[11] > _b[11] ? _a[11] : _b[11];
  //   d[12] = _a[12] > _b[12] ? _a[12] : _b[12];
  //   d[13] = _a[13] > _b[13] ? _a[13] : _b[13];
  //   d[14] = _a[14] > _b[14] ? _a[14] : _b[14];
  //   d[15] = _a[15] > _b[15] ? _a[15] : _b[15];
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //
  //   __m128i c = _mm_max_epi8(a, b);
  //   return VALIDATE_INT8_M128(c, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_max_epu16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const uint16_t *_a = (const uint16_t *)impl.test_cases_int_pointer1;
  //   const uint16_t *_b = (const uint16_t *)impl.test_cases_int_pointer2;
  //
  //   uint16_t d[8];
  //   d[0] = _a[0] > _b[0] ? _a[0] : _b[0];
  //   d[1] = _a[1] > _b[1] ? _a[1] : _b[1];
  //   d[2] = _a[2] > _b[2] ? _a[2] : _b[2];
  //   d[3] = _a[3] > _b[3] ? _a[3] : _b[3];
  //   d[4] = _a[4] > _b[4] ? _a[4] : _b[4];
  //   d[5] = _a[5] > _b[5] ? _a[5] : _b[5];
  //   d[6] = _a[6] > _b[6] ? _a[6] : _b[6];
  //   d[7] = _a[7] > _b[7] ? _a[7] : _b[7];
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //   __m128i c = _mm_max_epu16(a, b);
  //
  //   return VALIDATE_UINT16_M128(c, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_max_epu32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const uint32_t *_a = (const uint32_t *)impl.test_cases_int_pointer1;
  //   const uint32_t *_b = (const uint32_t *)impl.test_cases_int_pointer2;
  //
  //   uint32_t d[4];
  //   d[0] = _a[0] > _b[0] ? _a[0] : _b[0];
  //   d[1] = _a[1] > _b[1] ? _a[1] : _b[1];
  //   d[2] = _a[2] > _b[2] ? _a[2] : _b[2];
  //   d[3] = _a[3] > _b[3] ? _a[3] : _b[3];
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //   __m128i c = _mm_max_epu32(a, b);
  //
  //   return VALIDATE_UINT32_M128(c, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_min_epi32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int32_t *_a = (const int32_t *)impl.test_cases_int_pointer1;
  //   const int32_t *_b = (const int32_t *)impl.test_cases_int_pointer2;
  //
  //   int32_t d[4];
  //   d[0] = _a[0] < _b[0] ? _a[0] : _b[0];
  //   d[1] = _a[1] < _b[1] ? _a[1] : _b[1];
  //   d[2] = _a[2] < _b[2] ? _a[2] : _b[2];
  //   d[3] = _a[3] < _b[3] ? _a[3] : _b[3];
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //   __m128i c = _mm_min_epi32(a, b);
  //
  //   return VALIDATE_INT32_M128(c, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_min_epi8(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int8_t *_a = (const int8_t *)impl.test_cases_int_pointer1;
  //   const int8_t *_b = (const int8_t *)impl.test_cases_int_pointer2;
  //
  //   int8_t d[16];
  //   d[0] = _a[0] < _b[0] ? _a[0] : _b[0];
  //   d[1] = _a[1] < _b[1] ? _a[1] : _b[1];
  //   d[2] = _a[2] < _b[2] ? _a[2] : _b[2];
  //   d[3] = _a[3] < _b[3] ? _a[3] : _b[3];
  //   d[4] = _a[4] < _b[4] ? _a[4] : _b[4];
  //   d[5] = _a[5] < _b[5] ? _a[5] : _b[5];
  //   d[6] = _a[6] < _b[6] ? _a[6] : _b[6];
  //   d[7] = _a[7] < _b[7] ? _a[7] : _b[7];
  //   d[8] = _a[8] < _b[8] ? _a[8] : _b[8];
  //   d[9] = _a[9] < _b[9] ? _a[9] : _b[9];
  //   d[10] = _a[10] < _b[10] ? _a[10] : _b[10];
  //   d[11] = _a[11] < _b[11] ? _a[11] : _b[11];
  //   d[12] = _a[12] < _b[12] ? _a[12] : _b[12];
  //   d[13] = _a[13] < _b[13] ? _a[13] : _b[13];
  //   d[14] = _a[14] < _b[14] ? _a[14] : _b[14];
  //   d[15] = _a[15] < _b[15] ? _a[15] : _b[15];
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //
  //   __m128i c = _mm_min_epi8(a, b);
  //   return VALIDATE_INT8_M128(c, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_min_epu16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const uint16_t *_a = (const uint16_t *)impl.test_cases_int_pointer1;
  //   const uint16_t *_b = (const uint16_t *)impl.test_cases_int_pointer2;
  //
  //   uint16_t d[8];
  //   d[0] = _a[0] < _b[0] ? _a[0] : _b[0];
  //   d[1] = _a[1] < _b[1] ? _a[1] : _b[1];
  //   d[2] = _a[2] < _b[2] ? _a[2] : _b[2];
  //   d[3] = _a[3] < _b[3] ? _a[3] : _b[3];
  //   d[4] = _a[4] < _b[4] ? _a[4] : _b[4];
  //   d[5] = _a[5] < _b[5] ? _a[5] : _b[5];
  //   d[6] = _a[6] < _b[6] ? _a[6] : _b[6];
  //   d[7] = _a[7] < _b[7] ? _a[7] : _b[7];
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //   __m128i c = _mm_min_epu16(a, b);
  //
  //   return VALIDATE_UINT16_M128(c, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_min_epu32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const uint32_t *_a = (const uint32_t *)impl.test_cases_int_pointer1;
  //   const uint32_t *_b = (const uint32_t *)impl.test_cases_int_pointer2;
  //
  //   uint32_t d[4];
  //   d[0] = _a[0] < _b[0] ? _a[0] : _b[0];
  //   d[1] = _a[1] < _b[1] ? _a[1] : _b[1];
  //   d[2] = _a[2] < _b[2] ? _a[2] : _b[2];
  //   d[3] = _a[3] < _b[3] ? _a[3] : _b[3];
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //   __m128i c = _mm_min_epu32(a, b);
  //
  //   return VALIDATE_UINT32_M128(c, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_minpos_epu16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int16_t *_a = (const int16_t *)impl.test_cases_int_pointer1;
  //   uint16_t index = 0, min = (uint16_t)_a[0];
  //   for (int i = 0; i < 8; i++) {
  //     if ((uint16_t)_a[i] < min) {
  //       index = (uint16_t)i;
  //       min = (uint16_t)_a[i];
  //     }
  //   }
  //
  //   uint16_t d[8] = {min, index, 0, 0, 0, 0, 0, 0};
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i ret = _mm_minpos_epu16(a);
  //   return VALIDATE_UINT16_M128(ret, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_mpsadbw_epu8(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const uint8_t *_a = (const uint8_t *)impl.test_cases_int_pointer1;
  //   const uint8_t *_b = (const uint8_t *)impl.test_cases_int_pointer2;
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //   __m128i c;
  // #define TEST_IMPL(IDX)
  //   uint8_t a_offset##IDX = ((IDX >> 2) & 0x1) * 4;
  //   uint8_t b_offset##IDX = (IDX & 0x3) * 4;
  //
  //   uint16_t d##IDX[8] = {};
  //   for (int i = 0; i < 8; i++) {
  //     for (int j = 0; j < 4; j++) {
  //       d##IDX[i] += abs(_a[(a_offset##IDX + i) + j] - _b[b_offset##IDX +
  //       j]);
  //     }
  //   }
  //   c = _mm_mpsadbw_epu8(a, b, IDX);
  // CHECK_RESULT(VALIDATE_UINT16_M128(c,d##IDX));
  //
  //   IMM_8_ITER
  // #undef TEST_IMPL
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_mul_epi32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int32_t *_a = (const int32_t *)impl.test_cases_int_pointer1;
  //   const int32_t *_b = (const int32_t *)impl.test_cases_int_pointer2;
  //
  //   int64_t dx = (int64_t)(_a[0]) * (int64_t)(_b[0]);
  //   int64_t dy = (int64_t)(_a[2]) * (int64_t)(_b[2]);
  //
  //   __m128i a = _mm_loadu_si128((const __m128i *)_a);
  //   __m128i b = _mm_loadu_si128((const __m128i *)_b);
  //   __m128i r = _mm_mul_epi32(a, b);
  //
  //   return validate_int64(r, dx, dy);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_mullo_epi32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int32_t *_a = impl.test_cases_int_pointer1;
  //   const int32_t *_b = impl.test_cases_int_pointer2;
  //   int32_t d[4];
  //
  //   for (int i = 0; i < 4; i++) {
  //     d[i] = (int32_t)((int64_t)_a[i] * (int64_t)_b[i]);
  //   }
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //   __m128i c = _mm_mullo_epi32(a, b);
  //   return VALIDATE_INT32_M128(c, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_packus_epi32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   uint16_t max = UINT16_MAX;
  //   uint16_t min = 0;
  //   const int32_t *_a = (const int32_t *)impl.test_cases_int_pointer1;
  //   const int32_t *_b = (const int32_t *)impl.test_cases_int_pointer2;
  //
  //   uint16_t d[8];
  //   for (int i = 0; i < 4; i++) {
  //     if (_a[i] > (int32_t)max)
  //       d[i] = max;
  //     else if (_a[i] < (int32_t)min)
  //       d[i] = min;
  //     else
  //       d[i] = (uint16_t)_a[i];
  //   }
  //   for (int i = 0; i < 4; i++) {
  //     if (_b[i] > (int32_t)max)
  //       d[i + 4] = max;
  //     else if (_b[i] < (int32_t)min)
  //       d[i + 4] = min;
  //     else
  //       d[i + 4] = (uint16_t)_b[i];
  //   }
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //   __m128i c = _mm_packus_epi32(a, b);
  //
  //   return VALIDATE_UINT16_M128(c, d);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_round_pd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const double *_a = (double *)impl.test_cases_float_pointer1;
  //   double d[2];
  //   __m128d ret;
  //
  //   __m128d a = load_m128d(_a);
  //   switch (iter & 0x7) {
  //   case 0:
  //     d[0] = bankers_rounding(_a[0]);
  //     d[1] = bankers_rounding(_a[1]);
  //
  //     ret = _mm_round_pd(a, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
  //     break;
  //   case 1:
  //     d[0] = floor(_a[0]);
  //     d[1] = floor(_a[1]);
  //
  //     ret = _mm_round_pd(a, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
  //     break;
  //   case 2:
  //     d[0] = ceil(_a[0]);
  //     d[1] = ceil(_a[1]);
  //
  //     ret = _mm_round_pd(a, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
  //     break;
  //   case 3:
  //     d[0] = _a[0] > 0 ? floor(_a[0]) : ceil(_a[0]);
  //     d[1] = _a[1] > 0 ? floor(_a[1]) : ceil(_a[1]);
  //
  //     ret = _mm_round_pd(a, _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC);
  //     break;
  //   case 4:
  //     d[0] = bankers_rounding(_a[0]);
  //     d[1] = bankers_rounding(_a[1]);
  //
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_NEAREST);
  //     ret = _mm_round_pd(a, _MM_FROUND_CUR_DIRECTION);
  //     break;
  //   case 5:
  //     d[0] = floor(_a[0]);
  //     d[1] = floor(_a[1]);
  //
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_DOWN);
  //     ret = _mm_round_pd(a, _MM_FROUND_CUR_DIRECTION);
  //     break;
  //   case 6:
  //     d[0] = ceil(_a[0]);
  //     d[1] = ceil(_a[1]);
  //
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
  //     ret = _mm_round_pd(a, _MM_FROUND_CUR_DIRECTION);
  //     break;
  //   case 7:
  //     d[0] = _a[0] > 0 ? floor(_a[0]) : ceil(_a[0]);
  //     d[1] = _a[1] > 0 ? floor(_a[1]) : ceil(_a[1]);
  //
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_TOWARD_ZERO);
  //     ret = _mm_round_pd(a, _MM_FROUND_CUR_DIRECTION);
  //     break;
  //   }
  //
  //   return validate_double(ret, d[0], d[1]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_round_ps(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   float f[4];
  //   __m128 ret;
  //
  //   __m128 a = load_m128(_a);
  //   switch (iter & 0x7) {
  //   case 0:
  //     f[0] = bankers_rounding(_a[0]);
  //     f[1] = bankers_rounding(_a[1]);
  //     f[2] = bankers_rounding(_a[2]);
  //     f[3] = bankers_rounding(_a[3]);
  //
  //     ret = _mm_round_ps(a, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
  //     break;
  //   case 1:
  //     f[0] = floorf(_a[0]);
  //     f[1] = floorf(_a[1]);
  //     f[2] = floorf(_a[2]);
  //     f[3] = floorf(_a[3]);
  //
  //     ret = _mm_round_ps(a, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
  //     break;
  //   case 2:
  //     f[0] = ceilf(_a[0]);
  //     f[1] = ceilf(_a[1]);
  //     f[2] = ceilf(_a[2]);
  //     f[3] = ceilf(_a[3]);
  //
  //     ret = _mm_round_ps(a, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
  //     break;
  //   case 3:
  //     f[0] = _a[0] > 0 ? floorf(_a[0]) : ceilf(_a[0]);
  //     f[1] = _a[1] > 0 ? floorf(_a[1]) : ceilf(_a[1]);
  //     f[2] = _a[2] > 0 ? floorf(_a[2]) : ceilf(_a[2]);
  //     f[3] = _a[3] > 0 ? floorf(_a[3]) : ceilf(_a[3]);
  //
  //     ret = _mm_round_ps(a, _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC);
  //     break;
  //   case 4:
  //     f[0] = bankers_rounding(_a[0]);
  //     f[1] = bankers_rounding(_a[1]);
  //     f[2] = bankers_rounding(_a[2]);
  //     f[3] = bankers_rounding(_a[3]);
  //
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_NEAREST);
  //     ret = _mm_round_ps(a, _MM_FROUND_CUR_DIRECTION);
  //     break;
  //   case 5:
  //     f[0] = floorf(_a[0]);
  //     f[1] = floorf(_a[1]);
  //     f[2] = floorf(_a[2]);
  //     f[3] = floorf(_a[3]);
  //
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_DOWN);
  //     ret = _mm_round_ps(a, _MM_FROUND_CUR_DIRECTION);
  //     break;
  //   case 6:
  //     f[0] = ceilf(_a[0]);
  //     f[1] = ceilf(_a[1]);
  //     f[2] = ceilf(_a[2]);
  //     f[3] = ceilf(_a[3]);
  //
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
  //     ret = _mm_round_ps(a, _MM_FROUND_CUR_DIRECTION);
  //     break;
  //   case 7:
  //     f[0] = _a[0] > 0 ? floorf(_a[0]) : ceilf(_a[0]);
  //     f[1] = _a[1] > 0 ? floorf(_a[1]) : ceilf(_a[1]);
  //     f[2] = _a[2] > 0 ? floorf(_a[2]) : ceilf(_a[2]);
  //     f[3] = _a[3] > 0 ? floorf(_a[3]) : ceilf(_a[3]);
  //
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_TOWARD_ZERO);
  //     ret = _mm_round_ps(a, _MM_FROUND_CUR_DIRECTION);
  //     break;
  //   }
  //
  //   return validate_float(ret, f[0], f[1], f[2], f[3]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_round_sd(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const double *_a = (double *)impl.test_cases_float_pointer1;
  //   const double *_b = (double *)impl.test_cases_float_pointer2;
  //   double d[2];
  //   __m128d ret;
  //
  //   __m128d a = load_m128d(_a);
  //   __m128d b = load_m128d(_b);
  //   d[1] = _a[1];
  //   switch (iter & 0x7) {
  //   case 0:
  //     d[0] = bankers_rounding(_b[0]);
  //
  //     ret = _mm_round_sd(a, b, _MM_FROUND_TO_NEAREST_INT |
  //     _MM_FROUND_NO_EXC); break;
  //   case 1:
  //     d[0] = floor(_b[0]);
  //
  //     ret = _mm_round_sd(a, b, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
  //     break;
  //   case 2:
  //     d[0] = ceil(_b[0]);
  //
  //     ret = _mm_round_sd(a, b, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
  //     break;
  //   case 3:
  //     d[0] = _b[0] > 0 ? floor(_b[0]) : ceil(_b[0]);
  //
  //     ret = _mm_round_sd(a, b, _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC);
  //     break;
  //   case 4:
  //     d[0] = bankers_rounding(_b[0]);
  //
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_NEAREST);
  //     ret = _mm_round_sd(a, b, _MM_FROUND_CUR_DIRECTION);
  //     break;
  //   case 5:
  //     d[0] = floor(_b[0]);
  //
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_DOWN);
  //     ret = _mm_round_sd(a, b, _MM_FROUND_CUR_DIRECTION);
  //     break;
  //   case 6:
  //     d[0] = ceil(_b[0]);
  //
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
  //     ret = _mm_round_sd(a, b, _MM_FROUND_CUR_DIRECTION);
  //     break;
  //   case 7:
  //     d[0] = _b[0] > 0 ? floor(_b[0]) : ceil(_b[0]);
  //
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_TOWARD_ZERO);
  //     ret = _mm_round_sd(a, b, _MM_FROUND_CUR_DIRECTION);
  //     break;
  //   }
  //
  //   return validate_double(ret, d[0], d[1]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_round_ss(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const float *_a = impl.test_cases_float_pointer1;
  //   const float *_b = impl.test_cases_float_pointer2;
  //   float f[4];
  //   __m128 ret;
  //
  //   __m128 a = load_m128(_a);
  //   __m128 b = load_m128(_b);
  //   switch (iter & 0x7) {
  //   case 0:
  //     f[0] = bankers_rounding(_b[0]);
  //
  //     ret = _mm_round_ss(a, b, _MM_FROUND_TO_NEAREST_INT |
  //     _MM_FROUND_NO_EXC); break;
  //   case 1:
  //     f[0] = floorf(_b[0]);
  //
  //     ret = _mm_round_ss(a, b, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
  //     break;
  //   case 2:
  //     f[0] = ceilf(_b[0]);
  //
  //     ret = _mm_round_ss(a, b, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
  //     break;
  //   case 3:
  //     f[0] = _b[0] > 0 ? floorf(_b[0]) : ceilf(_b[0]);
  //
  //     ret = _mm_round_ss(a, b, _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC);
  //     break;
  //   case 4:
  //     f[0] = bankers_rounding(_b[0]);
  //
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_NEAREST);
  //     ret = _mm_round_ss(a, b, _MM_FROUND_CUR_DIRECTION);
  //     break;
  //   case 5:
  //     f[0] = floorf(_b[0]);
  //
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_DOWN);
  //     ret = _mm_round_ss(a, b, _MM_FROUND_CUR_DIRECTION);
  //     break;
  //   case 6:
  //     f[0] = ceilf(_b[0]);
  //
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
  //     ret = _mm_round_ss(a, b, _MM_FROUND_CUR_DIRECTION);
  //     break;
  //   case 7:
  //     f[0] = _b[0] > 0 ? floorf(_b[0]) : ceilf(_b[0]);
  //
  //     _MM_SET_ROUNDING_MODE(_MM_ROUND_TOWARD_ZERO);
  //     ret = _mm_round_ss(a, b, _MM_FROUND_CUR_DIRECTION);
  //     break;
  //   }
  //   f[1] = _a[1];
  //   f[2] = _a[2];
  //   f[3] = _a[3];
  //
  //   return validate_float(ret, f[0], f[1], f[2], f[3]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_stream_load_si128(const SSE2RVV_TEST_IMPL &impl,
                                   // #ifdef ENABLE_TEST_ALL
                                   uint32_t iter) {
  //   int32_t *addr = impl.test_cases_int_pointer1;
  //
  //   __m128i ret = _mm_stream_load_si128((__m128i *)addr);
  //
  //   return VALIDATE_INT32_M128(ret, addr);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_test_all_ones(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int32_t *_a = (const int32_t *)impl.test_cases_int_pointer1;
  __m128i a = load_m128i(_a);

  int32_t d0 = ~_a[0] & (~(uint32_t)0);
  int32_t d1 = ~_a[1] & (~(uint32_t)0);
  int32_t d2 = ~_a[2] & (~(uint32_t)0);
  int32_t d3 = ~_a[3] & (~(uint32_t)0);
  int32_t result = ((d0 | d1 | d2 | d3) == 0) ? 1 : 0;

  int32_t ret = _mm_test_all_ones(a);

  return result == ret ? TEST_SUCCESS : TEST_FAIL;
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_test_all_zeros(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int32_t *_a = (const int32_t *)impl.test_cases_int_pointer1;
  const int32_t *_mask = (const int32_t *)impl.test_cases_int_pointer2;
  __m128i a = load_m128i(_a);
  __m128i mask = load_m128i(_mask);

  int32_t _c[4];
  _c[0] = _a[0] & _mask[0];
  _c[1] = _a[1] & _mask[1];
  _c[2] = _a[2] & _mask[2];
  _c[3] = _a[3] & _mask[3];
  int32_t c = ((_c[0] | _c[1] | _c[2] | _c[3]) == 0) ? 1 : 0;
  int32_t ret = _mm_test_all_zeros(a, mask);
  return c == ret ? TEST_SUCCESS : TEST_FAIL;
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_test_mix_ones_zeros(const SSE2RVV_TEST_IMPL &impl,
                                     uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int32_t *_a = (const int32_t *)impl.test_cases_int_pointer1;
  const int32_t *_mask = (const int32_t *)impl.test_cases_int_pointer2;
  __m128i a = load_m128i(_a);
  __m128i mask = load_m128i(_mask);

  int32_t d[4];
  d[0] = !((_a[0]) & _mask[0]) & !((~_a[0]) & _mask[0]);
  d[1] = !((_a[1]) & _mask[1]) & !((~_a[1]) & _mask[1]);
  d[2] = !((_a[2]) & _mask[2]) & !((~_a[2]) & _mask[2]);
  d[3] = !((_a[3]) & _mask[3]) & !((~_a[3]) & _mask[3]);
  int32_t result = (!d[0] & !d[1] & !d[2] & !d[3]) ? 1 : 0;
  int32_t ret = _mm_test_mix_ones_zeros(a, mask);

  return result == ret ? TEST_SUCCESS : TEST_FAIL;
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_testc_si128(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int32_t *_a = impl.test_cases_int_pointer1;
  const int32_t *_b = impl.test_cases_int_pointer2;
  __m128i a = load_m128i(_a);
  __m128i b = load_m128i(_b);
  int testc = 1;
  for (int i = 0; i < 2; i++) {
    if ((~(((SIMDVec *)&a)->m128_u64[i]) & ((SIMDVec *)&b)->m128_u64[i])) {
      testc = 0;
      break;
    }
  }
  return _mm_testc_si128(a, b) == testc ? TEST_SUCCESS : TEST_FAIL;
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_testnzc_si128(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int32_t *_a = (const int32_t *)impl.test_cases_int_pointer1;
  const int32_t *_b = (const int32_t *)impl.test_cases_int_pointer2;
  __m128i a = load_m128i(_a);
  __m128i b = load_m128i(_b);

  int32_t d[4];
  d[0] = !((_a[0]) & _b[0]) & !((~_a[0]) & _b[0]);
  d[1] = !((_a[1]) & _b[1]) & !((~_a[1]) & _b[1]);
  d[2] = !((_a[2]) & _b[2]) & !((~_a[2]) & _b[2]);
  d[3] = !((_a[3]) & _b[3]) & !((~_a[3]) & _b[3]);
  int32_t result = (!d[0] & !d[1] & !d[2] & !d[3]) ? 1 : 0;
  int32_t ret = _mm_testnzc_si128(a, b);

  return result == ret ? TEST_SUCCESS : TEST_FAIL;
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

result_t test_mm_testz_si128(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
#ifdef ENABLE_TEST_ALL
  const int32_t *_a = impl.test_cases_int_pointer1;
  const int32_t *_b = impl.test_cases_int_pointer2;
  __m128i a = load_m128i(_a);
  __m128i b = load_m128i(_b);
  int testz = 1;
  for (int i = 0; i < 2; i++) {
    if ((((SIMDVec *)&a)->m128_u64[i] & ((SIMDVec *)&b)->m128_u64[i])) {
      testz = 0;
      break;
    }
  }
  return _mm_testz_si128(a, b) == testz ? TEST_SUCCESS : TEST_FAIL;
#else
  return TEST_UNIMPL;
#endif // ENABLE_TEST_ALL
}

/* SSE4.2 */

result_t test_mm_cmpestrc(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   GENERATE_MM_CMPESTRC_TEST_CASES
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cmpgt_epi64(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int64_t *_a = (const int64_t *)impl.test_cases_int_pointer1;
  //   const int64_t *_b = (const int64_t *)impl.test_cases_int_pointer2;
  //
  //   int64_t result[2];
  //   result[0] = _a[0] > _b[0] ? -1 : 0;
  //   result[1] = _a[1] > _b[1] ? -1 : 0;
  //
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //   __m128i iret = _mm_cmpgt_epi64(a, b);
  //
  //   return validate_int64(iret, result[0], result[1]);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cmpistrs(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   GENERATE_MM_CMPISTRS_TEST_CASES
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_cmpistrz(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   GENERATE_MM_CMPISTRZ_TEST_CASES
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_crc32_u16(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   uint32_t crc = *(const uint32_t *)impl.test_cases_int_pointer1;
  //   uint16_t v = iter;
  //   uint32_t result = _mm_crc32_u16(crc, v);
  //   ASSERT_RETURN(result == canonical_crc32_u16(crc, v));
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_crc32_u32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   uint32_t crc = *(const uint32_t *)impl.test_cases_int_pointer1;
  //   uint32_t v = *(const uint32_t *)impl.test_cases_int_pointer2;
  //   uint32_t result = _mm_crc32_u32(crc, v);
  //   ASSERT_RETURN(result == canonical_crc32_u32(crc, v));
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_crc32_u64(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   uint64_t crc = *(const uint64_t *)impl.test_cases_int_pointer1;
  //   uint64_t v = *(const uint64_t *)impl.test_cases_int_pointer2;
  //   uint64_t result = _mm_crc32_u64(crc, v);
  //   ASSERT_RETURN(result == canonical_crc32_u64(crc, v));
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_crc32_u8(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   uint32_t crc = *(const uint32_t *)impl.test_cases_int_pointer1;
  //   uint8_t v = iter;
  //   uint32_t result = _mm_crc32_u8(crc, v);
  //   ASSERT_RETURN(result == canonical_crc32_u8(crc, v));
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

/* AES */
result_t test_mm_aesenc_si128(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int32_t *a = (int32_t *)impl.test_cases_int_pointer1;
  //   const int32_t *b = (int32_t *)impl.test_cases_int_pointer2;
  //   __m128i data = _mm_loadu_si128((const __m128i *)a);
  //   __m128i rk = _mm_loadu_si128((const __m128i *)b);
  //
  //   __m128i resultReference = aesenc_128_reference(data, rk);
  //   __m128i resultIntrinsic = _mm_aesenc_si128(data, rk);
  //
  //   return validate_128bits(resultReference, resultIntrinsic);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_aesdec_si128(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const int32_t *a = (int32_t *)impl.test_cases_int_pointer1;
  //   const int32_t *b = (int32_t *)impl.test_cases_int_pointer2;
  //   __m128i data = _mm_loadu_si128((const __m128i *)a);
  //   __m128i rk = _mm_loadu_si128((const __m128i *)b);
  //
  //   __m128i resultReference = aesdec_128_reference(data, rk);
  //   __m128i resultIntrinsic = _mm_aesdec_si128(data, rk);
  //
  //   return validate_128bits(resultReference, resultIntrinsic);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_aesenclast_si128(const SSE2RVV_TEST_IMPL &impl,
                                  // #ifdef ENABLE_TEST_ALL
                                  uint32_t iter) {
  //   const int32_t *a = (const int32_t *)impl.test_cases_int_pointer1;
  //   const int32_t *b = (const int32_t *)impl.test_cases_int_pointer2;
  //   __m128i data = _mm_loadu_si128((const __m128i *)a);
  //   __m128i rk = _mm_loadu_si128((const __m128i *)b);
  //
  //   __m128i resultReference = aesenclast_128_reference(data, rk);
  //   __m128i resultIntrinsic = _mm_aesenclast_si128(data, rk);
  //
  //   return validate_128bits(resultReference, resultIntrinsic);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_aesdeclast_si128(const SSE2RVV_TEST_IMPL &impl,
                                  // #ifdef ENABLE_TEST_ALL
                                  uint32_t iter) {
  //   const uint8_t *a = (uint8_t *)impl.test_cases_int_pointer1;
  //   const uint8_t *rk = (uint8_t *)impl.test_cases_int_pointer2;
  //   __m128i _a = _mm_loadu_si128((const __m128i *)a);
  //   __m128i _rk = _mm_loadu_si128((const __m128i *)rk);
  //   uint8_t c[16] = {};
  //
  //   uint8_t v[4][4];
  //   for (int i = 0; i < 16; ++i) {
  //     v[((i / 4) + (i % 4)) % 4][i % 4] = crypto_aes_rsbox[a[i]];
  //   }
  //   for (int i = 0; i < 16; ++i) {
  //     c[i] = v[i / 4][i % 4] ^ rk[i];
  //   }
  //
  //   __m128i result_reference = _mm_loadu_si128((const __m128i *)c);
  //   __m128i result_intrinsic = _mm_aesdeclast_si128(_a, _rk);
  //
  //   return validate_128bits(result_reference, result_intrinsic);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_aesimc_si128(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const uint8_t *a = (uint8_t *)impl.test_cases_int_pointer1;
  //   __m128i _a = _mm_loadu_si128((const __m128i *)a);
  //
  //   uint8_t e, f, g, h, v[4][4];
  //   for (int i = 0; i < 16; ++i) {
  //     ((uint8_t *)v)[i] = a[i];
  //   }
  //   for (int i = 0; i < 4; ++i) {
  //     e = v[i][0];
  //     f = v[i][1];
  //     g = v[i][2];
  //     h = v[i][3];
  //
  //     v[i][0] = MULTIPLY(e, 0x0e) ^ MULTIPLY(f, 0x0b) ^ MULTIPLY(g, 0x0d) ^
  //               MULTIPLY(h, 0x09);
  //     v[i][1] = MULTIPLY(e, 0x09) ^ MULTIPLY(f, 0x0e) ^ MULTIPLY(g, 0x0b) ^
  //               MULTIPLY(h, 0x0d);
  //     v[i][2] = MULTIPLY(e, 0x0d) ^ MULTIPLY(f, 0x09) ^ MULTIPLY(g, 0x0e) ^
  //               MULTIPLY(h, 0x0b);
  //     v[i][3] = MULTIPLY(e, 0x0b) ^ MULTIPLY(f, 0x0d) ^ MULTIPLY(g, 0x09) ^
  //               MULTIPLY(h, 0x0e);
  //   }
  //
  //   __m128i result_reference = _mm_loadu_si128((const __m128i *)v);
  //   __m128i result_intrinsic = _mm_aesimc_si128(_a);
  //
  //   return validate_128bits(result_reference, result_intrinsic);
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

static inline uint32_t sub_word(uint32_t in) {
  return (crypto_aes_sbox[(in >> 24) & 0xff] << 24) |
         (crypto_aes_sbox[(in >> 16) & 0xff] << 16) |
         (crypto_aes_sbox[(in >> 8) & 0xff] << 8) |
         (crypto_aes_sbox[in & 0xff]);
}

// FIXME: improve the test case for AES-256 key expansion.
// Reference:
// https://github.com/randombit/botan/blob/master/src/lib/block/aes/aes_ni/aes_ni.cpp
result_t test_mm_aeskeygenassist_si128(const SSE2RVV_TEST_IMPL &impl,
                                       // #ifdef ENABLE_TEST_ALL
                                       uint32_t iter) {
  //   const uint32_t *a = (uint32_t *)impl.test_cases_int_pointer1;
  //   __m128i data = load_m128i(a);
  //   uint32_t sub_x1 = sub_word(a[1]);
  //   uint32_t sub_x3 = sub_word(a[3]);
  //   __m128i result_reference;
  //   __m128i result_intrinsic;
  // #define TEST_IMPL(IDX)
  //   uint32_t res##IDX[4] = {
  //       sub_x1,
  //       rotr(sub_x1, 8) ^ IDX,
  //       sub_x3,
  //       rotr(sub_x3, 8) ^ IDX,
  //   };
  //   result_reference = load_m128i(res##IDX);
  //   result_intrinsic = _mm_aeskeygenassist_si128(data, IDX);
  //   CHECK_RESULT(validate_128bits(result_reference, result_intrinsic));
  //
  //   IMM_256_ITER
  // #undef TEST_IMPL
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

/* Others */
result_t test_mm_clmulepi64_si128(const SSE2RVV_TEST_IMPL &impl,
                                  // #ifdef ENABLE_TEST_ALL
                                  uint32_t iter) {
  //   const uint64_t *_a = (const uint64_t *)impl.test_cases_int_pointer1;
  //   const uint64_t *_b = (const uint64_t *)impl.test_cases_int_pointer2;
  //   __m128i a = load_m128i(_a);
  //   __m128i b = load_m128i(_b);
  //   auto result = clmul_64(_a[0], _b[0]);
  //   if (!validate_uint64(_mm_clmulepi64_si128(a, b, 0x00), result.first,
  //                        result.second))
  //     return TEST_FAIL;
  //   result = clmul_64(_a[1], _b[0]);
  //   if (!validate_uint64(_mm_clmulepi64_si128(a, b, 0x01), result.first,
  //                        result.second))
  //     return TEST_FAIL;
  //   result = clmul_64(_a[0], _b[1]);
  //   if (!validate_uint64(_mm_clmulepi64_si128(a, b, 0x10), result.first,
  //                        result.second))
  //     return TEST_FAIL;
  //   result = clmul_64(_a[1], _b[1]);
  //   if (!validate_uint64(_mm_clmulepi64_si128(a, b, 0x11), result.first,
  //                        result.second))
  //     return TEST_FAIL;
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_get_denormals_zero_mode(const SSE2RVV_TEST_IMPL &impl,
                                         // #ifdef ENABLE_TEST_ALL
                                         uint32_t iter) {
  //   int res_denormals_zero_on, res_denormals_zero_off;
  //
  //   _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
  //   res_denormals_zero_on =
  //       _MM_GET_DENORMALS_ZERO_MODE() == _MM_DENORMALS_ZERO_ON;
  //
  //   _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_OFF);
  //   res_denormals_zero_off =
  //       _MM_GET_DENORMALS_ZERO_MODE() == _MM_DENORMALS_ZERO_OFF;
  //
  //   return (res_denormals_zero_on && res_denormals_zero_off) ? TEST_SUCCESS
  //                                                            : TEST_FAIL;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

// static int popcnt_reference(uint64_t a) {
// int count = 0;
// while (a != 0) {
//   count += a & 1;
//   a >>= 1;
// }
// return count;
// }

result_t test_mm_popcnt_u32(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const uint64_t *a = (const uint64_t *)impl.test_cases_int_pointer1;
  //   ASSERT_RETURN(popcnt_reference((uint32_t)a[0]) == _mm_popcnt_u32(a[0]));
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_popcnt_u64(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   const uint64_t *a = (const uint64_t *)impl.test_cases_int_pointer1;
  //   ASSERT_RETURN(popcnt_reference(a[0]) == _mm_popcnt_u64(a[0]));
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_mm_set_denormals_zero_mode(const SSE2RVV_TEST_IMPL &impl,
                                         // #ifdef ENABLE_TEST_ALL
                                         uint32_t iter) {
  //   result_t res_set_denormals_zero_on, res_set_denormals_zero_off;
  //   float factor = 2;
  //   float denormal = FLT_MIN / factor;
  //   float denormals[4] = {denormal, denormal, denormal, denormal};
  //   float factors[4] = {factor, factor, factor, factor};
  //   __m128 ret;
  //
  //   _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
  //   ret = _mm_mul_ps(load_m128(denormals), load_m128(factors));
  //   res_set_denormals_zero_on = validate_float(ret, 0, 0, 0, 0);
  //
  //   _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_OFF);
  //   ret = _mm_mul_ps(load_m128(denormals), load_m128(factors));
  //   res_set_denormals_zero_off =
  //       validate_float(ret, FLT_MIN, FLT_MIN, FLT_MIN, FLT_MIN);
  //
  //   if (res_set_denormals_zero_on == TEST_FAIL ||
  //       res_set_denormals_zero_off == TEST_FAIL)
  //     return TEST_FAIL;
  //   return TEST_SUCCESS;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

result_t test_rdtsc(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  //   uint64_t start = _rdtsc();
  //   for (int i = 0; i < 100000; i++) {
  // #if defined(_MSC_VER)
  //     _ReadWriteBarrier();
  // #else
  //     __asm__ __volatile__("" ::: "memory");
  // #endif
  //   }
  //   uint64_t end = _rdtsc();
  //   return end > start ? TEST_SUCCESS : TEST_FAIL;
  // #else
  return TEST_UNIMPL;
  // #endif  // ENABLE_TEST_ALL
}

#if defined(__riscv_v_elen)
#define REGISTER_SIZE __riscv_v_elen
#elif defined(__aarch64__)
#define REGISTER_SIZE 128
#elif (defined(__x86_64__) || defined(__i386__))
#define REGISTER_SIZE sizeof(__m128)
#endif

SSE2RVV_TEST_IMPL::SSE2RVV_TEST_IMPL(void) {
  test_cases_float_pointer1 = (float *)platform_aligned_alloc(REGISTER_SIZE);
  test_cases_float_pointer2 = (float *)platform_aligned_alloc(REGISTER_SIZE);
  test_cases_int_pointer1 = (int32_t *)platform_aligned_alloc(REGISTER_SIZE);
  test_cases_int_pointer2 = (int32_t *)platform_aligned_alloc(REGISTER_SIZE);
  SSE2RVV_INIT_RNG(123456);
  for (uint32_t i = 0; i < MAX_TEST_VALUE; i++) {
    test_cases_floats[i] = ranf(-100000, 100000);
    test_cases_ints[i] = (int32_t)ranf(-100000, 100000);
  }
}

// Dummy function to match the case label in run_single_test.
result_t test_last(const SSE2RVV_TEST_IMPL &impl, uint32_t iter) {
  // #ifdef ENABLE_TEST_ALL
  return TEST_SUCCESS;
}

result_t SSE2RVV_TEST_IMPL::load_test_float_pointers(uint32_t i) {
  result_t ret = do_mm_store_ps(
      test_cases_float_pointer1, test_cases_floats[i], test_cases_floats[i + 1],
      test_cases_floats[i + 2], test_cases_floats[i + 3]);
  if (ret == TEST_SUCCESS) {
    ret = do_mm_store_ps(test_cases_float_pointer2, test_cases_floats[i + 4],
                         test_cases_floats[i + 5], test_cases_floats[i + 6],
                         test_cases_floats[i + 7]);
  }
  return ret;
}

result_t SSE2RVV_TEST_IMPL::load_test_int_pointers(uint32_t i) {
  result_t ret = do_mm_store_ps(test_cases_int_pointer1, test_cases_ints[i],
                                test_cases_ints[i + 1], test_cases_ints[i + 2],
                                test_cases_ints[i + 3]);
  if (ret == TEST_SUCCESS) {
    ret = do_mm_store_ps(test_cases_int_pointer2, test_cases_ints[i + 4],
                         test_cases_ints[i + 5], test_cases_ints[i + 6],
                         test_cases_ints[i + 7]);
  }

  return ret;
}

result_t SSE2RVV_TEST_IMPL::run_single_test(INSTRUCTION_TEST test, uint32_t i) {
  result_t ret = TEST_SUCCESS;

  switch (test) {
#define _(x)                                                                   \
  case it_##x:                                                                 \
    ret = test_##x(*this, i);                                                  \
    break;
    INTRIN_LIST
#undef _
  }

  return ret;
}

const char *instruction_string[] = {
#define _(x) #x,
    INTRIN_LIST
#undef _
};

SSE2RVV_TEST *SSE2RVV_TEST::create(void) {
  SSE2RVV_TEST_IMPL *st = new SSE2RVV_TEST_IMPL;
  return static_cast<SSE2RVV_TEST *>(st);
}

} // namespace SSE2RVV
