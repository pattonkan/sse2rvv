// compat_test.cpp
//
// Focused regression tests for the RVV intrinsics compatibility fixes
// (old-spec compilers: Bianbu GCC 13/14).  Each test covers one of the three
// API differences that were patched:
//
//  1. __riscv_vreinterpret_v_<int>_b* / _b*_<int>  → vlm / vsm
//  2. __RISCV_VXRM_RNU / vaaddu / vnclip without rounding-mode param
//  3. __RISCV_FRM_* / vfcvt_x_f_v_i32m1_rm without explicit frm param
//
// Build (cross): riscv64-linux-gnu-g++ -I. -march=rv64gcv_zba -static \
//                  tests/compat_test.cpp -o tests/compat_test -lm

#include "sse2rvv.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

static int failures = 0;

#define CHECK(cond, msg)                                                       \
  do {                                                                         \
    if (!(cond)) {                                                             \
      printf("FAIL: %s\n", msg);                                              \
      ++failures;                                                              \
    } else {                                                                   \
      printf("PASS: %s\n", msg);                                              \
    }                                                                          \
  } while (0)

// ── helpers ──────────────────────────────────────────────────────────────────

static void get_f32(__m128 v, float out[4]) { memcpy(out, &v, 16); }
static void get_f64(__m128d v, double out[2]) { memcpy(out, &v, 16); }
static void get_i8(__m128i v, int8_t out[16]) { memcpy(out, &v, 16); }
static void get_i16(__m128i v, int16_t out[8]) { memcpy(out, &v, 16); }
static void get_i32(__m128i v, int32_t out[4]) { memcpy(out, &v, 16); }
static void get_u16(__m128i v, uint16_t out[8]) { memcpy(out, &v, 16); }
static void get_u8(__m128i v, uint8_t out[16]) { memcpy(out, &v, 16); }

// ── Group 1: vlm / vsm mask load-store ──────────────────────────────────────

// _mm_addsub_ps: even lanes subtract, odd lanes add
static void test_addsub_ps(void) {
  float a[4] = {1.0f, 2.0f, 3.0f, 4.0f};
  float b[4] = {0.5f, 0.5f, 0.5f, 0.5f};
  __m128 va, vb, vc;
  memcpy(&va, a, 16);
  memcpy(&vb, b, 16);
  vc = _mm_addsub_ps(va, vb);
  float r[4];
  get_f32(vc, r);
  // expected: [1-0.5, 2+0.5, 3-0.5, 4+0.5] = [0.5, 2.5, 2.5, 4.5]
  CHECK(r[0] == 0.5f && r[1] == 2.5f && r[2] == 2.5f && r[3] == 4.5f,
        "_mm_addsub_ps");
}

// _mm_blend_epi16: imm8=0b00001100 selects b for lanes 2,3
static void test_blend_epi16(void) {
  int16_t a[8] = {1, 2, 3, 4, 5, 6, 7, 8};
  int16_t b[8] = {10, 20, 30, 40, 50, 60, 70, 80};
  __m128i va, vb, vc;
  memcpy(&va, a, 16);
  memcpy(&vb, b, 16);
  vc = _mm_blend_epi16(va, vb, 0x0C); // lanes 2,3 from b
  int16_t r[8];
  get_i16(vc, r);
  CHECK(r[0] == 1 && r[1] == 2 && r[2] == 30 && r[3] == 40 && r[4] == 5 &&
            r[5] == 6 && r[6] == 7 && r[7] == 8,
        "_mm_blend_epi16(imm=0x0C)");
}

// _mm_blend_pd: imm8=0x1 → lane 0 from b
static void test_blend_pd(void) {
  double a[2] = {1.0, 2.0};
  double b[2] = {10.0, 20.0};
  __m128d va, vb, vc;
  memcpy(&va, a, 16);
  memcpy(&vb, b, 16);
  vc = _mm_blend_pd(va, vb, 0x1);
  double r[2];
  get_f64(vc, r);
  CHECK(r[0] == 10.0 && r[1] == 2.0, "_mm_blend_pd(imm=0x1)");
}

// _mm_blend_ps: imm8=0x5 (0101) → lanes 0,2 from b
static void test_blend_ps(void) {
  float a[4] = {1.0f, 2.0f, 3.0f, 4.0f};
  float b[4] = {10.0f, 20.0f, 30.0f, 40.0f};
  __m128 va, vb, vc;
  memcpy(&va, a, 16);
  memcpy(&vb, b, 16);
  vc = _mm_blend_ps(va, vb, 0x5);
  float r[4];
  get_f32(vc, r);
  CHECK(r[0] == 10.0f && r[1] == 2.0f && r[2] == 30.0f && r[3] == 4.0f,
        "_mm_blend_ps(imm=0x5)");
}

// _mm_hadd_pd: r[0]=a[0]+a[1], r[1]=b[0]+b[1]
static void test_hadd_pd(void) {
  double a[2] = {1.0, 2.0};
  double b[2] = {3.0, 4.0};
  __m128d va, vb, vc;
  memcpy(&va, a, 16);
  memcpy(&vb, b, 16);
  vc = _mm_hadd_pd(va, vb);
  double r[2];
  get_f64(vc, r);
  CHECK(r[0] == 3.0 && r[1] == 7.0, "_mm_hadd_pd");
}

// _mm_hsub_pd: r[0]=a[0]-a[1], r[1]=b[0]-b[1]
static void test_hsub_pd(void) {
  double a[2] = {5.0, 2.0};
  double b[2] = {10.0, 3.0};
  __m128d va, vb, vc;
  memcpy(&va, a, 16);
  memcpy(&vb, b, 16);
  vc = _mm_hsub_pd(va, vb);
  double r[2];
  get_f64(vc, r);
  CHECK(r[0] == 3.0 && r[1] == 7.0, "_mm_hsub_pd");
}

// _mm_movemask_epi8: sign bits of 16 bytes
static void test_movemask_epi8(void) {
  int8_t a[16] = {-1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1};
  __m128i va;
  memcpy(&va, a, 16);
  int mask = _mm_movemask_epi8(va);
  // lanes 0,2,15 are negative → bits 0,2,15 set = 0x8005
  CHECK(mask == 0x8005, "_mm_movemask_epi8");
}

// _mm_movemask_pd: sign bits of 2 doubles
static void test_movemask_pd(void) {
  double a[2] = {-1.0, 1.0};
  __m128d va;
  memcpy(&va, a, 16);
  int mask = _mm_movemask_pd(va);
  CHECK(mask == 0x1, "_mm_movemask_pd");
}

// _mm_movemask_ps: sign bits of 4 floats
static void test_movemask_ps(void) {
  float a[4] = {-1.0f, 1.0f, -1.0f, 1.0f};
  __m128 va;
  memcpy(&va, a, 16);
  int mask = _mm_movemask_ps(va);
  CHECK(mask == 0x5, "_mm_movemask_ps"); // bits 0,2 set
}

// _mm_movemask_pi8: sign bits of 8 bytes (__m64)
static void test_movemask_pi8(void) {
  int8_t a[8] = {-1, 0, 0, -1, 0, 0, 0, 0};
  __m64 va;
  memcpy(&va, a, 8);
  int mask = _mm_movemask_pi8(va);
  CHECK((mask & 0xff) == 0x09, "_mm_movemask_pi8"); // bits 0,3 set
}

// _mm_insert_epi32: insert value 99 at lane 2
static void test_insert_epi32(void) {
  int32_t a[4] = {1, 2, 3, 4};
  __m128i va;
  memcpy(&va, a, 16);
  __m128i vc = _mm_insert_epi32(va, 99, 2);
  int32_t r[4];
  get_i32(vc, r);
  CHECK(r[0] == 1 && r[1] == 2 && r[2] == 99 && r[3] == 4,
        "_mm_insert_epi32");
}

// _mm_insert_epi8: insert value 77 at lane 5
static void test_insert_epi8(void) {
  int8_t a[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
  __m128i va;
  memcpy(&va, a, 16);
  __m128i vc = _mm_insert_epi8(va, 77, 5);
  int8_t r[16];
  get_i8(vc, r);
  CHECK(r[5] == 77 && r[4] == 5 && r[6] == 7, "_mm_insert_epi8");
}

// _mm_loadu_si16: load a single 16-bit value into lane 0
static void test_loadu_si16(void) {
  int16_t src = 0x1234;
  __m128i v = _mm_loadu_si16(&src);
  int16_t r[8];
  get_i16(v, r);
  CHECK(r[0] == 0x1234 && r[1] == 0 && r[2] == 0, "_mm_loadu_si16");
}

// _mm_loadu_si32: load a single 32-bit value into lane 0
static void test_loadu_si32(void) {
  int32_t src = 0xDEADBEEF;
  __m128i v = _mm_loadu_si32(&src);
  int32_t r[4];
  get_i32(v, r);
  CHECK(r[0] == (int32_t)0xDEADBEEF && r[1] == 0 && r[2] == 0 && r[3] == 0,
        "_mm_loadu_si32");
}

// _mm_sub_ss: subtract scalar, keep upper lanes unchanged
static void test_sub_ss(void) {
  float a[4] = {10.0f, 2.0f, 3.0f, 4.0f};
  float b[4] = {1.0f, 99.0f, 99.0f, 99.0f};
  __m128 va, vb;
  memcpy(&va, a, 16);
  memcpy(&vb, b, 16);
  __m128 vc = _mm_sub_ss(va, vb);
  float r[4];
  get_f32(vc, r);
  CHECK(r[0] == 9.0f && r[1] == 2.0f && r[2] == 3.0f && r[3] == 4.0f,
        "_mm_sub_ss");
}

// _mm_sub_sd: subtract scalar double, keep upper lane unchanged
static void test_sub_sd(void) {
  double a[2] = {10.0, 20.0};
  double b[2] = {3.0, 99.0};
  __m128d va, vb;
  memcpy(&va, a, 16);
  memcpy(&vb, b, 16);
  __m128d vc = _mm_sub_sd(va, vb);
  double r[2];
  get_f64(vc, r);
  CHECK(r[0] == 7.0 && r[1] == 20.0, "_mm_sub_sd");
}

// ── Group 2: vaaddu / vnclip rounding-mode compat ────────────────────────────

// _mm_avg_epu16: (a+b+1)>>1 for 8 uint16 lanes
static void test_avg_epu16(void) {
  uint16_t a[8] = {0, 1, 100, 200, 1000, 0xFFFE, 0xFFFF, 300};
  uint16_t b[8] = {0, 1, 101, 201, 1001, 0xFFFF, 0xFFFF, 301};
  __m128i va, vb;
  memcpy(&va, a, 16);
  memcpy(&vb, b, 16);
  __m128i vc = _mm_avg_epu16(va, vb);
  uint16_t r[8];
  get_u16(vc, r);
  uint16_t expected[8];
  for (int i = 0; i < 8; i++)
    expected[i] = (uint16_t)(((uint32_t)a[i] + b[i] + 1) >> 1);
  int ok = 1;
  for (int i = 0; i < 8; i++)
    ok &= (r[i] == expected[i]);
  CHECK(ok, "_mm_avg_epu16");
}

// _mm_avg_epu8: (a+b+1)>>1 for 16 uint8 lanes
static void test_avg_epu8(void) {
  uint8_t a[16] = {0, 1, 100, 200, 255, 254, 3, 7, 10, 20, 30, 40, 50, 60, 70,
                   80};
  uint8_t b[16] = {0, 2, 101, 201, 255, 255, 4, 8, 11, 21, 31, 41, 51, 61, 71,
                   81};
  __m128i va, vb;
  memcpy(&va, a, 16);
  memcpy(&vb, b, 16);
  __m128i vc = _mm_avg_epu8(va, vb);
  uint8_t r[16];
  get_u8(vc, r);
  int ok = 1;
  for (int i = 0; i < 16; i++)
    ok &= (r[i] == (uint8_t)(((uint32_t)a[i] + b[i] + 1) >> 1));
  CHECK(ok, "_mm_avg_epu8");
}

// _mm_packs_epi16: saturate int16→int8
static void test_packs_epi16(void) {
  int16_t a[8] = {-200, -128, -1, 0, 127, 128, 200, 1000};
  int16_t b[8] = {-1000, -129, 0, 1, 126, 127, 200, 300};
  __m128i va, vb;
  memcpy(&va, a, 16);
  memcpy(&vb, b, 16);
  __m128i vc = _mm_packs_epi16(va, vb);
  int8_t r[16];
  get_i8(vc, r);
  // a: -128,-128,-1,0,127,127,127,127
  // b: -128,-128, 0,1,126,127,127,127
  CHECK(r[0] == -128 && r[1] == -128 && r[2] == -1 && r[3] == 0 &&
            r[4] == 127 && r[5] == 127 && r[6] == 127 && r[7] == 127 &&
            r[8] == -128 && r[9] == -128 && r[10] == 0 && r[11] == 1 &&
            r[12] == 126 && r[13] == 127 && r[14] == 127 && r[15] == 127,
        "_mm_packs_epi16");
}

// _mm_packs_epi32: saturate int32→int16
static void test_packs_epi32(void) {
  int32_t a[4] = {-100000, -32768, 32767, 100000};
  int32_t b[4] = {0, 1, -1, 32768};
  __m128i va, vb;
  memcpy(&va, a, 16);
  memcpy(&vb, b, 16);
  __m128i vc = _mm_packs_epi32(va, vb);
  int16_t r[8];
  get_i16(vc, r);
  CHECK(r[0] == -32768 && r[1] == -32768 && r[2] == 32767 && r[3] == 32767 &&
            r[4] == 0 && r[5] == 1 && r[6] == -1 && r[7] == 32767,
        "_mm_packs_epi32");
}

// ── Group 3: vfcvt_x_f_v_i32m1_rm explicit rounding mode ────────────────────

// _mm_cvt_ss2si: float to int (nearest-even by default)
static void test_cvt_ss2si(void) {
  float a[4] = {2.5f, 0.0f, 0.0f, 0.0f};
  __m128 va;
  memcpy(&va, a, 16);
  int r = _mm_cvt_ss2si(va);
  // 2.5 rounds to 2 (round-to-nearest-even)
  CHECK(r == 2, "_mm_cvt_ss2si(2.5f) rounds to 2");

  float b[4] = {-3.5f, 0.0f, 0.0f, 0.0f};
  memcpy(&va, b, 16);
  r = _mm_cvt_ss2si(va);
  // -3.5 rounds to -4 (round-to-nearest-even)
  CHECK(r == -4, "_mm_cvt_ss2si(-3.5f) rounds to -4");
}

// _mm_round_ps: explicit rounding modes
static void test_round_ps(void) {
  float a[4] = {1.6f, -1.6f, 2.4f, -2.4f};
  __m128 va;
  memcpy(&va, a, 16);

  // round to nearest
  __m128 r0 = _mm_round_ps(va, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
  float f0[4];
  get_f32(r0, f0);
  CHECK(f0[0] == 2.0f && f0[1] == -2.0f && f0[2] == 2.0f && f0[3] == -2.0f,
        "_mm_round_ps nearest");

  // round toward zero (truncate)
  __m128 r1 = _mm_round_ps(va, _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC);
  float f1[4];
  get_f32(r1, f1);
  CHECK(f1[0] == 1.0f && f1[1] == -1.0f && f1[2] == 2.0f && f1[3] == -2.0f,
        "_mm_round_ps truncate");

  // round toward negative infinity (floor)
  __m128 r2 = _mm_round_ps(va, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
  float f2[4];
  get_f32(r2, f2);
  CHECK(f2[0] == 1.0f && f2[1] == -2.0f && f2[2] == 2.0f && f2[3] == -3.0f,
        "_mm_round_ps floor");

  // round toward positive infinity (ceil)
  __m128 r3 = _mm_round_ps(va, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
  float f3[4];
  get_f32(r3, f3);
  CHECK(f3[0] == 2.0f && f3[1] == -1.0f && f3[2] == 3.0f && f3[3] == -2.0f,
        "_mm_round_ps ceil");
}

// ── main ─────────────────────────────────────────────────────────────────────

int main(void) {
  printf("=== sse2rvv compat regression tests ===\n\n");

  printf("-- vlm/vsm mask load/store --\n");
  test_addsub_ps();
  test_blend_epi16();
  test_blend_pd();
  test_blend_ps();
  test_hadd_pd();
  test_hsub_pd();
  test_movemask_epi8();
  test_movemask_pd();
  test_movemask_ps();
  test_movemask_pi8();
  test_insert_epi32();
  test_insert_epi8();
  test_loadu_si16();
  test_loadu_si32();
  test_sub_ss();
  test_sub_sd();

  printf("\n-- vaaddu/vnclip rounding-mode compat --\n");
  test_avg_epu16();
  test_avg_epu8();
  test_packs_epi16();
  test_packs_epi32();

  printf("\n-- vfcvt explicit rounding-mode compat --\n");
  test_cvt_ss2si();
  test_round_ps();

  printf("\n=== %s (%d failure(s)) ===\n",
         failures == 0 ? "ALL PASSED" : "SOME FAILED", failures);
  return failures ? 1 : 0;
}
