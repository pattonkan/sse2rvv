# sse2rvv
![Github Actions](https://github.com/DLTcollab/sse2rvv/workflows/Github%20Actions/badge.svg?branch=master)

A C/C++ header file that converts Intel SSE intrinsics to RISCV-V Extension intrinsics.

## Introduction

`sse2rvv` is a translator of Intel SSE (Streaming SIMD Extensions) intrinsics
to [RISCV-V Extension](https://github.com/riscv/riscv-v-spec),
shortening the time needed to get an RISCV working program that then can be used to
extract profiles and to identify hot paths in the code.
The header file `sse2rvv.h` contains several of the functions provided by Intel
intrinsic headers such as `<xmmintrin.h>`, only implemented with RISCV-based counterparts
to produce the exact semantics of the intrinsics.

## Mapping and Coverage

Header file | Extension |
---|---|
`<mmintrin.h>` | MMX |
`<xmmintrin.h>` | SSE |
`<emmintrin.h>` | SSE2 |
`<pmmintrin.h>` | SSE3 |
`<tmmintrin.h>` | SSSE3 |
`<smmintrin.h>` | SSE4.1 |
`<nmmintrin.h>` | SSE4.2 |
`<wmmintrin.h>` | AES  |

`sse2rvv` aims to support SSE, SSE2, SSE3, SSSE3, SSE4.1, SSE4.2 and AES extension.

In order to deliver RVV-equivalent intrinsics for all SSE intrinsics used widely,
please be aware that some SSE intrinsics exist a direct mapping with a concrete
NEON-equivalent intrinsic. Others, unfortunately, lack a 1:1 mapping, meaning that
their equivalents are built utilizing a number of NEON intrinsics.

For example, SSE intrinsic `_mm_loadu_si128` has a direct RVV mapping (`vld1q_s32`),
but SSE intrinsic `_mm_maddubs_epi16` has to be implemented with multiple RVV instructions.

### Floating-point compatibility

Some conversions require several RVV intrinsics, which may produce inconsistent results
compared to their SSE counterparts due to differences in the arithmetic rules of IEEE-754.

Taking a possible conversion of `_mm_rsqrt_ps` as example:

```c
__m128 _mm_rsqrt_ps(__m128 in)
{
    float32x4_t out = vrsqrteq_f32(vreinterpretq_f32_m128(in));

    out = vmulq_f32(
        out, vrsqrtsq_f32(vmulq_f32(vreinterpretq_f32_m128(in), out), out));

    return vreinterpretq_m128_f32(out);
}
```

The `_mm_rsqrt_ps` conversion will produce NaN if a source value is `0.0` (first INF for the
reciprocal square root of `0.0`, then INF * `0.0` using `vmulq_f32`). In contrast,
the SSE counterpart produces INF if a source value is `0.0`.
As a result, additional treatments should be applied to ensure consistency between the conversion and its SSE counterpart.

## Usage

- Put the file `sse2rvv.h` in to your source code directory.

- Locate the following SSE header files included in the code:
```C
#include <xmmintrin.h>
#include <emmintrin.h>
```
  {p,t,s,n,w}mmintrin.h could be replaceable as well.

- Replace them with:
```C
#include "sse2rvv.h"
```

- Explicitly specify platform-specific options to gcc/clang compilers.
  * On riscv64
  ```shell
  -march=r64gcv_zba
  ```

## Compile-time Configurations

Though floating-point operations in NEON use the IEEE single-precision format, NEON does not fully comply to the IEEE standard when inputs or results are denormal or NaN values for minimizing power consumption as well as maximizing performance.
Considering the balance between correctness and performance, `sse2rvv` recognizes the following compile-time configurations:
* `SSE2RVV_PRECISE_MINMAX`: Enable precise implementation of `_mm_min_{ps,pd}` and `_mm_max_{ps,pd}`. If you need consistent results such as handling with NaN values, enable it.
* `SSE2RVV_PRECISE_DIV` (deprecated): Enable precise implementation of `_mm_rcp_ps` and `_mm_div_ps` by additional Netwon-Raphson iteration for accuracy.
* `SSE2RVV_PRECISE_SQRT` (deprecated): Enable precise implementation of `_mm_sqrt_ps` and `_mm_rsqrt_ps` by additional Netwon-Raphson iteration for accuracy.
* `SSE2RVV_PRECISE_DP`: Enable precise implementation of `_mm_dp_pd`. When the conditional bit is not set, the corresponding multiplication would not be executed.

The above are turned off by default, and you should define the corresponding macro(s) as `1` before including `sse2rvv.h` if you need the precise implementations.

## Run Built-in Test Suite

`sse2rvv` provides a unified interface for developing test cases. These test
cases are located in `tests` directory, and the input data is specified at
runtime. Use the following commands to perform test cases:
```shell
$ make test
```

For running test with enabling features, you can use assign the features with `FEATURE` command.
If `none` is assigned, then the command will be the same as simply calling `make test`.
The following command enable `crypto` and `crc` features in the tests.
```
$ make FEATURE=crypto+crc test
```

You can specify GNU toolchain for cross compilation as well.
[QEMU](https://www.qemu.org/) should be installed in advance.
```shell
$ make CROSS_COMPILE=aarch64-linux-gnu- test # ARMv8-A
```
or
```shell
$ make CROSS_COMPILE=arm-linux-gnueabihf- test # ARMv7-A
```

Check the details via [Test Suite for SSE2RVV](tests/README.md).

## Reference
* [sse2neon](https://github.com/DLTcollab/sse2neon)
* [Intel Intrinsics Guide](https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html)
* [Microsoft: x86 intrinsics list](https://learn.microsoft.com/en-us/cpp/intrinsics/x86-intrinsics-list)
* [Arm Neon Intrinsics Reference](https://developer.arm.com/architectures/instruction-sets/simd-isas/neon/intrinsics)
* [Neon Programmer's Guide for Armv8-A](https://developer.arm.com/architectures/instruction-sets/simd-isas/neon/neon-programmers-guide-for-armv8-a)
* [NEON Programmer's Guide](https://static.docs.arm.com/den0018/a/DEN0018A_neon_programmers_guide_en.pdf)
* [qemu/target/i386/ops_sse.h](https://github.com/qemu/qemu/blob/master/target/i386/ops_sse.h): Comprehensive SSE instruction emulation in C. Ideal for semantic checks.
* [Porting Takua Renderer to 64-bit ARM- Part 1](https://blog.yiningkarlli.com/2021/05/porting-takua-to-arm-pt1.html)
* [Porting Takua Renderer to 64-bit ARM- Part 2](https://blog.yiningkarlli.com/2021/07/porting-takua-to-arm-pt2.html)
* [Comparing SIMD on x86-64 and arm64](https://blog.yiningkarlli.com/2021/09/neon-vs-sse.html)
* [Port with SSE2Neon and SIMDe](https://developer.arm.com/documentation/102581/0200/Port-with-SSE2Neon-and-SIMDe)
* [Genomics: Optimizing the BWA aligner for Arm Servers](https://community.arm.com/arm-community-blogs/b/high-performance-computing-blog/posts/optimizing-genomics-and-the-bwa-aligner-for-arm-servers)
* [Bit twiddling with Arm Neon: beating SSE movemasks, counting bits and more](https://community.arm.com/arm-community-blogs/b/infrastructure-solutions-blog/posts/porting-x86-vector-bitmask-optimizations-to-arm-neon)
* [C/C++ on Graviton](https://github.com/aws/aws-graviton-getting-started/blob/main/c-c%2B%2B.md)

## Licensing

`sse2rvv` is freely redistributable under the MIT License.
