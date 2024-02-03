# sse2rvv

A C/C++ header file that converts Intel SSE intrinsics to RISCV-V Extension intrinsics.

## Introduction

`sse2rvv` is a translator of Intel SSE (Streaming SIMD Extensions) intrinsics
to [RISCV-V Extension](https://github.com/riscv/riscv-v-spec),
shortening the time needed to get an RISCV working program that then can be used to
extract profiles and to identify hot paths in the code.
The header file `sse2rvv.h` contains several of the functions provided by Intel
intrinsic headers such as `<xmmintrin.h>`, only implemented with RISCV-based counterparts
to produce the exact semantics of the intrinsics.

This project is based on [sse2neon](https://github.com/DLTcollab/sse2neon), and modify it to RISCV version.

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
RVV-equivalent intrinsic. Others, unfortunately, lack a 1:1 mapping, meaning that
their equivalents are built utilizing a number of RVV intrinsics.

For example, SSE intrinsic `_mm_add_epi16` has a direct RVV mapping (`__riscv_vadd_vv_i16m1`),
but SSE intrinsic `_mm_maddubs_epi16` has to be implemented with multiple RVV instructions.

### Floating-point compatibility

Some conversions require several RVV intrinsics, which may produce inconsistent results
compared to their SSE counterparts due to differences in the arithmetic rules of IEEE-754.

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

## Run Built-in Test Suite

`sse2rvv` provides a unified interface for developing test cases. These test
cases are located in `tests` directory, and the input data is specified at
runtime. Use the following commands to perform test cases:
```shell
$ make test
```

## Reference

* [sse2neon](https://github.com/DLTcollab/sse2neon)
* [neon2rvv](https://github.com/howjmay/neon2rvv)
* [Intel Intrinsics Guide](https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html)
* [Microsoft: x86 intrinsics list](https://learn.microsoft.com/en-us/cpp/intrinsics/x86-intrinsics-list)
* [riscv-v-spec](https://github.com/riscv/riscv-v-spec)
* [rvv-intrinsic-doc](https://github.com/riscv-non-isa/rvv-intrinsic-doc)
* [riscv-c-api](https://github.com/riscv-non-isa/riscv-c-api-doc/blob/master/riscv-c-api.md)

## Licensing

`sse2rvv` is freely redistributable under the MIT License.
