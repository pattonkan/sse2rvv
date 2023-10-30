#ifndef SSE2RVV_BINDING_H
#define SSE2RVV_BINDING_H

#include <stdlib.h>

// The SSE2RVV unit tests run both within our own internal project
// as well as within the open source framework.
// This header file is used to abstract any distinctions between
// those two build environments.
//
// Initially, this is for how 16 byte aligned memory is allocated
namespace SSE2RVV {
void *platform_aligned_alloc(size_t size);
void platform_aligned_free(void *ptr);

} // namespace SSE2RVV

#endif
