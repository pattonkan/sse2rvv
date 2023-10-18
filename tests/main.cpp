#include <stdint.h>
#include <stdio.h>
#include "impl.h"

int main(int /*argc*/, const char ** /*argv*/)
{
    SSE2RVV::SSE2RVV_TEST *test = SSE2RVV::SSE2RVV_TEST::create();
    uint32_t pass_count = 0;
    uint32_t failed_count = 0;
    uint32_t ignore_count = 0;
    for (uint32_t i = 0; i < SSE2RVV::it_last; i++) {
        SSE2RVV::INSTRUCTION_TEST it = SSE2RVV::INSTRUCTION_TEST(i);
        SSE2RVV::result_t ret = test->run_test(it);
        // If the test fails, we will run it again so we can step into the
        // debugger and figure out why!
        if (ret == SSE2RVV::TEST_FAIL) {
            printf("Test %-30s failed\n", SSE2RVV::instruction_string[it]);
            failed_count++;
        } else if (ret == SSE2RVV::TEST_UNIMPL) {
            printf("Test %-30s skipped\n", SSE2RVV::instruction_string[it]);
            ignore_count++;
        } else {
            printf("Test %-30s passed\n", SSE2RVV::instruction_string[it]);
            pass_count++;
        }
    }
    test->release();
    printf(
        "SSE2RVV_TEST Complete!\n"
        "Passed:  %d\n"
        "Failed:  %d\n"
        "Ignored: %d\n"
        "Coverage rate: %.2f%%\n",
        pass_count, failed_count, ignore_count,
        (float) pass_count / (pass_count + failed_count + ignore_count) * 100);

    return failed_count ? -1 : 0;
}
