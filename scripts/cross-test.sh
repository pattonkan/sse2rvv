#!/usr/bin/env bash

. .ci/common.sh

check_platform

# Clang/LLVM is natively a cross-compiler.
# TODO: Do cross-compilation using Clang
# https://clang.llvm.org/docs/CrossCompilation.html
if [ $(printenv CXX | grep clang) ]; then
    exit
fi

set -x

make clean
make CROSS_COMPILE=riscv64-unknown-elf- build-test || exit 1 # riscv64
