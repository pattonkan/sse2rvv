#!/usr/bin/env bash

set -x

for file in ${SOURCES};
do
    clang-format ${file} > expected-format
    diff -u -p --label="${file}" --label="expected coding style" ${file} expected-format
done
exit $(clang-format --output-replacements-xml ${SOURCES} | egrep -c "</replacement>")
