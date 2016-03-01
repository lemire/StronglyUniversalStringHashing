#!/bin/bash

set -o errexit

declare -a TESTS

TESTS=(bitset_container_unit \
       array_container_unit \
       run_container_unit \
       toplevel_unit \
       unit)

function run_tests() {
  for test in ./test/correctness/*.exe; do $test; done
}

function main() {
  echo -e " \x1B[0;32mCompiling...\x1B[0m "
  make --silent  clean all
  echo -e " \x1B[0;32mRunning tests...\x1B[0m "
  run_tests

  echo -e "\n\n \x1B[0;31m[\x1B[0m \x1B[0;32mAll tests clear.\x1B[0m \x1B[0;31m]\x1B[0m \n\n"
}

main
