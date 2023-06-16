#!/usr/bin/env bash

echo "Generating C++ code..."
python3 generate_ndk_cpp.py test_configs/cry_baby_ndk_config.json
mv CryBabyNDK.* tests/cry_baby_test

echo "Configuring tests..."
cmake -Bbuild-tests -Stests

echo "Running Cry Baby test..."
cmake --build build-tests --target cry_baby_test --parallel
./build-tests/cry_baby_test/cry_baby_test
