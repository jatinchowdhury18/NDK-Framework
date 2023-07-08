#!/usr/bin/env bash

echo "Generating C++ code..."
python3 generate_ndk_cpp.py test_configs/cry_baby_ndk_config.json
mv CryBabyNDK.* tests/cry_baby_test
python3 generate_ndk_cpp.py test_configs/big_muff_2d_config.json
mv BigMuff2D.* tests/big_muff_test
python3 generate_ndk_cpp.py test_configs/big_muff_dp_config.json
mv BigMuffDP.* tests/big_muff_test
python3 generate_ndk_cpp.py test_configs/rc_lowpass_ndk_config.json
mv RCLowpassNDK.* tests/rc_lowpass_test
python3 generate_ndk_cpp.py test_configs/sk_lowpass_ndk_config.json
mv SKLowpassNDK.* tests/sk_lowpass_test

echo "Configuring tests..."
cmake -Bbuild-tests -Stests -DCMAKE_BUILD_TYPE=Release

echo "Running Cry Baby test..."
cmake --build build-tests --config Release --target cry_baby_test --parallel
./build-tests/cry_baby_test/cry_baby_test

echo "Running Big Muff test..."
cmake --build build-tests --config Release --target big_muff_test --parallel
./build-tests/big_muff_test/big_muff_test

echo "Running RC-Lowpass test..."
cmake --build build-tests --config Release --target rc_lowpass_test --parallel
./build-tests/rc_lowpass_test/rc_lowpass_test

echo "Running Sallen-Key Lowpass test..."
cmake --build build-tests --config Release --target sk_lowpass_test --parallel
./build-tests/sk_lowpass_test/sk_lowpass_test
