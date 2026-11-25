#!/usr/bin/env bash
# Build and run the 1D fin solver in every installed language, then regenerate figures/tables.
H=${1:-60}
command -v gcc   >/dev/null && gcc -O2 fin.c -o fin_c -lm && echo "C    : $(./fin_c $H) W"
command -v g++   >/dev/null && g++ -O2 fin.cpp -o fin_cpp && echo "C++  : $(./fin_cpp $H) W"
command -v rustc >/dev/null && rustc -O fin.rs -o fin_rs 2>/dev/null && echo "Rust : $(./fin_rs $H) W"
command -v go    >/dev/null && go run fin.go $H | sed 's/^/Go   : /'
python3 engine_fin.py
python3 engine_fin_advanced.py
python3 engine_fin_optimization.py
python3 test_engine_fin.py
