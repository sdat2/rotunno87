# !/usr/bin/bash
# usage: sh ./compile.sh
gfortran -std=f95 -Wextra -Wall -pedantic -frecursive -Ofast -fdump-parse-tree -O3 -floop-parallelize-all hurricane.f -o hurricane5.exe
# -Mipa=fast  -fsyntax-only
# Mconcur
# g77 hurricane.f -o hurricane2.exe
