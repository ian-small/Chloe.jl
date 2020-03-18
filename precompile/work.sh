#!/bin/bash
# from https://julialang.github.io/PackageCompiler.jl/dev/devdocs/binaries_part_2/
JULIA=/home/ianc/src/julia-1.3.1
time ${JULIA}/bin/julia --startup-file=no --trace-compile=app_precompile.jl main.jl "$@"
time ${JULIA}/bin/julia --startup-file=no -J"${JULIA}/lib/julia/sys.so" --output-o sys.o create_sysimage.jl
time gcc -shared -o sys.so -fPIC -Wl,--whole-archive sys.o -Wl,--no-whole-archive -L"${JULIA}/lib" -ljulia
rm sys.o
# can start julia with
#${JULIA}/bin/julia -Jsys.so
