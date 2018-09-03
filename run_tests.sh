#! /usr/bin/bash

cd src
julia run.jl examples/example1 1
julia run.jl examples/example2 1
julia run.jl examples/example3 1
julia run.jl examples/example4 1
julia run.jl examples/example5 1
julia run.jl examples/example6 1
julia run.jl examples/example7 1
julia run.jl examples/example8 1
julia run.jl examples/sp_example1 1
julia run.jl examples/sp_example2 1
julia run.jl examples/sp_example3 1
cd ../test
julia test_mutate.jl
julia test_fitness.jl
julia test_horiz.jl
julia test_stats.jl
cd ..
