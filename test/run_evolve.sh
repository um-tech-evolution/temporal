# bash script to run config examples
# These runs produce files  src/configs/example1.csv and src/configs/example2.csv.
cd ../src
julia -p 2 -L TemporalEvolution.jl run.jl configs/example1
julia -p 2 -L TemporalEvolution.jl run.jl configs/example2
julia -p 2 -L TemporalEvolution.jl run.jl configs/example3
julia -p 2 -L TemporalEvolution.jl run.jl configs/example4
julia -p 2 -L TemporalEvolution.jl run.jl configs/example5
