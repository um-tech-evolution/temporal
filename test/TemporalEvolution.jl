module TemporalEvolution
include("../src/types.jl")
include("../src/propsel.jl")
include("../src/evolve.jl")
include("../src/run_temporal.jl")
include("../test/test_fitness.jl")
include("../src/horiz.jl")
end
using TemporalEvolution
include("../src/ev_init.jl")
using EvInit


