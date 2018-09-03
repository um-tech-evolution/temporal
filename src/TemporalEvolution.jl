module TemporalEvolution
include("../src/types.jl")
include("../src/propsel.jl")
include("../src/evolve.jl")
include("../src/subpop_alive.jl")
include("../src/run_temporal.jl")
include("../src/fitness.jl")
include("../src/horiz.jl")
include("../src/mutate.jl")    # Also in evotech/unified/src as a symlink
include("../src/stats.jl")
include("../src/uni.jl")       # Also in evotech/unified/src as a symlink
include("../src/spatial.jl")
end
using TemporalEvolution


