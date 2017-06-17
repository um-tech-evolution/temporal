module Fidelity
include("types.jl")
include("propsel_retention.jl")
include("run_fidelity.jl")
end
using Fidelity
module FdInit
include("fd_init.jl")
end
using FdInit


