# Sets a number of global variables for use in debugging, such as unit tests.
module EvInit
export ev_init,N,num_subpops,num_attributes,ngens,num_emmigrants,mutation_stddev,move_range,move_time_interval,opt_loss_cutoff,
    horiz_select,uniform_start,min_fit,linear_fitness
#=
run_evolve(N,num_attributes,num_subpops, ngens,mutation_stddev,num_emmigrants,move_range,move_time_interval,opt_loss_cutoff,uniform_start,min_fit=min_fit,
      linear_fitness=linear_fitness)
=#
using TemporalEvolution
function ev_init()
  #include("types.jl")
  #include("../src/propsel.jl")
  #include("../src/fitness.jl")
  global num_trials=1
  global N=27
  global num_subpops = 9
  global num_attributes = 3
  global ngens = 10
  global num_emmigrants = 2
  global mutation_stddev = 0.04
  global move_range = 0.1
  global move_time_interval = 5
  #global opt_loss_cutoff = 0.31
  global horiz_select = false
  global uniform_start = false
  global min_fit = 0.3
  global topology="vonneumann"   # must be one of "circular", "ring", "vonneumann", "moore", or "global" 
  global linear_fitness=true
  global linfit_slope = 1.0
  global tr = TemporalEvolution.temporal_result(num_trials, N, num_attributes, num_subpops, ngens, mutation_stddev, num_emmigrants, 
    move_range, move_time_interval, horiz_select, min_fit, topology=topology, uniform_start=uniform_start, linear_fitness=linear_fitness, linfit_slope=linfit_slope )
  global ideal = fill( tr.ideal_init, tr.num_attributes )
  global vt = Dict{Int64,variant_type}()
  tr
end  

#evolve(tr)
 
end

