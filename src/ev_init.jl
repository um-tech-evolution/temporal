# Sets a number of global variables for use in debugging, such as unit tests.
module EvInit
export ev_init,N,num_subpops,num_attributes,ngens,num_emmigrants,mutStddev,move_range,move_time_interval,opt_loss_cutoff,
    horiz_select,uniform_start,minFit,linear_fitness
#=
run_evolve(N,num_attributes,num_subpops, ngens,mutStddev,num_emmigrants,move_range,move_time_interval,opt_loss_cutoff,uniform_start,minFit=minFit,
      linear_fitness=linear_fitness)
=#
using TemporalEvolution
function ev_init()
  #include("types.jl")
  #include("../src/propsel.jl")
  #include("../src/fitness.jl")
  global simtype=2
  global num_trials=1
  global N=64
  global num_subpops = 8
  global num_attributes = 1
  global ngens = 60
  global num_emmigrants = 0
  global mutStddev = 0.015
  global move_range = 0.07
  global move_time_interval = 5
  global horiz_select = false
  global uniform_start = true
  global minFit = 0.46
  global topology="ring"   # must be one of "circular", "ring", "vonneumann", "moore", or "global" 
  global linear_fitness=true
  global linfit_slope = 1.0
  global burn_in = 0.0
  global tr = TemporalEvolution.temporal_result(simtype, num_trials, N, num_attributes, num_subpops, ngens, mutStddev, num_emmigrants, 
    move_range, move_time_interval, horiz_select, minFit, topology=topology, uniform_start=uniform_start, linear_fitness=linear_fitness, linfit_slope=linfit_slope, burn_in=burn_in )
  global ideal = fill( tr.ideal_init, tr.num_attributes )
  global vt = Dict{Int64,variant_type}()
  tr
end  

#evolve(tr)
 
end
#srand(1)

