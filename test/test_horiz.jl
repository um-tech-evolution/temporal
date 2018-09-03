#! /usr/local/julia/bin/julia
# Unit tests for functions in horiz.jl
# This test program enables running horiz_transfer() from horiz.jl and the component functions
#  horiz_transfer_by_fitness() and new_emigrants_funct() separately.
#  Other functions from horiz.jl could also be called.
#  Tests are done that the results are equal to my current results.
#  Some tests may depend on the random number generator generating the same sequence of 
#     random numbers when started with the same seed.

include("../src/TemporalEvolution.jl")
using FactCheck

function test_horiz( paramd::param_type )
  #println("N: ",paramd[:N])
  ideal = fill( paramd[:ideal_init], paramd[:num_attributes] )
  subpop_size = Int(floor(paramd[:N]/paramd[:num_subpops]))
  vt = Dict{Int64,temporal_variant_type}()
  id = [0]
  #meta_pop = init_meta_pop( paramd, vt, ideal, id )  # 
  meta_pop = init_meta_pop_test( vt, ideal, id )
  #print_meta_pop( paramd, meta_pop, vt )
  #print_meta_pop_attributes( paramd, meta_pop, vt )
  mmeans = fit_means( meta_pop, vt )
  #println("mmeans: ",mmeans)
  facts("testing mmeans") do
    @fact [0.43964466094067267,0.41242365978762213,0.43858807184546883,0.4135994505535974] --> roughly(mmeans)
  end
  source_subpop_list = horiz_transfer_by_fitness!( meta_pop, paramd, vt, ideal, mmeans, id, topology=paramd[:topology], neg_select=false, emigrant_select=false )
  #println("source subpop_list: ",source_subpop_list)
  id = [0]
  meta_pop = init_meta_pop_test( vt, ideal, id )
  source_subpop_list = [2,1,1,1]
  srand(1)
  new_emigrants = new_emigrants_funct(meta_pop, paramd, vt, source_subpop_list, ideal, id, emigrant_select=false )
  #println("new emigrants: ",new_emigrants)
  facts("testing new emigrants") do
    @fact new_emigrants --> Array{Int64,1}[[8], [9], [10], [11]]
  end
  g = 1
  srand(1)
  horiz_transfer( meta_pop, paramd, vt, ideal, mmeans, id, g )
  mmeans = fit_means( meta_pop, vt )
  #println("mmeans: ",mmeans)
  #=
  print("[")
  for m in mmeans
    print(m,",")
  end
  println("]")
  =#
  facts("testing mmeans") do
    @fact [0.43964466094067267,0.4369211344706805,0.42559341256113936,0.4265657954113531] --> roughly(mmeans)
  end
  #print_meta_pop( paramd, meta_pop, vt )
  #print_meta_pop_attributes( paramd, meta_pop, vt )
  facts("testing meta_pop") do
    @fact meta_pop --> Array{Int64,1}[[1, 2], [3, 12], [6, 13], [8, 14]]
  end
end

# Initialize meta_pop with specified attribute values for testing
function init_meta_pop_test(vt::Dict{Int64,temporal_variant_type}, ideal::Vector{Float64}, id::Vector{Int64} )
  N = 8
  num_subpops = 4
  subpop_size = Int(floor(N/num_subpops))
  attr = [[0.5,0.45],[0.45,0.45],[0.43,0.47],[0.43,0.43],[0.48,0.46],[0.44,0.45],[0.43,0.48],[0.42,0.43]]
  meta_pop = [ (j-1)*subpop_size+collect(1:subpop_size)  for j = 1:num_subpops ]
  for i = 1:N
    vt[i] = temporal_variant_type( fitness( attr[i], ideal, minFit=0.4, linfit_slope=1.0 ), attr[i], 0 )
  end
  id[1] += N
  return meta_pop
end
  
paramd = init_dictionary( TemporalEvolution.temporal_param_fields )
paramd[:simname] = "test_horiz"
paramd = read_parameter_file( "horiz_params.jl", paramd )
context("Testing horiz.jl") do
  test_horiz( paramd )
end
