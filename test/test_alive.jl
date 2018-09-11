using FactCheck
include("../src/TemporalEvolution.jl")
# Sample run:  julia test_alive.jl

function test_alive(paramd::param_type)
  id = [0]
  paramd[:mutStddev] = 0.03
  paramd[:N] = 16
  paramd[:num_subpops] = 4
  subpop_size = Int(floor(paramd[:N]/paramd[:num_subpops]))
  ideal = fill( paramd[:ideal_init], paramd[:num_attributes] )
  vt = Dict{Int64,temporal_variant_type}()
  meta_pop = init_meta_pop( paramd, vt, ideal, id )
  for s in meta_pop
    r = count_alive( s, paramd, vt )
    println("subpop: ",s,"  result: ",r)
  end
  alive = subpops_alive( meta_pop, paramd, vt )
  println("alive: ",alive)
  prev_gen =      [true,true,false,false]
  after_propsel = [true,false,true,false]
  after_horiz =   [false,false,true,true]
  (propsel_loss,propsel_gain,horiz_loss,horiz_gain) = 
      opt_gained_lost_counts( prev_gen, after_propsel, after_horiz )
  println("(propsel_loss,propsel_gain,horiz_loss,horiz_gain):",
           (propsel_loss,propsel_gain,horiz_loss,horiz_gain))
end

srand(1)   # set random number seed to 1
paramd = init_dictionary( TemporalEvolution.temporal_param_fields )
paramd[:simname] = "test_mutate"
paramd = read_parameter_file( "mutate_params.jl", paramd )

context("running test_alive()") do
  test_alive(paramd)
end

