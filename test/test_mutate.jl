using FactCheck
include("../src/TemporalEvolution.jl")
# Test the mutate_meta_pop!() function
# Sample run:  julia5 test_mutate.jl
# Sample run:  julia -L ../src/TemporalEvolution.jl test_mutate.jl  # Failed 1/29/18

function test_mutate(tr)
  id = [0]
  #tr = ev_init()
  tr[:mutStddev] = 0.03
  tr[:N] = 256
  tr[:num_subpops] = 1
  subpop_size = Int(floor(tr[:N]/tr[:num_subpops]))
  ideal = fill( tr[:ideal_init], tr[:num_attributes] )
  vt = Dict{Int64,temporal_variant_type}()
  meta_pop = init_meta_pop( tr, vt, ideal, id )
  gen_innov_counts = mutate_meta_pop!( meta_pop, vt, ideal, id, tr )
  mmeans = fit_means( meta_pop, vt )
  facts("  >  testing mutate_meta_pop!") do
    println("mmeans: ",mmeans)
    @fact mmeans[1] --> roughly(0.40100530022116637)
  end
  # Test gen_innov_counts
  num_pos = 0
  num_half_pos = 0
  num_neg = 0
  num_half_neg = 0
  for j = 1:tr[:num_subpops]
    s = zeros(Float64,subpop_size)
    for i = 1:subpop_size
      fitness_i = fitness( vt[meta_pop[j][i]].attributes, ideal, minFit=tr[:minFit],
          linfit_slope=tr[:linfit_slope] )
      s[i] = fitness_i/mmeans[j]
      #println("j: ",j,"  i: ",i,"  fit_i: ",fitness_i,"  s[i]: ",s[i])
      num_pos += s[i] > 1.0+1.0/subpop_size ? 1 : 0
      num_half_pos += s[i] > 1.0+0.5/subpop_size ? 1 : 0
      num_neg += s[i] < 1.0-1.0/subpop_size ? 1 : 0
      num_half_neg += s[i] < 1.0-0.5/subpop_size ? 1 : 0
    end
  end    
  #println("     pos: ",gen_innov_counts.pos,"  ",num_pos)
  #println("half pos: ",gen_innov_counts.half_pos,"  ",num_half_pos)
  #println("     neg: ",gen_innov_counts.neg,"  ",num_neg)
  #println("half neg: ",gen_innov_counts.half_neg,"  ",num_half_neg)
  facts("  >  testing gen_innov_counts") do
    @fact gen_innov_counts.pos --> num_pos
    @fact gen_innov_counts.half_pos --> num_half_pos
    @fact gen_innov_counts.neg --> num_neg
    @fact gen_innov_counts.half_neg --> num_half_neg
  end
end

srand(1)   # set random number seed to 1
paramd = init_dictionary( TemporalEvolution.temporal_param_fields )
paramd[:simname] = "test_mutate"
paramd = read_parameter_file( "mutate_params.jl", paramd )

context("Testing mutation in evolve.jl") do
  test_mutate(paramd)
end

