# Unit test of src/stats.jl
# This file give examples of the functions fit_means_vars_stddevs_cfvars
#    and function attr_means_vars_stddevs_cfvars
using Base.Test
using FactCheck
include("../src/TemporalEvolution.jl")

function setup_test()
  vt = Dict{Int64,temporal_variant_type}()
  subpops = [[1,2,3],[4,5,6]]
  vt[1] = temporal_variant_type( 3.0, [1.0, 2.0], 0 )
  vt[2] = temporal_variant_type( 4.0, [3.0, 1.0], 0 )
  vt[3] = temporal_variant_type( 8.0, [3.0, 5.0], 0 )
  vt[4] = temporal_variant_type( 3.0, [0.0, 3.0], 0 )
  vt[5] = temporal_variant_type( 9.0, [5.0, 4.0], 0 )
  vt[6] = temporal_variant_type( 6.0, [4.0, 2.0], 0 )
  return subpops, vt
end

# Same as the built-in function var
variance(lst) = sum(map(x->x^2,lst .- mean(lst)))/(length(lst)-1)

# Same as the built-in function std
std_dev(lst) = sqrt(variance(lst))

coefvar(lst) = std_dev(lst)/mean(lst)

context("Testing stats.jl") do
  (subpops, vt) = setup_test()
  facts("testing fitness stats") do
    ( mmeans, vvars, stddevs, cfvars ) = fit_means_vars_stddevs_cfvars( subpops, vt )
    #println([variance([vt[i].fitness for i in s]) for s in subpops])
    #println("fitness: means: ", mmeans, "  variances: ", vvars, "  stddevs: ", stddevs, "  coef vars: ", cfvars  )
    @fact [mean([vt[i].fitness for i in s]) for s in subpops] --> roughly(mmeans)
    @fact [variance([vt[i].fitness for i in s]) for s in subpops] --> roughly(vvars)
    @fact [std_dev([vt[i].fitness for i in s]) for s in subpops] --> roughly(stddevs)
    @fact [coefvar([vt[i].fitness for i in s]) for s in subpops] --> roughly(cfvars)
  end
  
  facts("testing attribute stats") do
    ( mmeans, vvars, stddevs, cfvars ) = attr_means_vars_stddevs_cfvars( subpops, vt )
    #println("attribs: means: ", mmeans, "  variances: ", vvars, "  stddevs: ", stddevs, "  coef vars: ", cfvars  )
    a1_1 = [1, 3, 3]; m1_1 = mean(a1_1); v1_1 = variance(a1_1); s1_1 = std_dev(a1_1); c1_1 = coefvar(a1_1)
    a1_2 = [0, 5, 4]; m1_2 = mean(a1_2); v1_2 = variance(a1_2); s1_2 = std_dev(a1_2); c1_2 = coefvar(a1_2)
    a2_1 = [2, 1, 5]; m2_1 = mean(a2_1); v2_1 = variance(a2_1); s2_1 = std_dev(a2_1); c2_1 = coefvar(a2_1)
    a2_2 = [3, 4, 2]; m2_2 = mean(a2_2); v2_2 = variance(a2_2); s2_2 = std_dev(a2_2); c2_2 = coefvar(a2_2)
    #println( [(m1_1+m2_1)/2, (m1_2+m2_2)/2] )
    #println( [(v1_1+v2_1)/2, (v1_2+v2_2)/2] )
    #println( [(s1_1+s2_1)/2, (s1_2+s2_2)/2] )
    #println( [(c1_1+c2_1)/2, (c1_2+c2_2)/2] )
    @fact [(m1_1+m2_1)/2, (m1_2+m2_2)/2] --> roughly(mmeans)
    @fact [(v1_1+v2_1)/2, (v1_2+v2_2)/2] --> roughly(vvars)
    @fact [(s1_1+s2_1)/2, (s1_2+s2_2)/2] --> roughly(stddevs)
    @fact [(c1_1+c2_1)/2, (c1_2+c2_2)/2] --> roughly(cfvars)
  end
end

