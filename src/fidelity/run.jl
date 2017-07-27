#include("types.jl")
#include("run_fidelity.jl")

function run_trials( simname::AbstractString )
  stream = open("$(simname).csv","w")
  println("stream: ",stream)
  fd = fidelity( T, N_list[1], q_list[1], s_list[1], ngens, init_high_freq, psel_first, deterministic_mutation_list[1] )
  fd_list_run = fidelity_type[]
  for s in s_list
    for q in q_list
      for N in N_list
        for deterministic_mutation in deterministic_mutation_list
          fd = fidelity( T, N, q, s, ngens, init_high_freq, psel_first, deterministic_mutation )
          Base.push!( fd_list_run, fd )
        end
      end
    end
  end
  fd_list_result = pmap( count_high_low_fixed, fd_list_run )
  writeheader( STDOUT, fd, N_list )
  writeheader( stream, fd, N_list )
  trial = 1
  for fd_result in fd_list_result
    writerow(STDOUT, trial, fd_result)
    writerow(stream, trial, fd_result)
    trial += 1
  end
end


if length(ARGS) > 0
  simname = ARGS[1]
else
  simname = "examples/example1"
end
include("$(simname).jl")
println("simname: ",simname)
#fd = fidelity(T,N,q,s,ngens,init_high_freq,psel_first)
#println("fd: ",fd)
run_trials( simname )

#=
#(g,k) = propsel_first(N,q,s,ngens,init_high_freq)
#println("final g: ",g,"  k: ",k)
#(hf,lf,gr) = count_high_low_fixed( T, N,q,s,ngens,init_high_freq,psel_first)
fd = count_high_low_fixed( fd )
println( "high fixed_count: ",fd.high_fixed_count)
println( "low fixed_count: ",fd.low_fixed_count)
println( "gen_reached count: ",fd.gen_reached_count)
if fd.gen_reached_count > 0
  println("average high count: ",fd.average_high)
end
=#
