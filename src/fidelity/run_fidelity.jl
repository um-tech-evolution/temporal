export print_fidelity, writeheader, writerow, horiz_param_check
#=
Recommended command line to run:
>  julia -L TemporalEvolution.jl run_fidelity.jl examples/example1
=#
export  writeheader, writerow
#include("types.jl")
  
function print_fidelity_type( fd::fidelity_type )
  println("num_trials: ", fd.T)
  println("N: ", fd.N)
  println("s: ", fd.s)
  println("q: ", fd.q)
  println("ngens: ", fd.ngens)
  println("init_high_freq: ", fd.init_high_freq)
  println("psel_first: ", fd.psel_first)
  println("high_fixed_count: ", fd.high_fixed_count)
  println("low_fixed_count: ", fd.low_fixed_count)
  println("gen_reached_count: ", fd.gen_reached_count)
  println("average_high: ", fd.average_high)
end

function writeheader( stream::IO, fd::fidelity_type, N_list )
  param_strings = [
    "# $(string(Dates.today()))",
    "# num_trials=$(fd.T)",
    "# N_list=$(N_list)",
    "# s=$(fd.s)",
    "# q=$(fd.q)",
    "# ngens=$(fd.ngens)",
    "# init_high_freq=$(fd.init_high_freq)",
    "# psel_first=$(fd.psel_first)",
  ]
  write(stream,join(param_strings,"\n"),"\n")
  heads = [
    "trial",
    "N",
    "high_fixed_counts",
    "low_fixed_counts",
    "gen_reached_count",
    "average_high"
    ]
  write(stream,join(heads,","),"\n")
end
    
function writerow( stream::IO, trial::Int64, fd::fidelity_type )
  line = Any[
    trial,
    fd.N,
    fd.high_fixed_count,
    fd.low_fixed_count,
    fd.gen_reached_count,
    fd.average_high
  ]
  write(stream,join(line,","),"\n")
end

