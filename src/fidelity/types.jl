export fidelity_type, fidelity

type fidelity_type 
  T::Int64
  N::Int64
  s::Float64
  q::Float64
  ngens::Int64
  init_high_freq::Float64
  psel_first::Bool               # prop selection first
  deterministic_mutation::Bool  # do mutation using q without sampling 
  high_fixed_count::Int64
  low_fixed_count::Int64
  gen_reached_count::Int64
  predicted_p::Float64
  average_high::Float64
end

function fidelity( T::Int64, N::Int64, q::Float64, s::Float64, ngens::Int64, init_high_freq::Float64, 
    psel_first::Bool, deterministic_mutation::Bool )
  return fidelity_type(T,N,s,q,ngens,init_high_freq,psel_first,deterministic_mutation,0,0,0,0.0,0.0)
end
 
