export fidelity_type, fidelity

type fidelity_type 
  T::Int64
  N::Int64
  s::Float64
  q::Float64
  ngens::Int64
  init_high_freq::Float64
  psel_first::Bool
  high_fixed_count::Int64
  low_fixed_count::Int64
  gen_reached_count::Int64
  predicted_high::Float64
  average_high::Float64
end

function fidelity( T::Int64, N::Int64, q::Float64, s::Float64, ngens::Int64, init_high_freq::Float64, psel_first::Bool )
  return fidelity_type(T,N,s,q,ngens,init_high_freq,psel_first,0,0,0,0.0,0.0)
end
 
