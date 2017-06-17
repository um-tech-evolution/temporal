export fd_init

function fd_init()
  global T = 10
  global N = 20
  global s = 0.25
  global q = 0.9
  global ngens = 40
  global init_high_freq=0.5
  global psel_first=false
  global fd = fidelity(T,N,q,s,ngens,init_high_freq,psel_first)
  return fd
end
