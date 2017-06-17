#!/usr/local/julia/bin/julia

global s, q
s = q = 0.0

if length(ARGS)>=2
  s = float(ARGS[1])
  q = float(ARGS[2])
else
  println("not enough args")
end

p1 = (s*q+q-1.0)/s
p2 = (s*q+q-1.0)/(q*s)
println("s: ",s,"  q: ",q,"  p1: ",p1,"  p2:",p2)
