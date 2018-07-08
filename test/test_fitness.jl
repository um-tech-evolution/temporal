# to run:  julia test_fitness.jl
include("../src/fitness.jl")
attr0 = fill(0.4,3)
ideal0 = attr0
attr1 = [0.3, 0.35, 0.4]
ideal5 = fill(0.5,3)
minFit0 = 0.20
minFit1 = 0.25
attribs = [attr0,attr1]
ideals = [ideal0,ideal5]
minFits = [minFit0,minFit1]
fitnesses = [0.5,0.5,
0.3881966011250105,0.3881966011250105,
0.3267949192431123,0.3267949192431123,
0.23074175964327476,0.25]

function test_fitness()
    i = 1
    for ideal in ideals
      for attrib in attribs
        for minFit in minFits
          fit = fitness(attrib,ideal,minFit=minFit,linfit_slope=1.0)
          println("ideal: ",ideal,"  attr: ",attrib,"  minFit: ",minFit,"  fit: ",fit,"  ans: ",fitnesses[i])
          @assert isapprox(fit,fitnesses[i])
          i+=1
        end
      end
    end
end

test_fitness()
