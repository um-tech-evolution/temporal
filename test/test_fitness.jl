# to run:  julia test_fitness.jl
include("../src/fitness.jl")
attr0 = fill(0.0,3)
ideal0 = attr0
attr1 = [0.0, 0.2, 0.4]
ideal5 = fill(0.5,3)
minFit0 = 0.0
minFit1 = 0.25
attribs = [attr0,attr1]
ideals = [ideal0,ideal5]
minFits = [minFit0,minFit1]

function test_fitness()
    for ideal in ideals
      for attrib in attribs
        for minFit in minFits
          fit = fitness(attrib,ideal,minFit=minFit,linfit_slope=1.0)
          println("ideal: ",ideal,"  attr: ",attrib,"  minFit: ",minFit,"  fit: ",fit)
        end
      end
    end
end

test_fitness()
