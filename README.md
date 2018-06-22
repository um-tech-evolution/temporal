Julia code for cultural evolution simulation of a quantitative trait which undergoes temporal 
change.  There is a metapopulation of individuals with variants of the trait.  The metapopulation
is subdivided into subpopulations, and there may or may not be horizontal transfer of traits 
between subpopulations.

In support of a future cultural evolution paper on fidelity, populaiton size, complexity,
and temporal varying environment.

There are two types of simulation:
  simtype == 2:  Equilibrium results on mean fitness, mean attribute standard deviation,
      fraction of subpops at minimum fitness.
  simtype == 1:  Results on the number of generations or move updates until all individuals at at minimum fitness.

src/complexity.jl  is a program that takes the output csv file from a run of simtype==1 and generates a table (for a plot)
  of complexity (number of attributes) as a function of population size (N) and mutation standard deviation (mutStddev).

See src/types.jl for documentation of parameters.
