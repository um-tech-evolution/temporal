#using DataStructures
export propsel, propsel!, reverse_propsel

@doc """ function propsel()
Apply proportional selection to Population pop using fitness, 
and return the result.  
"""
function propsel( pop::Population, variant_table::Dict{Int64,temporal_variant_type} )
  new_pop = deepcopy(pop)
  propsel!( new_pop, variant_table )
  new_pop
end

@doc """ function propsel()
Apply proportional selection to Population pop to generate a new population of size n.
Assumes that all elements (integers) in pop have their fitnesses defined in variant_table.
"""
function propsel( pop::Population, n::Int64, variant_table::Dict{Int64,temporal_variant_type} )
  #println("propsel n: ",n,"  pop: ",pop)
  fmax = 0.0
  for v in pop
    fitv = variant_table[v].fitness
    if fitv < 0.0
      error("negative fitness in propsel")
    end
    if fitv > fmax
      fmax = fitv
    end
  end
  if fmax == 0.0
    # all elements have fitness zero
    println("all elements have fitness zero")
    deepcopy(pop)
  end

  new_pop = zeros(Int64,n)
  N = length(pop)
  k = 0
  while k < n
    i = rand(1:N)
    w = variant_table[pop[i]].fitness/fmax
    if rand() < w
      new_pop[k + 1] = pop[i]
      k += 1
    end
  end
  #println("end propsel: new_pop: ",new_pop)
  new_pop
end

@doc """function propsel!()
Conduct proportional selection in-place.  
"""
function propsel!( pop::Population, variant_table::Dict{Int64,temporal_variant_type} )
  fmax = 0.0
  for v in pop
    fitv = variant_table[v].fitness
    if fitv < 0.0
      error("negative fitness in propsel!")
    end
    if fitv > fmax
      fmax = fitv
    end
  end
  if fmax == 0.0
    # all elements have fitness zero
    println("all elements have fitness zero")
    return
  end

  n = length(pop)
  selected = zeros(Int64, n)
  k = 0
  while k < n
    i = rand(1:n)
    w = variant_table[pop[i]].fitness/fmax
    if rand() < w
      selected[k + 1] = i
      k += 1
    end
  end
  pop[:] = [ pop[selected[i]] for i = 1:n ]
end

@doc """function reverse_propsel()
Apply reverse proportional selection to pop to return the population without the selected elements.
Reverse proportional selection uses rescaled netative fitness to favor lower fitness individuals.
The rescaling is  fit -> -fit + fmax + fmin.  
Note that under this rescaling, fmax -> fmin,  fmin -> fmax.
Assumes that all elements (integers) in pop have their fitnesses defined in variant_table.
n  must be less than N/2.
"""
function reverse_propsel( pop::Population, n::Int64, variant_table::Dict{Int64,temporal_variant_type} )
  N = length(pop)
  #println("reverse_propsel N: ",N,"  n: ",n)
  if n > Int(ceil(N/2))
    error("n  (number selected) is greater than N/2 in reverse_propsel!")
  end
  fmax = 0.0
  fmin = typemax(Float64)
  for v in pop
    fitv = variant_table[v].fitness
    if fitv < 0.0
      error("negative fitness in reverse_propsel!")
    end
    if fitv > fmax
      fmax = fitv
    end
    if fitv < fmin
      fmin = fitv
    end
  end
  if fmax == 0.0
    # all elements have fitness zero
    return
  end

  selected = trues(N)
  k = 0
  while k < n
    i = rand(1:N)
    #println("i: ",i," s[i]: ",selected[i],"  fit: ",variant_table[pop[i]].fitness)
    while !selected[i]
      i = rand(1:N)
      #println("  i: ",i," s[i]: ",selected[i])
    end
    w = (-variant_table[pop[i]].fitness+fmax+fmin)/fmax
    #println("w: ",w)
    if rand() < w
      selected[i] = false
      k += 1
      #println("selected")
    end
  end
  pop[selected]
end

@doc """ function propsel_alt()
  Roulette wheel version of proportial selection
  Less efficient in almost all circumstances than the above version.
  Written as check of correctness for the above version
"""
function propsel_alt( pop::Population, dfe::Function)
  fitness_table = Dict{Int64,Float64}()
  N = length(pop)
  new_pop = zeros(Int64,N)
  fitdict = Dict{Int64,Float64}()
  popctr = pop_counter( pop )
  for p in pop
    fit = dfe_fitness(p, dfe, fitness_table )  # set fitness of p to be dfe(p).
  end
  sum_fitness = 0.0
  for p in keys(popctr)
    dfit = dfe_fitness( p, dfe, fitness_table )*popctr[p]
    sum_fitness += dfit
  end
  for p in keys(popctr)
    fitdict[p] = dfe_fitness( p, dfe, fitness_table )*popctr[p]/sum_fitness
  end
  for i = 1:N
    r = rand()
    sumfit = 0.0
    for p in keys(fitdict)
      sumfit += fitdict[p]
      if r < sumfit
        new_pop[i] = p
        break
      end
    end
  end
  new_pop
end

@doc """ function negative_propsel()
Apply proportional selection to Population pop to generate the population of size n to be deleted
Use rescaled netative fitness to favor lower fitness individuals.
The rescaling is  fit -> -fit + fmax + fmin.  
Under this rescaling, fmax -> fmin,  fmin -> fmax.
Assumes that all elements (integers) in pop have their fitnesses defined in variant_table.
Not as useful as reverse_propsel() defined above.
"""
function negative_propsel( pop::Population, n::Int64, variant_table::Dict{Int64,temporal_variant_type} )
  new_pop = zeros(Int64,n)
  fmax = 0.0
  fmin = typemax(Float64)
  for v in pop
    fitv = variant_table[v].fitness
    if fitv < 0.0
      error("negative fitness in negative_propsel")
    end
    if fitv > fmax
      fmax = fitv
    end
    if fitv < fmin
      fmin = fitv
    end
  end
  if fmax == 0.0
    # all elements have fitness zero
    return
  end

  N = length(pop)
  k = 0
  while k < n
    i = rand(1:N)
    w = (-variant_table[pop[i]].fitness+fmax+fmin)/fmax
    if rand() < w
      new_pop[k + 1] = pop[i]
      k += 1
    end
  end
  new_pop
end
