#= horiz.jl
Do horizontal transmission between subpopulations.
The source subpopulation is a neighbor of the destination subpopulation in the given topology.
Possible topologies are:  circular, ring, global, vonneuman, moore.
There are 3 ways in which horizontal selection can be done:
  1.  If probHSelect > 0, with probability probHSelect, the source population is chosen to have maximum mean fitness over the neighborhood.
  2.  If emigrant_select == true, then emigrants from the source population are chosen by proportional selection.
  3.  If neg_select == true, then the individuals to be replaced from the destination population are chosen by reverse proportional selection.
Currently (5/2018), if the parameter horiz_select == true, both emigrant_select and neg_select are set to true.
If linear fitness is used with a large value of minFit (such as 0.4), then proportional selection and reverse prop. sel. are weak.
Thus, probHSelect > 0 give much stronger horizontal selection.
=#

using Base.Test
export horiz_transfer, horiz_transfer_circular!, new_emigrants_funct, add_emigrants, horiz_transfer_by_fitness!,
  new_emigrants_funct

function horiz_transfer( meta_pop::PopList, tp::param_type, vt::Dict{Int64,temporal_variant_type}, ideal::Vector{Float64}, mmeans::Vector{Float64},
      id::Vector{Int64}, generation::Int64  )
  emigration = tp[:topology] != "none" && tp[:num_subpops] != 1 && (tp[:num_emigrants] > 0 || (tp[:migration_rate] != :null && tp[:migration_rate] > 0.0))
  #println("emigration: ",emigration)
  # Set num_emigrants from migration_rate
  if ! emigration
    return
  end
  # Check that tp[:num_emigrants] and tp[:migration_rate] are not both nonzero.
  if tp[:migration_rate] != :null && tp[:migration_rate] != 0.0 && tp[:num_emigrants] > 0 
    error("parameters num_emigrants and migration_rate cannot not both be nonzero.")
  end
  if tp[:topology]=="circular"   # Circular is not "by fitness"
    source_subpop_list = horiz_transfer_circular!( meta_pop, tp, vt, ideal, id, generation, neg_select=tp[:horiz_select], emigrant_select=tp[:horiz_select] )
  elseif  tp[:topology]=="ring" || tp[:topology]=="global"
    source_subpop_list = horiz_transfer_by_fitness!( meta_pop, tp, vt, ideal, mmeans, id, neg_select=tp[:horiz_select], topology=tp[:topology], emigrant_select=tp[:horiz_select] )
  elseif  tp[:num_subpops] >= 9 && (tp[:topology]=="vonneumann" || tp[:topology]=="moore")
    source_subpop_list = horiz_transfer_by_fitness!( meta_pop, tp, vt, ideal, mmeans, id, neg_select=tp[:horiz_select], topology=tp[:topology], emigrant_select=tp[:horiz_select] )
  elseif tp[:num_subpops] > 1 && tp[:topology]!="none" 
    println("nsbp: ",tp[:num_subpops],"  topo: ",tp[:topology],"  ne: ",tp[:num_emigrants],"  test: ",(tp[:num_subpops] >= 9 && tp[:num_emigrants] > 0 && (tp[:topology]=="vonneumann" || tp[:topology]=="moore")))
    println("Warning! no horizontal transfer done with tp[:num_subpops]=",tp[:num_subpops]," and topology=",tp[:topology])
  end
  #println("source subpop list: ",source_subpop_list)
  #print_meta_pop_attributes( tp, meta_pop, vt )
  new_emigrants = new_emigrants_funct( meta_pop, tp, vt, source_subpop_list, ideal, id, emigrant_select=tp[:horiz_select] )
  #print_meta_pop_attributes( tp, meta_pop, vt )
  #println("new_emigrants: ",new_emigrants)
  add_emigrants( meta_pop, tp, vt, new_emigrants, neg_select=false)
end

#=
# Check topology and horizontal transfer settings for correctness
function horiz_param_check( topology_list::Vector{String}, num_subpops_list::Vector{Int64}, num_emigrants_list::Vector{Int64} )
  println("topology_list: ",topology_list)
  println("num_subpops_list: ",num_subpops_list)
  println("num_emigrants_list: ",num_emigrants_list)
  for topology in topology_list
    for num_subpops in num_subpops_list
      for ne in num_emigrants_list
        print("t: ",topology,"  nsbp: ",num_subpops,"  ne: ",ne,"  ")
        if topology == "none" || num_subpops == 1
          println("OK")
        elseif ne > 0 && topology=="circular"
          println("OK")
        elseif  ne > 0 && topology=="ring" || topology=="global"
          println("OK")
        elseif  num_subpops >= 9 && ne > 0 && topology=="vonneumann" || topology!="moore"
          println("OK")
        elseif num_subpops > 1 && topology!="none" && ne > 0
          error("Warning!!! no horizontal transfer done with num_subpops=",num_subpops," and topology=",topology)
        end
      end
    end
  end
end
=#

@doc """ horiz_transfer_circular!()
  Transfers variants between subpopulations in a circular fashion (either forward or backward).
  Elements to be transfered are selected by proportional selection.
  Elements to be replaced can be random or selected by reverse proportional selection depending on the flag neg_select.
  metapop is modified by this function (as a side effect)
"""
function horiz_transfer_circular!(  meta_pop::PopList, tp::param_type, vt::Dict{Int64,temporal_variant_type}, ideal::Vector{Float64}, id::Vector{Int64},
     generation::Int64; neg_select::Bool=true, emigrant_select::Bool=true )
  #println("horiz_transfer_circular! forward: ",forward,"  num_attributes: ",tp[:num_attributes])
  forward = generation % 2 == 0 ? true : false  
  subpop_size = Int(floor(tp[:N]/tp[:num_subpops]))
  source_subpop_list = zeros(Int64,tp[:num_subpops])
  for j = 1:tp[:num_subpops]
    #println("j: ",j,"  j%tp[:num_subpops]+1: ",j%tp[:num_subpops]+1,"  (j+tp[:num_subpops]-2)%tp[:num_subpops]+1: ",(j+tp[:num_subpops]-2)%tp[:num_subpops]+1)
    if forward
      k = (j+tp[:num_subpops]-2)%tp[:num_subpops]+1
    else
      k = j%tp[:num_subpops]+1
    end
    source_subpop_list[j] = k
  end
  return source_subpop_list
end

@doc """ function horiz_transfer_by_fitness!( )
  Do horizontal transfer where, with probability tp[:probHSelect], the source subpopulation has lower fitness than the destination.
  The source population will be a neighbor of the destination population in the topology.
  Options for topology are:
    "ring":    neighbors are determined by subpop index with wraparound (as in horis_transfer_circular)
    "moore":   each subpop has grid coordinates, the neighborhood is a Moore neighborhood (with 8 neighbors) in this grid
    "vonneumann":   each subpop has grid coordinates, the neighborhood is a von Neumann neighborhood (with 4 neighbors) in this grid
    "global":    each subpop is a neighbor of all other subpops
  For grid options, the metapop size (tp[:N]) should be a multiple of the number of subpops (tp[:num_subpops]).
"""    
function horiz_transfer_by_fitness!(  meta_pop::PopList, tp::param_type, vt::Dict{Int64,temporal_variant_type}, ideal::Vector{Float64}, means::Vector{Float64},
      id::Vector{Int64}; topology::String="ring", neg_select::Bool=true, emigrant_select::Bool=true )
  #println("horiz_transfer_by_fitness!  topology: ",tp[:topology],"  means: ",means,"  tp[:probHSelect]: ",tp[:probHSelect]) 
  source_subpop_list = zeros(Int64,tp[:num_subpops])
  # Note:  source is the source population, dest is the destination population
  if topology == "ring"
    for dest = 1:tp[:num_subpops]
      source_forward = (dest+tp[:num_subpops]-2)%tp[:num_subpops]+1
      source_backward = dest%tp[:num_subpops]+1
      if means[source_forward] > means[dest] 
        #source = (dest+tp[:num_subpops]-2)%tp[:num_subpops]+1
        source = source_forward
      else
        #source = dest%tp[:num_subpops]+1
        source = source_backward
      end
      source_subpop_list[dest] = source
      #println("dest: ",dest,"  source: ",source,"  means[dest]: ",means[dest],"  means[source]: ",means[source])
    end
  else 
    ncols,nrows = factorize( tp[:num_subpops] )
    #ncols = Int(floor(sqrt(tp[:num_subpops])))
    #nrows = Int(floor(tp[:num_subpops]/ncols))
    #println("horiz_transfer_by_fitness!  num_subpops: ",tp[:num_subpops],"  ncols: ",ncols,"  nrows: ",nrows) 
    #println("horiz_transfer_by_fitness!  means: ",means)
    for dest = 1:tp[:num_subpops]
      #println("hzf dest: ",dest)
      if topology == "vonneumann" 
        if ncols < 3 || nrows < 3
          println("ncols: ",ncols,"  nrows: ",nrows)
          error("ncols and nrows must be at least 3 in function horiz_transfer_by_fitness() with topology=='vonneumann'")
        end
        nbd = von_neumann_nbd(ncols,nrows,dest)
      elseif topology == "moore" 
        nbd = moore_nbd(ncols,nrows,dest)
        if ncols < 3 || nrows < 3
          println("ncols: ",ncols,"  nrows: ",nrows)
          error("ncols and nrows must be at least 3 in function horiz_transfer_by_fitness() with topology=='moore'")
        end
      elseif topology == "global" 
        nbd = append!(collect(1:(dest-1)),collect((dest+1):tp[:num_subpops]))
      else
        println("topology: ",topology)
        error("topology in horiz_transfer_by_fitness! must be one of 'ring', 'moore', 'vonneumann',or 'global'.")
      end
      #println("horiz:  ncols: ",ncols,"  nrows: ",nrows,"  nbd:",nbd,"  means[nbd]: ",[means[ii] for ii in nbd])
      #println("horiz: dest: ",dest,"  nbd:",nbd,"  means[nbd]: ",[means[ii] for ii in nbd])
      rnd = rand()
      if rnd < tp[:probHSelect]   # source  is the index of the source population
        max_index = nbd[indmax(means[nbd])]
        if means[max_index] >= means[dest]
          source_subpop_list[dest] = max_index
        else
          source_subpop_list[dest] = dest
        end
        source = source_subpop_list[dest]
        @assert means[source] >= means[dest]    # The source subpop has greater mean fitness than the dest subpop
      else   # choose a random subpop
        source = rand(1:length(nbd))  # source  is the index of a random subpop from the neighborhood
        source_subpop_list[dest] = nbd[source]
      end
      #println("horiz:  dest: ",dest,"  source: ",source,"  means[dest]: ",means[dest],"  means[source]: ",means[source])
    end
  end
  #println("source_subpop_list: ",source_subpop_list)
  return source_subpop_list
end

@doc """ new_emigrants_funct()
  Choose emigrants from source population 
  Chosen using fitness proportional selection if emigrant_select==true
  Chosen randonly if emigrant_select==false
  Emigrants are mutated if tp[:horiz_mutate] is true
  Returns new_emigrants which is a list of the emigrants for each subpopulation
"""
function new_emigrants_funct( meta_pop::PopList, tp::param_type, vt::Dict{Int64,temporal_variant_type}, source_subpop_list::Vector{Int64}, 
    ideal::Vector{Float64}, id::Vector{Int64}; emigrant_select::Bool=true )
  subpop_size = Int(floor(tp[:N]/tp[:num_subpops]))
  new_emigrants = Population[ Population() for j = 1:tp[:num_subpops] ]  # new_emigrants[j] is the list of immigrants into subpop j
  for j = 1:tp[:num_subpops]
    if tp[:num_emigrants] > 0
      num_emigrants = tp[:num_emigrants]
    elseif tp[:migration_rate] > 0  # Set num_emigrants from tp[:migration_rate] and subpop_size
      d = Distributions.Poisson( tp[:migration_rate]*subpop_size )   # d is a Poisson distribution with mean tp[:migration_rate]*subpop_size]
      num_emigrants = min(subpop_size,rand(d))
    else
      error("one of tp[:num_emigrants] or tp[:migration_rate] should be positive.")
    end
    #println("new_emigrants_funct: num_emigrants: ",num_emigrants)
    emigrants = Vector{Int64}()  # emigrants is the list of individuals that will emmigrate from subpop j
    if emigrant_select
      #Base.push!( emigrants, propsel( meta_pop[neighborlist[j]], tp[:num_emigrants], vt ) )
      emigrants = propsel( meta_pop[source_subpop_list[j]], num_emigrants, vt ) 
    else
      s = StatsBase.sample(collect(1:subpop_size),num_emigrants,replace=false,ordered=true) # random sample of indices
      emigrants = meta_pop[source_subpop_list[j]][s]    # Neutral
    end
    #println("j: ",j,"  emigrants: ",emigrants)
    k = source_subpop_list[j]
    if k != j   # only create new_emigrants if the source is different from the destination
      # Create new variants for the emigrants in the new subpop
      for e in emigrants   # meta_pop[k] is the source, meta_pop[j] is the destination
        i = id[1]
        #println("e: ",i,"  vt[e]: ",vt[e])
        vt[i] = deepcopy(vt[e])
        if tp[:horiz_mutate]
          mutate_attributes!( vt[i].attributes, tp[:mutStddev], tp[:additive_error] )   
        end
        fit = fitness( vt[i].attributes, ideal, minFit=tp[:minFit] )  
        #println("new emigrant e: ",e,"  i: ",i,"  fitness: ",fit,"  from subpop: ",k)
        vt[i].fitness = fitness( vt[i].attributes, ideal, minFit=tp[:minFit]  )  
        #println("vt[",e,"]: ",vt[e])
        #println("vt[",i,"]: ",vt[i])
        Base.push!( new_emigrants[j], i )
        id[1] += 1
      end
    end
  end
  #println("new_emigrants: ",new_emigrants)
  new_emigrants
end

@doc """ add_emigrants()
For each j, removes length(new_emigrants[j]) individuals from subpop[j] and replaces them with new_emigrants[j].
If neg_select==true, the the individuals to be removed are chosen by reverse proportional selection, otherwise randomly.
"""
function add_emigrants( meta_pop::PopList, tp::param_type, vt::Dict{Int64,temporal_variant_type}, new_emigrants::PopList;
      neg_select::Bool=true )
  subpop_size = Int(floor(tp[:N]/tp[:num_subpops]))
  #println("B add_emigrants: meta_pop: ",meta_pop)
  for j = 1:tp[:num_subpops]
    #println("add emigrants subpop j: ",j,"  new_emigrants[j]: ",new_emigrants[j])
    if length(new_emigrants[j]) > 0 
      pop_after_deletion = Population[]
      #println("j: ",j,"  j%tp[:num_subpops]+1: ",j%tp[:num_subpops]+1,"  (j+tp[:num_subpops]-2)%tp[:num_subpops]+1: ",(j+tp[:num_subpops]-2)%tp[:num_subpops]+1)
      if neg_select  # use reverse proportional selection to delete elements by negative fitness
        pop_after_deletion = reverse_propsel(meta_pop[j],length(new_emigrants[j]),vt)
      else  # delete random elements to delete
        s = StatsBase.sample(collect(1:subpop_size),subpop_size-length(new_emigrants[j]),replace=false,ordered=true) # random sample of indices
        pop_after_deletion = meta_pop[j][s]
        #println("length(pop_after_deletion): ",length(pop_after_deletion))
        #println("j: ",j,"  pop_after_deletion: ",pop_after_deletion)
      end
      meta_pop[j] = append!( pop_after_deletion, new_emigrants[j] )
    end
  end
  #println("A add_emigrants: meta_pop: ",meta_pop)
  return meta_pop
end

# Computes a factorization of n:   n == ncols*nrows   with ncols <= nrows.
# It may be that ncols < 3, in which case horiz_transfer_by_fitness will fail later.
function factorize( n::Int64 )
  ncols = Int(floor(sqrt(n)))
  nrows = Int(floor(n/ncols))
  while ncols > 1  && ncols*nrows != n
    ncols -= 1
    nrows = Int(floor(n/ncols))
  end
  return ncols,nrows
end


# convert 1-based linear coordinates to zero-based grid coordinates
row_col(ncols,i) = ( Int(floor((i-1)/ncols)), (i-1) %ncols )

# convert zero-based grid coordinates to 1-based linear coordinates
indx(ncols,z) = z[1]*ncols + z[2] + 1

add_1_col(ncols,z) = z[2]+1 >= ncols ? (z[1],0) : (z[1], z[2]+1)
sub_1_col(ncols,z) = z[2] == 0 ? (z[1],ncols-1) : (z[1], z[2]-1)
add_1_row(nrows,z) = z[1]+1 >= nrows ? (0,z[2]) : (z[1]+1, z[2])
sub_1_row(nrows,z) = z[1] == 0 ? (nrows-1,z[2]) : (z[1]-1, z[2])

function von_neumann_nbd(nrows,ncols,i)
  nbd = Int64[]
  rc = row_col(ncols,i)
  #println("rc: ",rc)
  #println("a1c: ",add_1_col(ncols, rc ))
  #println("s1c: ",sub_1_col(ncols, rc ))
  #println("a1r: ",add_1_row(ncols, rc ))
  #println("s1r: ",sub_1_row(ncols, rc ))
  Base.push!( nbd, indx( ncols, add_1_col(ncols, rc ) ))
  Base.push!( nbd, indx( ncols, sub_1_col(ncols, rc ) ))
  Base.push!( nbd, indx( nrows, add_1_row(nrows, rc ) ))
  Base.push!( nbd, indx( nrows, sub_1_row(nrows, rc ) ))
  return nbd
end

function moore_nbd(nrows,ncols,i)
  nbd = Int64[]
  rc = row_col(ncols,i)
  #println("rc: ",rc)
  #println("a1c: ",add_1_col(ncols, rc ))
  #println("a1c: ",add_1_row(nrows,add_1_col(ncols, rc )))
  #println("a1c: ",sub_1_row(nrows,add_1_col(ncols, rc )))
  #println("s1c: ",sub_1_col(ncols, rc ))
  #println("s1c: ",add_1_row(nrows,sub_1_col(ncols, rc )))
  #println("s1c: ",sub_1_row(nrows,sub_1_col(ncols, rc )))
  #println("a1r: ",add_1_row(ncols, rc ))
  #println("s1r: ",sub_1_row(ncols, rc ))
  Base.push!( nbd, indx( ncols, add_1_col(ncols, rc ) ))
  Base.push!( nbd, indx( ncols, add_1_row(nrows,add_1_col(ncols, rc ) )))
  Base.push!( nbd, indx( ncols, sub_1_row(nrows,add_1_col(ncols, rc ) )))
  Base.push!( nbd, indx( ncols, sub_1_col(ncols, rc ) ))
  Base.push!( nbd, indx( ncols, add_1_row(nrows,sub_1_col(ncols, rc ) )))
  Base.push!( nbd, indx( ncols, sub_1_row(nrows,sub_1_col(ncols, rc ) )))
  Base.push!( nbd, indx( nrows, add_1_row(nrows, rc ) ))
  Base.push!( nbd, indx( nrows, sub_1_row(nrows, rc ) ))
  return nbd
end

