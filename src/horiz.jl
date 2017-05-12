export horiz_transfer, horiz_transfer_circular!, new_emmigrants_funct, add_emmigrants, horiz_transfer_by_fitness!

function horiz_transfer( meta_pop::PopList, tr::temporal_result_type, vt::Dict{Int64,variant_type}, ideal::Vector{Float64}, mmeans::Vector{Float64},
      id::Vector{Int64}, generation::Int64  )
  if tr.topology == "none" || tr.num_subpops == 1 || tr.ne == 0
    return
  elseif tr.ne > 0 && tr.topology=="circular"
    horiz_transfer_circular!( meta_pop, tr, vt, ideal, id, generation, neg_select=tr.horiz_select, emmigrant_select=tr.horiz_select )
  elseif  tr.ne > 0 && tr.topology=="ring" || tr.topology=="global"
    horiz_transfer_by_fitness!( meta_pop, tr, vt, ideal, mmeans, id, neg_select=tr.horiz_select, topology=tr.topology, emmigrant_select=tr.horiz_select )
  elseif  tr.num_subpops >= 9 && tr.ne > 0 && (tr.topology=="vonneumann" || tr.topology=="moore")
    horiz_transfer_by_fitness!( meta_pop, tr, vt, ideal, mmeans, id, neg_select=tr.horiz_select, topology=tr.topology, emmigrant_select=tr.horiz_select )
  elseif tr.num_subpops > 1 && tr.topology!="none" && tr.ne > 0
    println("nsbp: ",tr.num_subpops,"  topo: ",tr.topology,"  ne: ",tr.ne,"  test: ",(tr.num_subpops >= 9 && tr.ne > 0 && (tr.topology=="vonneumann" || tr.topology=="moore")))
    println("Warning! no horizontal transfer done with tr.num_subpops=",tr.num_subpops," and topology=",tr.topology)
  end
end

#=
# Check topology and horizontal transfer settings for correctness
function horiz_param_check( topology_list::Vector{String}, num_subpops_list::Vector{Int64}, num_emmigrants_list::Vector{Int64} )
  println("topology_list: ",topology_list)
  println("num_subpops_list: ",num_subpops_list)
  println("num_emmigrants_list: ",num_emmigrants_list)
  for topology in topology_list
    for num_subpops in num_subpops_list
      for ne in num_emmigrants_list
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
function horiz_transfer_circular!(  meta_pop::PopList, tr::temporal_result_type, vt::Dict{Int64,variant_type}, ideal::Vector{Float64}, id::Vector{Int64},
     generation::Int64; neg_select::Bool=true, emmigrant_select::Bool=true )
  #println("horiz_transfer_circular! forward: ",forward,"  num_attributes: ",tr.num_attributes)
  forward = generation % 2 == 0 ? true : false  
  subpop_size = Int(floor(tr.N/tr.num_subpops))
  neighbor_list = zeros(Int64,tr.num_subpops)
  for j = 1:tr.num_subpops
    #println("j: ",j,"  j%tr.num_subpops+1: ",j%tr.num_subpops+1,"  (j+tr.num_subpops-2)%tr.num_subpops+1: ",(j+tr.num_subpops-2)%tr.num_subpops+1)
    if forward
      k = (j+tr.num_subpops-2)%tr.num_subpops+1
    else
      k = j%tr.num_subpops+1
    end
    neighbor_list[j] = k
  end
  new_emmigrants = new_emmigrants_funct( meta_pop, tr, vt, neighbor_list, ideal, id, emmigrant_select=emmigrant_select )
  #println("new emmigrants: ",new_emmigrants)
  add_emmigrants( meta_pop, tr, vt, new_emmigrants, neg_select=neg_select)
end

@doc """ function horiz_transfer_by_fitness!( )
  Do horizontal transfer where the source subpopulation has lower fitness than the destination.
  The source population will be a neighbor of the destination population in the topology.
  Options for topology are:
    "ring":    neighbors are determined by subpop index with wraparound (as in horis_transfer_circular)
    "moore":   each subpop has grid coordinates, the neighborhood is a Moore neighborhood (with 8 neighbors) in this grid
    "vonneumann":   each subpop has grid coordinates, the neighborhood is a von Neumann neighborhood (with 4 neighbors) in this grid
    "global":    each subpop is a neighbor of all other subpops
  For grid options, the metapop size (tr.N) should be a multiple of the number of subpops (tr.num_subpops).
"""    
function horiz_transfer_by_fitness!(  meta_pop::PopList, tr::temporal_result_type, vt::Dict{Int64,variant_type}, ideal::Vector{Float64}, means::Vector{Float64},
      id::Vector{Int64}; topology::String="ring", neg_select::Bool=true, emmigrant_select::Bool=true )
  #println("horiz_transfer_by_fitness!  topology: ",tr.topology,"  means: ",means) 
  neighbor_list = zeros(Int64,tr.num_subpops)
  # Note:  k is the source population, j is the destination population
  if topology == "ring"
    for j = 1:tr.num_subpops
      k_forward = (j+tr.num_subpops-2)%tr.num_subpops+1
      k_backward = j%tr.num_subpops+1
      if means[k_forward] > means[j] 
        #k = (j+tr.num_subpops-2)%tr.num_subpops+1
        k = k_forward
      else
        #k = j%tr.num_subpops+1
        k = k_backward
      end
      neighbor_list[j] = k
      #println("j: ",j,"  k: ",k,"  means[j]: ",means[j],"  means[k]: ",means[k])
    end
  else 
    ncols,nrows = factorize( tr.num_subpops )
    #ncols = Int(floor(sqrt(tr.num_subpops)))
    #nrows = Int(floor(tr.num_subpops/ncols))
    #println("horiz_transfer_by_fitness!  num_subpops: ",tr.num_subpops,"  ncols: ",ncols,"  nrows: ",nrows) 
    for j = 1:tr.num_subpops
      #println("hzf j: ",j)
      if topology == "vonneumann" 
        if ncols < 3 || nrows < 3
          println("ncols: ",ncols,"  nrows: ",nrows)
          error("ncols and nrows must be at least 3 in function horiz_transfer_by_fitness() with topology=='vonneumann'")
        end
        nbd = von_neumann_nbd(ncols,nrows,j)
      elseif topology == "moore" 
        nbd = moore_nbd(ncols,nrows,j)
        if ncols < 3 || nrows < 3
          println("ncols: ",ncols,"  nrows: ",nrows)
          error("ncols and nrows must be at least 3 in function horiz_transfer_by_fitness() with topology=='moore'")
        end
      elseif topology == "global" 
        nbd = append!(collect(1:(j-1)),collect((j+1):tr.num_subpops))
      else
        println("topology: ",topology)
        error("topology in horiz_transfer_by_fitness! must be one of 'ring', 'moore', 'vonneumann',or 'global'.")
      end
      #println("ncols: ",ncols,"  nrows: ",nrows,"  nbd:",nbd,"  means[nbd]: ",[means[ii] for ii in nbd])
      #max_mean = means[nbd[1]]
      max_index = nbd[1]
      for i = 1:length(nbd)
        if means[nbd[i]] > means[max_index]
          max_index = nbd[i]
        end
      end
      if means[max_index] >= means[j]
        neighbor_list[j] = max_index
      else
        neighbor_list[j] = j
      end
      if rand() < tr.probHSelect
        k = max_index   # k is the index of the source population
      else
        k = rand(length(nbd))
      end
      #println("j: ",j,"  k: ",k,"  means[j]: ",means[j],"  means[k]: ",means[k])
    end
  end
  #println("neighbor_list: ",neighbor_list)
  new_emmigrants = new_emmigrants_funct( meta_pop, tr, vt, neighbor_list, ideal, id, emmigrant_select=emmigrant_select )
  add_emmigrants( meta_pop, tr, vt, new_emmigrants, neg_select=neg_select)
end

@doc """ new_emmigrants_funct()
  Choose emmigrants from source population 
  Chosen using fitness proportional selection if emmigrant_select==true
  Chosen randonly if emmigrant_select==false
"""
function new_emmigrants_funct( meta_pop::PopList, tr::temporal_result_type, vt::Dict{Int64,variant_type}, neighbor_list::Vector{Int64}, 
    ideal::Vector{Float64}, id::Vector{Int64}; emmigrant_select::Bool=true )
  subpop_size = Int(floor(tr.N/tr.num_subpops))
  emmigrants = PopList()  # emmigrants is the list of individuals that will emmigrate from subpop k
  for k = 1:tr.num_subpops
    if emmigrant_select
      Base.push!( emmigrants, propsel( meta_pop[k], tr.ne, vt ) )
    else
      s = StatsBase.sample(collect(1:subpop_size),tr.ne,replace=false,ordered=true) # random sample of indices
      Base.push!( emmigrants, meta_pop[k][s] )   # Neutral
    end
  end
  #println("emmigrants: ",emmigrants)
  new_emmigrants = Population[ Population() for j = 1:tr.num_subpops ]  # new_emmigrants[j] is the list of immigrants into subpop j
  for j = 1:tr.num_subpops
    k = neighbor_list[j]
    if k != j   # only create new_emmigrants if the source is different from the destination
      # Create new variants for the emmigrants in the new subpop
      for e in emmigrants[k]   # meta_pop[k] is the source, meta_pop[j] is the destination
        i = id[1]
        #println("e: ",e,"  i: ",i)
        #println("new emmigrant i: ",i,"  subpop_index:",k,"  num_attributes: ",tr.num_attributes )
        vt[i] = deepcopy(vt[e])
        # The following line is not needed if fitness is static
        # Not needed in the current program because fitness is calculated during mutation, and fitness shift does not happen between mutation and horiz transfer
        #vt[i].fitness = fitness( vt[i].attributes, ideal, minFit=tr.minFit, linear_fitness=tr.linear_fitness )  
        #println("vt[",e,"]: ",vt[e])
        #println("vt[",i,"]: ",vt[i])
        Base.push!( new_emmigrants[j], i )
        id[1] += 1
      end
    end
  end
  #println("new_emmigrants: ",new_emmigrants)
  new_emmigrants
end

@doc """ add_emmitrants()
For each j, removes length(new_emmigrants[j]) individuals from subpop[j] and replaces them with new_emmigrants[j].
If neg_select==true, the the individuals to be removed are chosen by reverse proportional selection, otherwise randomly.
"""
function add_emmigrants( meta_pop::PopList, tr::temporal_result_type, vt::Dict{Int64,variant_type}, new_emmigrants::PopList;
      neg_select::Bool=true )
  subpop_size = Int(floor(tr.N/tr.num_subpops))
  for j = 1:tr.num_subpops
    #println("add emmigrants j: ",j)
    if length(new_emmigrants[j]) > 0 
      pop_after_deletion = Population[]
      #println("j: ",j,"  j%tr.num_subpops+1: ",j%tr.num_subpops+1,"  (j+tr.num_subpops-2)%tr.num_subpops+1: ",(j+tr.num_subpops-2)%tr.num_subpops+1)
      if neg_select  # use reverse proportional selection to delete elements by negative fitness
        pop_after_deletion = reverse_propsel(meta_pop[j],tr.ne,vt)
      else  # delete random elements to delete
        s = StatsBase.sample(collect(1:subpop_size),subpop_size-tr.ne,replace=false,ordered=true) # random sample of indices
        pop_after_deletion = meta_pop[j][s]
        #println("length(pop_after_deletion): ",length(pop_after_deletion))
        #println("pop_after_deletion: ",pop_after_deletion)
      end
      meta_pop[j] = append!( pop_after_deletion, new_emmigrants[j] )
    end
  end
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

