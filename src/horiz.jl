export horiz_transfer_circular!, new_emmigrants_funct, add_emmigrants, horiz_transfer_by_fitness!
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
  if topology == "ring"
    for j = 1:tr.num_subpops
      k_forward = (j+tr.num_subpops-2)%tr.num_subpops+1
      k_backward = j%tr.num_subpops+1
      if means[k_forward] > means[j] 
        k = (j+tr.num_subpops-2)%tr.num_subpops+1
      else
        k = j%tr.num_subpops+1
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
      k = max_index
      #println("j: ",j,"  k: ",k,"  means[j]: ",means[j],"  means[k]: ",means[k])
    end
  end
  #println("neighbor_list: ",neighbor_list)
  new_emmigrants = new_emmigrants_funct( meta_pop, tr, vt, neighbor_list, ideal, id, emmigrant_select=emmigrant_select )
  add_emmigrants( meta_pop, tr, vt, new_emmigrants, neg_select=neg_select)
end

@doc """ new_emmigrants_funct()
  Choose emmigrants from source population.
  Chosen using fitness proportional selection if emmigrant_select==true
  Chosen randonly if emmigrant_select==false
"""
function new_emmigrants_funct( meta_pop::PopList, tr::temporal_result_type, vt::Dict{Int64,variant_type}, neighbor_list::Vector{Int64}, 
    ideal::Vector{Float64}, id::Vector{Int64}; emmigrant_select::Bool=true )
  subpop_size = Int(floor(tr.N/tr.num_subpops))
  emmigrants = PopList()
  for j = 1:tr.num_subpops
    if emmigrant_select
      Base.push!( emmigrants, propsel( meta_pop[j], tr.ne, vt ) )
    else
      s = StatsBase.sample(collect(1:subpop_size),tr.ne,replace=false,ordered=true) # random sample of indices
      Base.push!( emmigrants, meta_pop[j][s] )   # Neutral
    end
  end
  #println("emmigrants: ",emmigrants)
  new_emmigrants = Population[ Population() for j = 1:tr.num_subpops ]
  for j = 1:tr.num_subpops
    k = neighbor_list[j]
    if k != j
      # Create new variants for the emmigrants in the new subpop
      for e in emmigrants[k]   # meta_pop[k] is the source, meta_pop[j] is the destination
        i = id[1]
        #println("e: ",e,"  i: ",i)
        #println("new emmigrant i: ",i,"  subpop_index:",k,"  num_attributes: ",tr.num_attributes )
        vt[i] = deepcopy(vt[e])
        # The following line is not needed if fitness is static
        #vt[i].fitness = fitness( vt[i].attributes, ideal, min_fit=tr.min_fit, linear_fitness=tr.linear_fitness )  
        #println("vt[",e,"]: ",vt[e])
        #println("vt[",i,"]: ",vt[i])
        Base.push!( new_emmigrants[j], i )
        id[1] += 1
      end
    end
  end
  new_emmigrants
end

function add_emmigrants( meta_pop::PopList, tr::temporal_result_type, vt::Dict{Int64,variant_type}, new_emmigrants::PopList;
      neg_select::Bool=true )
  subpop_size = Int(floor(tr.N/tr.num_subpops))
  for j = 1:tr.num_subpops
    if length(new_emmigrants[j]) > 0 
      pop_after_deletion = Population[]
      #println("j: ",j,"  j%tr.num_subpops+1: ",j%tr.num_subpops+1,"  (j+tr.num_subpops-2)%tr.num_subpops+1: ",(j+tr.num_subpops-2)%tr.num_subpops+1)
      if neg_select  # use reverse proportional selection to delete elements by negative fitness
        pop_after_deletion = reverse_propsel(meta_pop[j],tr.ne,vt)
      else  # delete random elements to delete
        s = StatsBase.sample(collect(1:subpop_size),subpop_size-tr.ne,replace=false,ordered=true) # random sample of indices
        pop_after_deletion = meta_pop[j][s]
        #println("length(pop_after_deletion): ",length(pop_after_deletion))
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

