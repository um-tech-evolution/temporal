# Input:  csv file produced by temporal run with simtype==1 which runs the simulation until all subpops
#   lose all individuals with fitness greater than minFit.  
# This run must be done with a list of values for the paramters N, num_attributes, and mutStddev.: (These lists may be of length 1.)
using DataFrames, CSV
include("dataframe_io.jl")   # load functions to read and write dataframes
# cutoff should be set to be slightly less than ngens in the input .csv file.  
global cutoff = 980

# For each pop size N and mutStddev, find the maximum number of attributes where the generational lifetime is at least cutoff.
#    Note that the max possible lifetime is ngens (as set in the parameter file that generated the input CSV file).
# The options for lifetime_field are :generational_lifetime and :move_update_lifetime
function find_max_attr( df::DataFrame, cutoff::Number, lifetime_field::Symbol)
  println("lifetime_field: ",lifetime_field)
  N_list = Int[]
  try
    N_list = unique( df[:,:N] )
  catch
    error("The parameter N must be specified as a list in the configuration file for the temporal run.")
  end
  num_attributes_list = Int[]
  try
    num_attributes_list = unique( df[:,:num_attributes] )
  catch
    error("The parameter num_attributes must be specified as a list in the configuration file for the temporal run.")
  end
  mut_stddev_list = Float64[]
  try
    mut_stddev_list = unique( df[:,:mutStddev] )
  catch
    error("The parameter mutStddev must be specified as a list in the configuration file for the temporal run.")
  end
  #println("num_attributes_list: ",num_attributes_list)
  num_rows = size( df[(df[:N].==N_list[1]) .& (df[:mutStddev].==mut_stddev_list[1]),:])[1]
  result_df = DataFrame()
  result_df[:N] = N_list
  for stddev in mut_stddev_list
    k = 1
    result = fill(num_attributes_list[1],length(N_list))
    for N in N_list
      #print("N: ",N)
      df_N = df[(df[:N].==N) .& (df[:mutStddev].==stddev),:]
      
      #if N == N_list[1] 
      #  println("df_N; ",df_N)
      #end
      
      j = findfirst(x->x<=cutoff,df_N[:,lifetime_field])
      if j > 1
        result[k] = num_attributes_list[j-1] 
      elseif j==1
        result[k] = 0
      else
        result[k] = num_attributes_list[end]
        println("Warning: All attributes have fidelity for N = ",N," and mutStddev = ",stddev) 
        println("Consider running simulation with more attributes.")
      end
      #=
      for i = num_rows:-1:1
        #println("i: ",i,"  df_N[i,lifetime_field]: ",df_N[i,lifetime_field])
        if df_N[i,lifetime_field] > cutoff
          result[k] = num_attributes_list[i]
          println("  result[",k,"]: ",result[k])
          break
        end
      end
      =#
      k += 1
    end
    if lifetime_field == :generational_lifetime
      result_df[:lifetime_field] = :gen
    elseif lifetime_field == :move_update_lifetime
      result_df[:lifetime_field] = :mov
    end
    result_df[Symbol("stddev=$stddev")] = result
  end
  result_df 
end



fname = ARGS[1]
if fname[end-3:end] != ".csv"
  fname = "$(fname).csv"
end
if length(ARGS) >= 2 && ARGS[2][1] == 'g'
    lifetime_field = :generational_lifetime
  elseif length(ARGS) >= 2 && ARGS[2][1] == 'm'
    lifetime_field = :move_update_lifetime
  else  # default is now generational lifetime
    lifetime_field = :generational_lifetime
end
if !(lifetime_field == :generational_lifetime || lifetime_field == :move_update_lifetime)
  error("lifetime_field must be either generational_lifetime or move_update_lifetime")
else
  if lifetime_field == :generational_lifetime
    println("Using generational lifetime")
  end
  if lifetime_field == :move_update_lifetime
    println("Using move update lifetime")
  end
end
df = read_dataframe( fname )
println("read dataframe from file: ",fname," of size: ",size(df))
# delete column that is not used
#delete!(df,:gen_limit_reached_count)  # delete unused column
resultdf = find_max_attr( df, cutoff, lifetime_field )
println(resultdf)
if lifetime_field == :generational_lifetime 
  write_dataframe(fname,resultdf,"$(fname[1:end-4]).cmplx_gen.csv")
elseif lifetime_field == :move_update_lifetime
  write_dataframe(fname,resultdf,"$(fname[1:end-4]).cmplx_mu.csv")
end

