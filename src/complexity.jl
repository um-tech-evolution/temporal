# Input:  csv file produced by Temporal.jl run with simtype==1 which runs the simulation until all subpops
#   lose all individuals with fitness greater than minFit.  
#  This program uses generational update results rather than move-update results.
#  For now, there is no temporal variation in the environment, and only a single subpopulation (and thus no horizontal transmission).
# For each population
using DataFrames, CSV
include("dataframe_io.jl")

function find_max_attr( df::DataFrame, cutoff::Number  )
  N_list = unique( df[:,:N] )
  num_attributes_list = unique( df[:,:num_attributes] )
  mut_stddev_list = unique( df[:,:mutStddev] )
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
      
      j = findfirst(x->x<=cutoff,df_N[:,:generational_lifetime])
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
        #println("i: ",i,"  df_N[i,:generational_lifetime]: ",df_N[i,:generational_lifetime])
        if df_N[i,:generational_lifetime] > cutoff
          result[k] = num_attributes_list[i]
          println("  result[",k,"]: ",result[k])
          break
        end
      end
      =#
      k += 1
    end
    result_df[Symbol("stddev=$stddev")] = result
  end
  result_df 
end



fname = ARGS[1]
# use julia v. 5 to avoid depwarn messages from readtabltddev = 
#df = readtable(fname, makefactors=true, allowcomments=true)
df = read_dataframe( fname )
println("size: ",size(df))
delete!(df,:move_update_lifetime)
delete!(df,:gen_limit_reached_count)
resultdf = find_max_attr( df, 980 )
println(resultdf)
write_dataframe(fname,resultdf,"$(fname[1:end-4]).cmplx.csv")

