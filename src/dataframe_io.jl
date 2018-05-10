using DataFrames
# copied from spatial_df/src/dataframe_utils.jl
function read_params( filename::AbstractString )
  lines = []
  open(filename) do f
    while !eof(f)
      line = readline(f)
      if line[1] == '#'
        Base.push!(lines,line)
      else
        break
      end
    end
  end
  lines
end


# Does not quote strings
function df_row( df::DataFrame, row::Int64 )
  result = []
  for name in names(df)
    Base.push!(result,df[row,name])
  end
  join(result,",")
end

function my_parse( str::AbstractString )
  str = try 
          tmp = str
          parse(Int,tmp) 
        catch 
          str 
        end
  str = try 
          tmp = str
          parse(Float64,tmp) 
        catch 
          str 
        end
  str
end

@doc """ function read_dataframe( in_fname::AbstractString )
  Reads a CSV file into a dataframe.
  Comment lines in the CSV file are preceded by hash marks 
  Comment lines must precede the header line and the data lines.
"""
function read_dataframe( in_fname::AbstractString )
  f = open(in_fname) 
  line = ""
  lnum = 0
  while !eof(f)
    line = readline(f)
    lnum += 1
    if strip(line)[1] != '#'
      break
    end
  end
  headers = map(Symbol,split(line,","))
  ncols = length(headers)
  df = DataFrame()
  for hdr in headers
    df[hdr] = []
  end
  while !eof(f)
    line = readline(f)
    lnum += 1
    vals = split(line,",")
    if length(vals) != ncols
      error("incorrect number of values in line ",lnum," of file: ",in_fname)
    end
    i = 1
    for hdr in headers
      Base.push!( df[hdr], my_parse(vals[i]) )
      i += 1
    end
  end
  df
end
    

function write_dataframe( in_fname::AbstractString, df::DataFrame, out_fname::AbstractString )
  params = read_params( in_fname )
  open(out_fname,"w") do f 
    for param in params
      println(f,param)
    end
    println(f,join(map(string,names(df)),","))
    for i =1:size(df)[1]
      println(f,df_row( df, i ) )
    end
  end
end
