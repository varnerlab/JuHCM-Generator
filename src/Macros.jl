macro include_function(src_file_name)

  # create src_buffer -
  src_buffer::Array{AbstractString} = AbstractString[]

  # path to distrubtion -
  path_to_src_file = "./include/"*src_file_name*".jl"
  open(path_to_src_file,"r") do src_file
    for line in eachline(src_file)
      push!(src_buffer,line)
    end
  end

  string_value = ""
  for line in src_buffer
    string_value *= line
  end

  return string_value
end

macro include_function_matlab(src_file_name)

  # create src_buffer -
  src_buffer::Array{AbstractString} = AbstractString[]

  # path to distrubtion -
  if (is_windows() == true)
    path_to_src_file = "$(pwd())\\include\\"*src_file_name*".m"
  else
    path_to_src_file = "./include/"*src_file_name*".m"
  end

  open(path_to_src_file,"r") do src_file
    for line in eachline(src_file)
      push!(src_buffer,line)
    end
  end

  string_value = ""
  for line in src_buffer
    string_value *= line
  end

  return string_value
end
