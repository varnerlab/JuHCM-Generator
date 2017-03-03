# ----------------------------------------------------------------------------------- #
# Copyright (c) 2017 Varnerlab
# Robert Frederick Smith School of Chemical and Biomolecular Engineering
# Cornell University, Ithaca NY 14850
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# ----------------------------------------------------------------------------------- #
function Control(time::Float64,state_array::Array{Float64,1},kinetic_rate_array::Array{Float64,1},data_dictionary::Dict{AbstractString,Any})

  # Initialize the control variable array -
  cybernetic_variable_array = Float64[]

  # Get the mode_matrix
  mode_matrix = data_dictionary["flux_balance_modes_matrix"]
  (number_of_reactions,number_of_modes) = size(mode_matrix)

  # Calculate the cybernetic_Variables
  fill!(cybernetic_variable_array,0.0)
  for mode_index = 1:number_of_modes

    # create the 'cybernetic' set -
    push!(cybernetic_variable_array,mode_matrix[1,mode_index]*kinetic_rate_array[mode_index])
  end

  # what is ths maximum rate?
  maximum_rate = maximum(cybernetic_variable_array)

  # Initialize cybernetic variables u and v -
  v = zeros(number_of_modes)
  for mode_index = 1:number_of_modes
    v[mode_index] = cybernetic_variable_array[mode_index]/maximum_rate
  end

  # return the cybernetic v-variable -
  return v
end
