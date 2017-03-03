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
function Balances(time::Float64,state_array::Array{Float64,1},dxdt_vector::Array{Float64,1},data_dictionary::Dict{AbstractString,Any})

  # Correct for smalls -
  idx = find(state_array.<1e-9)
  state_array[idx] = 1e-9

  # get parameters fron the data dictionary -
  beta = data_dictionary["degradation_constant_array"]
  mode_matrix = data_dictionary["flux_balance_modes_matrix"]
  stoichiometric_matrix = data_dictionary["external_stoichiometric_matrix"]

  # what is my system size?
  (number_of_reactions,number_of_modes) = size(mode_matrix)

  # Define rate vector
  kinetic_rate_array = Kinetics(time,state_array,data_dictionary);

  # calculate the control vector -
  cybernetic_variable_array = Control(time,state_array,kinetic_rate_array,data_dictionary)

  # correct the rates (rate*control) -
  kinetic_rate_array = kinetic_rate_array.*cybernetic_variable_array;

  # calculate the enzyme balances -
  dxdt_vector[1:number_of_modes] = -(beta).*state_array[1:number_of_modes]

  # calculate the metabolic models -
  dxdt_vector[(number_of_modes+1):end] = stoichiometric_matrix*mode_matrix*kinetic_rate_array
end
