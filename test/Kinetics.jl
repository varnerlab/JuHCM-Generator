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
#
# ----------------------------------------------------------------------------------- #
# Function: Kinetics
# Description: Calculate the flux array at time t
# Generated on: 2017-03-03T20:52:45.099
#
# Input arguments:
# t::Float64 => Current time value (scalar) 
# x::Array{Float64,1} => State array (number_of_species x 1) 
# data_dictionary::Dict{AbstractString,Any} => Dictionary holding model parameters 
#
# Output arguments:
# flux_array::Array{Float64,1} => Flux array (number_of_rates x 1) at time t 
# ----------------------------------------------------------------------------------- #
function Kinetics(time,state_array,data_dictionary)

	# Get data/parameters from the data_dictionary - 
	flux_balance_modes_matrix = data_dictionary["flux_balance_modes_matrix"]
	K = data_dictionary["saturation_constant_array"]
	VMax = data_dictionary["rate_constant_array"]

	# initialize the kinetic_rate_array - 
	(number_of_reactions,number_of_modes) = size(flux_balance_modes_matrix)
	kinetic_rate_array = zeros(number_of_modes)

	# alias the state_array - 
	E_M1 = state_array[1]
	E_M2 = state_array[2]
	E_M3 = state_array[3]
	E_M4 = state_array[4]

	A_e = state_array[5]
	B_e = state_array[6]
	C_e = state_array[7]

	# build the kinetic rates - 
	kinetic_rate_array[1] = VMax[1]*E_M1*(A_e/(K[1]+A_e))
	kinetic_rate_array[2] = VMax[2]*E_M2*(C_e/(K[2]+C_e))
	kinetic_rate_array[3] = VMax[3]*E_M3*(A_e/(K[3]+A_e))*(C_e/(K[4]+C_e))
	kinetic_rate_array[4] = VMax[4]*E_M4*(A_e/(K[5]+A_e))

	return kinetic_rate_array
end
