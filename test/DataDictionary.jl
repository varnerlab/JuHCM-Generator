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
# Function: DataDictionary
# Description: Holds simulation and model parameters as key => value pairs in a Julia Dict()
# Generated on: 2017-03-03T07:00:56.023
#
# Input arguments:
# time_start::Float64 => Simulation start time value (scalar) 
# time_stop::Float64 => Simulation stop time value (scalar) 
# time_step::Float64 => Simulation time step (scalar) 
#
# Output arguments:
# data_dictionary::Dict{AbstractString,Any} => Dictionary holding model and simulation parameters as key => value pairs 
# ----------------------------------------------------------------------------------- #
function DataDictionary(time_start,time_stop,time_step)

	# Load the networks from disk - 
	stoichiometric_matrix = readdlm("Network.dat");
	flux_balance_modes_matrix = readdlm("Modes.dat");

	# How many modes do we have? - 
	(number_of_reactions,number_of_modes) = size(flux_balance_modes_matrix)

	# Setup the initial_condition_array - 
	enzyme_initial_conditions = ones(number_of_modes)
	initial_condition_array = [
		enzyme_initial_conditions	;	# enzymes
		0.0	;	# 1 A_e
		0.0	;	# 2 B_e
		0.0	;	# 3 C_e
	]

	# Setup the index_vector_external_species - 
	index_vector_external_species = [
		4	;	# 1 A_e
		5	;	# 2 B_e
		6	;	# 3 C_e
	]

	# Setup the degradation_constant_array - 
	default_protein_half_life = 24.0	# units:hr
	default_degrdation_constant = -(1/protein_half_life)*log(0.5)	# units: hr^-1
	degradation_constant_array = default_degrdation_constant*ones(number_of_modes)

	# Setup the external stoichiometric matrix - 
	external_stoichiometric_matrix = stoichiometric_matrix[index_vector_external_species,:]

	# =============================== DO NOT EDIT BELOW THIS LINE ============================== #
	data_dictionary = Dict{AbstractString,Any}()
	data_dictionary["stoichiometric_matrix"] = stoichiometric_matrix
	data_dictionary["external_stoichiometric_matrix"] = external_stoichiometric_matrix
	data_dictionary["initial_condition_array"] = initial_condition_array
	data_dictionary["degradation_constant_array"] = degradation_constant_array
	data_dictionary["flux_balance_modes_matrix"] = flux_balance_modes_matrix
	data_dictionary["index_vector_external_species"] = index_vector_external_species
	# =============================== DO NOT EDIT ABOVE THIS LINE ============================== #
	return data_dictionary
end
