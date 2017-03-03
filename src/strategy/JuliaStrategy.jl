function build_function_header_buffer(comment_dictionary)

  # initialize -
  buffer = ""

  # get some data from the comment_dictionary -
  function_name = comment_dictionary["function_name"]
  function_description = comment_dictionary["function_description"]
  input_arg_array = comment_dictionary["input_args"]
  output_arg_array = comment_dictionary["output_args"]

  buffer *= "# ----------------------------------------------------------------------------------- #\n"
  buffer *= "# Function: $(function_name)\n"
  buffer *= "# Description: $(function_description)\n"
  buffer *= "# Generated on: $(now())\n"
  buffer *= "#\n"
  buffer *= "# Input arguments:\n"

  for argument_dictionary in input_arg_array

    arg_symbol = argument_dictionary["symbol"]
    arg_description = argument_dictionary["description"]

    # write the buffer -
    buffer *= "# $(arg_symbol) => $(arg_description) \n"
  end

  buffer *= "#\n"
  buffer *= "# Output arguments:\n"
  for argument_dictionary in output_arg_array

    arg_symbol = argument_dictionary["symbol"]
    arg_description = argument_dictionary["description"]

    # write the buffer -
    buffer *= "# $(arg_symbol) => $(arg_description) \n"
  end
  buffer *= "# ----------------------------------------------------------------------------------- #\n"

  # return the buffer -
  return buffer
end


function build_copyright_header_buffer(problem_object::ProblemObject)

  # What is the current year?
  current_year = string(Dates.year(now()))

  buffer = ""
  buffer*= "# ----------------------------------------------------------------------------------- #\n"
  buffer*= "# Copyright (c) $(current_year) Varnerlab\n"
  buffer*= "# Robert Frederick Smith School of Chemical and Biomolecular Engineering\n"
  buffer*= "# Cornell University, Ithaca NY 14850\n"
  buffer*= "#\n"
  buffer*= "# Permission is hereby granted, free of charge, to any person obtaining a copy\n"
  buffer*= "# of this software and associated documentation files (the \"Software\"), to deal\n"
  buffer*= "# in the Software without restriction, including without limitation the rights\n"
  buffer*= "# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell\n"
  buffer*= "# copies of the Software, and to permit persons to whom the Software is\n"
  buffer*= "# furnished to do so, subject to the following conditions:\n"
  buffer*= "#\n"
  buffer*= "# The above copyright notice and this permission notice shall be included in\n"
  buffer*= "# all copies or substantial portions of the Software.\n"
  buffer*= "#\n"
  buffer*= "# THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR\n"
  buffer*= "# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,\n"
  buffer*= "# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE\n"
  buffer*= "# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER\n"
  buffer*= "# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,\n"
  buffer*= "# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN\n"
  buffer*= "# THE SOFTWARE.\n"
  buffer*= "# ----------------------------------------------------------------------------------- #\n"

  # return -
  return buffer
end

function build_debug_buffer(problem_object::ProblemObject)
  filename = "Debug.txt"

  buffer = ""

  # write the list of reactions -
  list_of_reactions::Array{ReactionObject} = problem_object.list_of_reactions
  counter = 1
  for (index,reaction_object) in enumerate(list_of_reactions)

    reaction_string = reaction_object.reaction_name
    reaction_type = reaction_object.reaction_type

    # Build comment string -
    comment_string = build_reaction_comment_string(reaction_object)
    buffer *= "$(counter) $(reaction_string)::$(comment_string)\n"

    # update the counter -
    counter = counter + 1;
  end

  # type SpeciesObject
  #
  #   species_type::Symbol
  #   species_bound_type::Symbol
  #   species_symbol::AbstractString
  #   species_index::Int
  #   stoichiometric_coefficient::Float64
  #   species_compartment::Symbol
  #
  #   function SpeciesObject()
  #     this = new()
  #   end
  # end

  buffer *= "\n"

  # write the list of species -
  list_of_species::Array{SpeciesObject} = problem_object.list_of_species
  counter = 1
  for (index,species_object) in enumerate(list_of_species)

    species_symbol = species_object.species_symbol
    buffer *= "$(index) $(species_symbol)\n"
  end

  # build the component -
  program_component::ProgramComponent = ProgramComponent()
  program_component.filename = filename
  program_component.buffer = buffer

  # return -
  return (program_component)
end

function build_kinetics_buffer(problem_object::ProblemObject)

  filename = "Kinetics.jl"

  # build the header -
  header_buffer = build_copyright_header_buffer(problem_object)

  # get the comment buffer -
  comment_header_dictionary = problem_object.configuration_dictionary["function_comment_dictionary"]["kinetics_function"]
  function_comment_buffer = build_function_header_buffer(comment_header_dictionary)

  # initialize the buffer -
  buffer = ""
  buffer *= header_buffer
  buffer *= "#\n"
  buffer *= function_comment_buffer
  buffer *= "function Kinetics(time,state_array,data_dictionary)\n"

  buffer *= "\n"
  buffer *= "\t# Get data/parameters from the data_dictionary - \n"
  buffer *= "\tflux_balance_modes_matrix = data_dictionary[\"flux_balance_modes_matrix\"]\n"
  buffer *= "\t(number_of_reactions,number_of_modes) = size(flux_balance_modes_matrix)\n"

  # initialize -
  buffer *= "\n"
  buffer *= "\t# initialize the kinetic_rate_array - \n"
  buffer *= "\tkinetic_rate_array = zeros(number_of_modes)\n"
  
  buffer *= "\n"
  buffer *= "\treturn kinetic_rate_array\n"
  buffer *= "end\n"

  # build the component -
  program_component::ProgramComponent = ProgramComponent()
  program_component.filename = filename
  program_component.buffer = buffer

  # return -
  return (program_component)
end

function build_data_dictionary_buffer(problem_object::ProblemObject)

  filename = "DataDictionary.jl"

  # build the header -
  header_buffer = build_copyright_header_buffer(problem_object)

  # get the comment buffer -
  comment_header_dictionary = problem_object.configuration_dictionary["function_comment_dictionary"]["data_dictionary_function"]
  function_comment_buffer = build_function_header_buffer(comment_header_dictionary)

  # Get the default -
  # default_parameter_dictionary = problem_object.configuration_dictionary["default_parameter_dictionary"]
  # enzyme_initial_condition = parse(Float64,default_parameter_dictionary["default_protein_initial_condition"])
  # default_rate_constant = parse(Float64,default_parameter_dictionary["default_enzyme_kcat"])
  # default_upper_bound = default_rate_constant*enzyme_initial_condition

  # initialize the buffer -
  buffer = ""
  buffer *= header_buffer
  buffer *= "#\n"
  buffer *= function_comment_buffer
  buffer *= "function DataDictionary(time_start,time_stop,time_step)\n"
  buffer *= "\n"
  buffer *= "\t# Load the networks from disk - \n"
  buffer *= "\tstoichiometric_matrix = readdlm(\"Network.dat\");\n"
  buffer *= "\tflux_balance_modes_matrix = readdlm(\"Modes.dat\");\n"

  buffer *= "\n"
  buffer *= "\t# How many modes do we have? - \n"
  buffer *= "\t(number_of_reactions,number_of_modes) = size(flux_balance_modes_matrix)\n"

  # setup initial condition array -
  buffer *= "\n"
  buffer *= "\t# Setup the initial_condition_array - \n"
  buffer *= "\tenzyme_initial_conditions = ones(number_of_modes)\n"
  buffer *= "\tinitial_condition_array = [\n"
  buffer *= "\t\tenzyme_initial_conditions\t;\t# enzymes\n"

  list_of_species::Array{SpeciesObject} = problem_object.list_of_species
  counter = 1
  for (index,species_object) in enumerate(list_of_species)

    # Get the bound type, and species -
    species_bound_type = species_object.species_bound_type
    species_symbol = species_object.species_symbol

    if (species_bound_type == :unbalanced)
      buffer *= "\t\t0.0\t;\t# $(counter) $(species_symbol)\n"
      counter = counter + 1
    end
  end

  buffer *= "\t]\n"

  # setup the indexes of external species -
  buffer *= "\n"
  buffer *= "\t# Setup the index_vector_external_species - \n"
  buffer *= "\tindex_vector_external_species = [\n"
  counter = 1
  for (index,species_object) in enumerate(list_of_species)

    # Get the bound type, and species -
    species_bound_type = species_object.species_bound_type
    species_symbol = species_object.species_symbol

    if (species_bound_type == :unbalanced)
      buffer *= "\t\t$(index)\t;\t# $(counter) $(species_symbol)\n"
      counter = counter + 1
    end
  end
  buffer *= "\t]\n"

  # setup the degrdation constant array -
  buffer *= "\n"
  buffer *= "\t# Setup the degradation_constant_array - \n"
  buffer *= "\tdefault_protein_half_life = 24.0\t# units:hr\n"
  buffer *= "\tdefault_degrdation_constant = -(1/protein_half_life)*log(0.5)\t# units: hr^-1\n"
  buffer *= "\tdegradation_constant_array = default_degrdation_constant*ones(number_of_modes)\n"

  # setup the external stoichiometric_matrix -
  buffer *= "\n"
  buffer *= "\t# Setup the external stoichiometric matrix - \n"
  buffer *= "\texternal_stoichiometric_matrix = stoichiometric_matrix[index_vector_external_species,:]\n"

  # return block -
  buffer *= "\n"
  buffer *= "\t# =============================== DO NOT EDIT BELOW THIS LINE ============================== #\n"
  buffer *= "\tdata_dictionary = Dict{AbstractString,Any}()\n"
  buffer *= "\tdata_dictionary[\"stoichiometric_matrix\"] = stoichiometric_matrix\n"
  buffer *= "\tdata_dictionary[\"external_stoichiometric_matrix\"] = external_stoichiometric_matrix\n"
  buffer *= "\tdata_dictionary[\"initial_condition_array\"] = initial_condition_array\n"
  buffer *= "\tdata_dictionary[\"degradation_constant_array\"] = degradation_constant_array\n"
  buffer *= "\tdata_dictionary[\"flux_balance_modes_matrix\"] = flux_balance_modes_matrix\n"
  buffer *= "\tdata_dictionary[\"index_vector_external_species\"] = index_vector_external_species\n"
  buffer *= "\t# =============================== DO NOT EDIT ABOVE THIS LINE ============================== #\n"
  buffer *= "\treturn data_dictionary\n"
  buffer *= "end\n"

  # build the component -
  program_component::ProgramComponent = ProgramComponent()
  program_component.filename = filename
  program_component.buffer = buffer

  # return -
  return (program_component)
end
