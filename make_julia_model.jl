include("./src/Include.jl")

# Parse the command line options -
function parse_commandline()
    settings_object = ArgParseSettings()
    @add_arg_table settings_object begin
      "-o"
        help = "Directory where the Matlab model files will be written."
        arg_type = AbstractString
        default = "."

      "-n"
        help = "Path to the biochemical reaction file written in the vff-reaction format."
        arg_type = AbstractString
        required = true

      "-m"
        help = "Path to the flux model file written in the vff-mode format."
        arg_type = AbstractString
        required = true
    end

    # return a dictionary w/args -
    return parse_args(settings_object)
end

function main()

  # Build the arguement dictionary -
  parsed_args = parse_commandline()

  # Load the JSON configuration file -
  config_dict = JSON.parsefile("./config/Configuration.json")

  # hcm specific logic goes here -
  component_set = Set{ProgramComponent}()

  # Load the statement_vector -
  path_to_model_file = parsed_args["n"]
  metabolic_statement_vector::Array{VFFSentence} = parse_vff_metabolic_statements(path_to_model_file)

  # load the flux model statement vector -
  path_to_flux_mode_file = parsed_args["m"]
  flux_mode_statement_vector::Array{VFFFluxModeSentence} = parse_vff_flux_mode_statements(path_to_flux_mode_file)

  # Generate the problem object -
  problem_object = generate_problem_object(metabolic_statement_vector,flux_mode_statement_vector,config_dict)
  problem_object.configuration_dictionary = config_dict

  # Write the stoichiometric_matrix --
  program_component_stoichiometric_matrix = generate_stoichiomteric_matrix_buffer(problem_object)
  push!(component_set,program_component_stoichiometric_matrix)

  # Write the kinetics buffer --
  program_component_kinetics_buffer = build_kinetics_buffer(problem_object)
  push!(component_set,program_component_kinetics_buffer)

  # Write the data dictionary -
  program_component_data_dictionary = build_data_dictionary_buffer(problem_object)
  push!(component_set,program_component_data_dictionary)

  # write debug buffer -
  program_component_debug = build_debug_buffer(problem_object)
  push!(component_set,program_component_debug)

  # Dump the component_set to disk -
  path_to_output_file = parsed_args["o"]
  write_program_components_to_disk(path_to_output_file,component_set)

  # Transfer distrubtion jl files to the output -
  transfer_distribution_files("./distribution",path_to_output_file,".jl")
end

# call main ... and go ...
main()
