function write_program_components_to_disk(file_path::AbstractString,set_of_program_components::Set{ProgramComponent})

  # go through each component, and dump the buffer to disk -
  for program_component in set_of_program_components

    # Get the data -
    filename = program_component.filename
    program_buffer = program_component.buffer

    # build the path -
    path_to_program_file = file_path*"/"*filename

    # Write the file -
    outfile = open(path_to_program_file, "w")
    write(outfile,program_buffer);
    close(outfile);
  end
end

function build_reaction_comment_string(reaction_object::ReactionObject)

  # Ok, let's build a comment from the reactioj object {Reactants} -> {Products}
  list_of_reactants = reaction_object.list_of_reactants
  list_of_products = reaction_object.list_of_products

  # Build reactant string -
  reactant_buffer = ""
  for species_object::SpeciesObject in list_of_reactants

    # symbol etc -
    species_symbol = species_object.species_symbol
    stoichiometric_coefficient = species_object.stoichiometric_coefficient

    if (stoichiometric_coefficient == 1.0)
      # fill reactant buffer -
      reactant_buffer *= "$(species_symbol)+"
    else
      # fill reactant buffer -
      reactant_buffer *= "$(stoichiometric_coefficient)\*$(species_symbol)+"
    end
  end

  # Cutoff the trailing +
  reactant_buffer = reactant_buffer[1:end-1]

  # Check for empty ...
  if (length(list_of_reactants) == 0)
    reactant_buffer = "[]"
  end

  # Build product string -
  product_buffer = ""
  for species_object::SpeciesObject in list_of_products

    # symbol etc -
    species_symbol = species_object.species_symbol
    stoichiometric_coefficient = species_object.stoichiometric_coefficient

    if (stoichiometric_coefficient == 1.0)
      # fill product buffer -
      product_buffer *= "$(species_symbol)+"
    else
      # fill product buffer -
      product_buffer *= "$(stoichiometric_coefficient)\*$(species_symbol)+"
    end
  end

  # Cutoff the trailing +
  product_buffer = product_buffer[1:end-1]

  # Check for empty ...
  if (length(list_of_products) == 0)
    product_buffer = "[]"
  end

  # comment string -
  comment_buffer = ""
  comment_buffer = reactant_buffer*" --> "*product_buffer
  return comment_buffer
end

function transfer_distribution_file(path_to_distribution_files::AbstractString,
                                      input_file_name_with_extension::AbstractString,
                                      path_to_output_files::AbstractString,
                                      output_file_name_with_extension::AbstractString)

  # Load the specific file -
  # create src_buffer -
  src_buffer::Array{AbstractString} = AbstractString[]

  # path to distrubtion -
  path_to_src_file = path_to_distribution_files*"/"*input_file_name_with_extension
  open(path_to_src_file,"r") do src_file
    for line in eachline(src_file)
      push!(src_buffer,line)
    end
  end

  # Write the file to the output -
  path_to_program_file = path_to_output_files*"/"*output_file_name_with_extension
  outfile = open(path_to_program_file, "w")
  write(outfile,src_buffer);
  close(outfile);

end

function transfer_distribution_files(path_to_distribution_files::AbstractString,
                                      path_to_output_files::AbstractString,
                                      file_extension::AbstractString)


  # Search the directory for src files -
  # load the files -
  searchdir(path,key) = filter(x->contains(x,key),readdir(path))

  # build src file list -
  list_of_src_files = searchdir(path_to_distribution_files,file_extension)

  # go thru the src file list, and copy the files to the output path -
  for src_file in list_of_src_files

    # create src_buffer -
    src_buffer::Array{AbstractString} = AbstractString[]

    # path to distrubtion -
    path_to_src_file = path_to_distribution_files*"/"*src_file
    open(path_to_src_file,"r") do src_file
      for line in eachline(src_file)
        push!(src_buffer,line)
      end
    end

    # Write the file to the output -
    path_to_program_file = path_to_output_files*"/"*src_file
    outfile = open(path_to_program_file, "w")
    write(outfile,src_buffer);
    close(outfile);
  end
end

function update_reaction_list!(list_of_reactions::Array{ReactionObject},
  list_of_target_species::Array{SpeciesObject},
  target_bound_symbol::Symbol)

  # Create a set of free species symbols -
  set_target_species_symbol = Set{String}()
  for (index,target_species_object) in enumerate(list_of_target_species)

    # Grab the symbol -
    species_symbol = target_species_object.species_symbol
    push!(set_target_species_symbol,species_symbol)
  end

  # Iterate -
  for (index,reaction_object) in enumerate(list_of_reactions)

    # grab the reactants -
    list_of_reactants::Array{SpeciesObject} = reaction_object.list_of_reactants
    for (reactant_index,reactant_object) in enumerate(list_of_reactants)

      # Get the symbol for the reactants -
      reactant_symbol = reactant_object.species_symbol
      if (in(reactant_symbol,set_target_species_symbol) == true)
        reactant_object.species_bound_type = target_bound_symbol
      end
    end

    # grab the products -
    list_of_products::Array{SpeciesObject} = reaction_object.list_of_products
    for (product_index,product_object) in enumerate(list_of_products)

      # Get the symbol for the reactants -
      product_symbol = product_object.species_symbol
      if (in(product_symbol,set_target_species_symbol) == true)
        product_object.species_bound_type = target_bound_symbol
      end
    end
  end
end

function classify_species_bounds!(list_of_species::Array{SpeciesObject},
  list_of_target_species::Array{SpeciesObject},
  target_bound_symbol::Symbol)

  # Create a set of free species symbols -
  set_target_species_symbol = Set{String}()
  for (index,target_species_object) in enumerate(list_of_target_species)

    # Grab the symbol -
    species_symbol = target_species_object.species_symbol
    push!(set_target_species_symbol,species_symbol)
  end

  # iterate -
  for (index,local_species_object) in enumerate(list_of_species)

    # what is species symbol -
    local_species_symbol = local_species_object.species_symbol
    if (in(local_species_symbol,set_target_species_symbol) == true)
      local_species_object.species_bound_type = target_bound_symbol
    end
  end
end

function partition!(list_of_reactions::Array{ReactionObject})

  # ok, we need to split the reaction list into solved and kinetic -
  list_of_solved_indexes::Array{Int} = Int[]
  list_of_kinetic_indexes::Array{Int} = Int[]

  for (index,reaction_object) in enumerate(list_of_reactions)

    is_kinetic_reaction::Bool = false

    # get the list of reactants -
    list_of_reactants::Array{SpeciesObject} = reaction_object.list_of_reactants
    for (local_index,species_object::SpeciesObject) in enumerate(list_of_reactants)

      if (species_object.species_bound_type == :free)
          is_kinetic_reaction = true
      end
    end

    if (is_kinetic_reaction == true)
      push!(list_of_kinetic_indexes,index)
    else
      push!(list_of_solved_indexes,index)
    end
  end # end-for

  # Ok, so we have two index lists, cat them -
  permutation_index_array = vcat(list_of_solved_indexes,list_of_kinetic_indexes)

  # permute the array -
  permute!(list_of_reactions,permutation_index_array)

  # how many solved reactions do we have?
  number_of_solved_reactions = length(list_of_solved_indexes)
  for index in 1:number_of_solved_reactions
    reaction_object = list_of_reactions[index]
    reaction_object.reaction_type = :solved
  end

  # how many kinetic rates do we have?
  number_of_kinetic_rates = length(list_of_kinetic_indexes)
  for index in 1:number_of_kinetic_rates
    translated_index = index+number_of_solved_reactions
    reaction_object = list_of_reactions[translated_index]
    reaction_object.reaction_type = :kinetic
  end
end

function partition!(list_of_species::Array{SpeciesObject})

  # ok, frist, we need to split into balanced and unbalanced lists -
  list_of_balanced_indexes::Array{Int} = Int[]
  list_of_free_indexes::Array{Int} = Int[]

  for (index,species_object) in enumerate(list_of_species)

    # what is bounds symbol -
    species_bound_type::Symbol = species_object.species_bound_type
    if (species_bound_type == :balanced)
      push!(list_of_balanced_indexes,index)
    elseif (species_bound_type == :unbalanced)
      push!(list_of_free_indexes,index)
    end
  end

  # combine -
  permutation_index_array = vcat(list_of_balanced_indexes,list_of_free_indexes)

  # permute the array -
  permute!(list_of_species,permutation_index_array)
end


function extract_reaction_name_set(list_of_reactions::Array{ReactionObject})

  reaction_name_set = Set{String}()
  for (index,reaction_object) in enumerate(list_of_reactions)

    local_reaction_name = reaction_object.reaction_name
    push!(reaction_name_set,local_reaction_name)
  end

  return reaction_name_set
end


function find_index_of_reaction_with_name(reaction_array::Array{ReactionObject},reaction_name::String)

  desired_reaction_index = -1
  for (reaction_index,reaction_object) in enumerate(reaction_array)

    # get reaction name -
    test_reaction_name = reaction_object.reaction_name
    if (test_reaction_name == reaction_name)
      desired_reaction_index = reaction_index
      break
    end
  end

  return desired_reaction_index
end

function copy(sentence::VFFSentence)

  # create a new sentence -
  sentence_copy::VFFSentence = VFFSentence()

  # sentence_name::AbstractString
  # sentence_reactant_clause::AbstractString
  # sentence_product_clause::AbstractString
  # sentence_reverse_bound::Float64
  # sentence_forward_bound::Float64
  # sentence_delimiter::Char
  sentence_copy.original_sentence = sentence.original_sentence
  sentence_copy.sentence_name = sentence.sentence_name
  sentence_copy.sentence_reactant_clause = sentence.sentence_reactant_clause
  sentence_copy.sentence_product_clause = sentence.sentence_product_clause
  sentence_copy.sentence_reverse_bound = sentence.sentence_reverse_bound
  sentence_copy.sentence_forward_bound = sentence.sentence_forward_bound
  sentence_copy.sentence_delimiter = sentence.sentence_delimiter

  return sentence_copy
end

function generate_stoichiomteric_matrix_buffer(problem_object::ProblemObject)

  # From the problem object - get the data for the stoichiometric_matrix -
  # ...
  # ...

  # Build the buffer =
  buffer = ""

  # build the buffer -
  list_of_species::Array{SpeciesObject} = problem_object.list_of_species
  list_of_reactions::Array{ReactionObject} = problem_object.list_of_reactions
  for species_object in list_of_species

    #@show species_object

      for reaction_object in list_of_reactions

        #@show reaction_object

        # is this species involved in this reaction?
        if (is_species_a_reactant_in_reaction(species_object,reaction_object))

          # ok, we have a *reactant* - get the st coefficient, multiply by -1 and add it to the buffer
          stoichiometric_coefficient = get_stoichiometric_coefficient(species_object,reaction_object,:reactant)
          buffer *= " -$(stoichiometric_coefficient) "

        elseif (is_species_a_product_in_reaction(species_object,reaction_object))

          # ok, we have a *product* - get the st coefficient, and add it to the buffer
          stoichiometric_coefficient = get_stoichiometric_coefficient(species_object,reaction_object,:product)
          buffer *= " $(stoichiometric_coefficient) "
        else

          # this species is *not* involved in this reaction - coefficient is zero
          buffer *= " 0.0 "
        end
      end

      # add the new line -
      buffer *= "\n"
  end


  # build the component -
  filename = "Network.dat"
  program_component::ProgramComponent = ProgramComponent()
  program_component.filename = filename
  program_component.buffer = buffer

  # return -
  return (program_component)
end

function get_stoichiometric_coefficient(species_object::SpeciesObject,reaction_object::ReactionObject,direction::Symbol)

  if (direction == :reactant)

    # what is the symbol?
    test_symbol = species_object.species_symbol

    # Find symbol in list, return coeff
    list_of_reactants::Array{SpeciesObject} = reaction_object.list_of_reactants
    for reactant in list_of_reactants

      if (test_symbol == reactant.species_symbol)
          return reactant.stoichiometric_coefficient
      end
    end

  else

    # what is the symbol?
    test_symbol = species_object.species_symbol

    # Get list of species -
    list_of_products::Array{SpeciesObject} = reaction_object.list_of_products
    for product in list_of_products

      if (test_symbol == product.species_symbol)
        return product.stoichiometric_coefficient
      end
    end
  end

  return 0.0
end

function is_species_a_reactant_in_reaction(species_object::SpeciesObject,reaction_object::ReactionObject)

  # Default: false -
  is_species_in_reaction::Bool = false

  # what is the symbol?
  test_symbol = species_object.species_symbol

  # Get list of species -
  list_of_reactants::Array{SpeciesObject} = reaction_object.list_of_reactants
  for reactant in list_of_reactants

    if (test_symbol == reactant.species_symbol)
      return true
    end
  end

  return is_species_in_reaction
end

function is_species_a_product_in_reaction(species_object::SpeciesObject,reaction_object::ReactionObject)

  # Default: false -
  is_species_in_reaction::Bool = false

  # what is the symbol?
  test_symbol = species_object.species_symbol

  # Get list of species -
  list_of_products::Array{SpeciesObject} = reaction_object.list_of_products
  for product in list_of_products

    if (test_symbol == product.species_symbol)
      return true
    end
  end

  return is_species_in_reaction
end
