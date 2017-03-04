function parse_vff_measured_species_statements(path_to_model_file::AbstractString)

  # initialize -
  measured_species_vector = SpeciesObject[]
  tmp_array::Array{AbstractString} = AbstractString[]
  desired_handler_symbol::Symbol = :measured_species_handler
  current_handler_symbol::Symbol = :metabolic_reaction_handler

  try

    # Open the model file, and read each line into a vector -
    open(path_to_model_file,"r") do model_file
      for line in eachline(model_file)

          if (contains(line,"//") == false && search(line,"\n")[1] != 1)
            push!(tmp_array,chomp(line))
          end
      end
    end

    counter = 1
    for sentence in tmp_array

      # Check the current pragma -
      if (contains(sentence,"#pragma") == true)

        # Extract the current handler -
        current_handler_symbol = extract_vff_handler_symbol(sentence)
      end

      # Depending upon the handler_symbol, we will process diff types of records -
      if (current_handler_symbol == desired_handler_symbol && contains(sentence,"#pragma") == false)

        measured_species_object = vff_measured_species_object_factory(sentence,current_handler_symbol)
        measured_species_object.species_index = counter
        push!(measured_species_vector,measured_species_object)

        # Update the counter -
        counter = counter + 1
      end
    end

  catch err
    showerror(STDOUT, err, backtrace());println()
  end

  return measured_species_vector
end

function parse_vff_flux_mode_statements(path_to_flux_mode_file::AbstractString)

  # Create an array of mode sentences -
  sentence_vector = VFFFluxModeSentence[]
  tmp_array::Array{AbstractString} = AbstractString[]

  # open the file -
  try

    # Open the model file, and read each line into a vector -
    open(path_to_flux_mode_file,"r") do model_file
      for line in eachline(model_file)

          if (contains(line,"//") == false && search(line,"\n")[1] != 1)
            push!(tmp_array,chomp(line))
          end
      end
    end

    for sentence in tmp_array

      vff_flux_model_sentence = vff_flux_mode_sentence_factory(sentence)
      push!(sentence_vector,vff_flux_model_sentence)
    end

  catch err
    showerror(STDOUT, err, backtrace());println()
  end

  return sentence_vector
end

function parse_vff_metabolic_statements(path_to_model_file::AbstractString)

  # We are going to load the sentences in the file into a vector
  # if not a valid model file, then throw an error -
  sentence_vector = VFFSentence[]
  expanded_sentence_vector::Array{VFFSentence} = VFFSentence[]
  tmp_array::Array{AbstractString} = AbstractString[]
  desired_handler_symbol::Symbol = :metabolic_reaction_handler
  current_handler_symbol::Symbol = :metabolic_reaction_handler

  try

    # Open the model file, and read each line into a vector -
    open(path_to_model_file,"r") do model_file
      for line in eachline(model_file)

          if (contains(line,"//") == false && search(line,"\n")[1] != 1)
            push!(tmp_array,chomp(line))
          end
      end
    end

    for sentence in tmp_array

      # Check the current pragma -
      if (contains(sentence,"#pragma") == true)

        # Extract the current handler -
        current_handler_symbol = extract_vff_handler_symbol(sentence)
      end

      # Depending upon the handler_symbol, we will process diff types of records -
      if (current_handler_symbol == desired_handler_symbol && contains(sentence,"#pragma") == false)

        local_vff_sentence_array = vff_metabolic_sentence_factory(sentence,current_handler_symbol)
        for local_vff_sentence in collect(local_vff_sentence_array)
          push!(expanded_sentence_vector,local_vff_sentence)
        end
      end
    end

  catch err
    showerror(STDOUT, err, backtrace());println()
  end

  return expanded_sentence_vector
end

# ========================================================================================= #
# Helper functions -
# ========================================================================================= #
function extract_vff_handler_symbol(sentence::String)

  # Default - reaction handler -
  handler_symbol = :metabolic_reaction_handler

  # split the sentence -
  split_array = split(sentence,"::")

  # grab handler -
  handler_string = split_array[2]
  if (handler_string == "metabolic_reaction_handler")
    handler_symbol = :metabolic_reaction_handler
  end

  return handler_symbol

end

function vff_flux_mode_sentence_factory(sentence::String)

  vff_flux_model_sentence = VFFFluxModeSentence()
  vff_flux_model_sentence.original_sentence = sentence

  @show sentence

  # split the sentence -
  split_array = split(sentence,",")
  metabolite_source_array = String[]
  for (index,token) in enumerate(split_array)

    if (index == 1)
      vff_flux_model_sentence.sentence_name = token
    elseif (index == 2)
      vff_flux_model_sentence.pivot_index = parse(Int64,token)
    else
      push!(metabolite_source_array,token)
    end
  end

  # grab the array -
  vff_flux_model_sentence.species_symbol_array = metabolite_source_array

  # return -
  return vff_flux_model_sentence
end

function vff_metabolic_sentence_factory(sentence::String,handler_symbol::Symbol)

  sentence_vector = VFFSentence[]

  # Ok, so now we have the array for sentences -
  vff_sentence = VFFSentence()
  vff_sentence.original_sentence = sentence

  @show sentence

  # split the sentence -
  split_array = split(sentence,",")

  # sentence_name::AbstractString
  # sentence_reactant_clause::AbstractString
  # sentence_product_clause::AbstractString
  # sentence_reverse_bound::Float64
  # sentence_forward_bound::Float64
  # sentence_delimiter::Char
  vff_sentence.sentence_name = split_array[1]
  vff_sentence.sentence_reactant_clause = split_array[2]
  vff_sentence.sentence_product_clause = split_array[3]
  vff_sentence.sentence_reverse_bound = parse(Float64,split_array[4])
  vff_sentence.sentence_forward_bound = parse(Float64,split_array[5])
  vff_sentence.sentence_handler = handler_symbol
  vff_sentence.sentence_delimiter = ','

  # add sentence to sentence_vector -
  push!(sentence_vector,vff_sentence)

  # Check - is this reversible?
  if (vff_sentence.sentence_reverse_bound == -Inf)

    # ok, so we have a reversible reaction -
    # first change lower bound to 0 -
    vff_sentence.sentence_reverse_bound = 0.0

    # create a new copy of sentence -
    vff_sentence_copy = copy(vff_sentence)
    vff_sentence_copy.sentence_name = (vff_sentence_copy.sentence_name)*"_reverse"
    vff_sentence_copy.sentence_reactant_clause = vff_sentence.sentence_product_clause
    vff_sentence_copy.sentence_product_clause = vff_sentence.sentence_reactant_clause
    vff_sentence_copy.sentence_handler = vff_sentence.sentence_handler
    vff_sentence_copy.sentence_type_flag = vff_sentence.sentence_type_flag

    # add sentence and sentence copy to the expanded_sentence_vector -
    push!(sentence_vector,vff_sentence_copy)
  end

  return sentence_vector
end
