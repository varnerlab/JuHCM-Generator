
type VFFFluxModeSentence

  original_sentence::AbstractString
  sentence_name::AbstractString
  species_symbol_array::Array{String}
  pivot_index::Int

  function VFFFluxModeSentence()
    this = new()
  end

end

type VFFSentence

  original_sentence::AbstractString
  sentence_name::AbstractString
  sentence_type_flag::Int
  sentence_reactant_clause::AbstractString
  sentence_product_clause::AbstractString
  sentence_reverse_bound::Float64
  sentence_forward_bound::Float64
  sentence_delimiter::Char
  sentence_handler::Symbol

  function VFFSentence()
    this = new()
  end
end


type ProgramComponent

  filename::AbstractString
  buffer::AbstractString

  function ProgramComponent()
    this = new()
  end

end

type SpeciesObject

  species_type::Symbol
  species_bound_type::Symbol
  species_symbol::AbstractString
  species_index::Int
  stoichiometric_coefficient::Float64
  species_compartment::Symbol

  function SpeciesObject()
    this = new()
  end
end

type ReactionObject

  enyzme_generation_flag::Int
  reaction_type::Symbol
  reaction_index::Int
  is_reaction_reversible::Bool
  reaction_name::AbstractString
  list_of_reactants::Array{SpeciesObject}
  list_of_products::Array{SpeciesObject}

  function ReactionObject()
    this = new()
  end
end


type ProblemObject

  configuration_dictionary::Dict{AbstractString,Any}
  list_of_species::Array{SpeciesObject}
  list_of_reactions::Array{ReactionObject}
  list_of_flux_modes::Array{VFFFluxModeSentence}

  function ProblemObject()
    this = new()
  end
end
