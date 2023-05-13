
export format_camelcase, format_snakecase

"""
    format_camelcase(str::AbstractString)

Format `str` in CamelCase, such that the first letter of each word 
in the sentence is capitalized and spaces are removed.
"""
function format_camelcase(str::S) where {S <: AbstractString}
    format_camelcase(S, str)
end

function format_camelcase(::Type{T}, s::S) where {T,S<:AbstractString}
    words = split(s, r"[\_, \s, \-]")
    if length(words) > 1
        T(join(uppercasefirst.(lowercase.(words))))
    elseif length(words) == 1
        T(uppercasefirst(s))
    else
        throw(error("$s cannot be formatted in CamelCase."))
    end
end


"""
    format_snakecase(str::AbstractString)

Format `str` in SnakeCase, such that all the letters are in lower case and 
spaces are replaced with underscores.
"""
function format_snakecase(str::S) where {S <: AbstractString}
    format_snakecase(S, str)
end

function format_snakecase(::Type{T}, s::S) where {T,S<:AbstractString}
    return T(join(lowercase.(split(replace(s, r"[\-\.\s]" => "_"), r"(?=[A-Z])")), "_"))
end
