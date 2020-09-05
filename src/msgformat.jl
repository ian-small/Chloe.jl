
# using JSON

using JuliaWebAPI
import JuliaWebAPI: juliaformat, AbstractMsgFormat, juliaformat, wireformat, JSONMsgFormat

# overide JuliaWebAPI interface

struct JSONMsgFormatEx <: AbstractMsgFormat
end

const IJSONMsgFormat = JSONMsgFormat()

JuliaWebAPI.wireformat(fmt::JSONMsgFormatEx, 
    cmd::String, args...; data...) = wireformat(IJSONMsgFormat, cmd, args... ; data...)
JuliaWebAPI.wireformat(fmt::JSONMsgFormatEx, code::Int, 
    headers::Dict{String,String}, resp, id=nothing) = wireformat(IJSONMsgFormat, code, headers, resp, id)


function JuliaWebAPI.juliaformat(fmt::JSONMsgFormatEx, msgstr)
    @info "JSONMsgFormatEx: calling $(msgstr)"
    juliaformat(IJSONMsgFormat, msgstr)
end
