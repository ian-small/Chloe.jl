
module WebAPI
export TerminatingJSONMsgFormat

import JuliaWebAPI
import JuliaWebAPI: JSONMsgFormat, AbstractMsgFormat, juliaformat, wireformat

# overide JuliaWebAPI interface

struct TerminatingJSONMsgFormat <: AbstractMsgFormat
end

const IJSONMsgFormat = JSONMsgFormat()

JuliaWebAPI.wireformat(fmt::TerminatingJSONMsgFormat, 
    cmd::String, args...; data...) = wireformat(IJSONMsgFormat, cmd, args... ; data...)
JuliaWebAPI.wireformat(fmt::TerminatingJSONMsgFormat, code::Int, 
    headers::Dict{String,String}, resp, id=nothing) = wireformat(IJSONMsgFormat, code, headers, resp, id)


function JuliaWebAPI.juliaformat(fmt::TerminatingJSONMsgFormat, msgstr)
    ret = juliaformat(IJSONMsgFormat, msgstr)
    if get(ret, "cmd", "") == ":exit"
        pid = get(ret, "args", [])
        if length(pid) != 1
            error("wrong arguments")
        end
        if getpid() != pid[1]
            error("exit: incorrect process")
        end
        ret["cmd"] = ":terminate"
    end
    return ret
end

end
