
import ArgParse: ArgParseSettings, @add_arg_table!, parse_args
using JuliaWebAPI

const ADDRESS = "ipc:///tmp/chloe-client"
 
function chloe_client(;fasta = String[], address = ADDRESS, output::String)
    invoker = APIInvoker(address)
    res = apicall(invoker, "chloe", fasta[1], output)
    if res["code"] !== 200
        println("failed")
    else
        println("result: $(res["data"])")
    end

end


# const julia_v07 = VERSION > v"0.7-"

client_args = ArgParseSettings(prog = "Chloë", autofix_names = true)  # turn "-" into "_" for arg names.

@add_arg_table! client_args begin
    "fasta"
        arg_type = String
        nargs = 1
        required = true
        help = "fasta file to process"
    "--output", "-o"
        arg_type = String
        required = true
        help = "output .sff filename"
    "--address", "-a"
        arg_type = String
        default = ADDRESS
        help = "ZMQ address to connect to"

end

client_args.epilog = """
Annotate a fasta file unsing Chloe server
"""

function client_main() 
    parsed_args = parse_args(ARGS, client_args; as_symbols = true)
    # filter!(kv->kv.second ∉ (nothing, false), parsed_args)
    chloe_client(;parsed_args...)
end


if abspath(PROGRAM_FILE) == @__FILE__
    client_main()
end