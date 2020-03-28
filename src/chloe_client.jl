
using ArgParse
using JuliaWebAPI

const ADDRESS = "tcp://127.0.0.1:9999"
 
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

args = ArgParseSettings(prog = "Chloë", autofix_names = true)  # turn "-" into "_" for arg names.

@add_arg_table! args begin
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
        help = "ZMQ address to listen on"

end

args.epilog = """
Annotate a fasta file unsing Chloe server
"""

function real_main() 
    parsed_args = parse_args(ARGS, args; as_symbols = true)
    # filter!(kv->kv.second ∉ (nothing, false), parsed_args)
    chloe_client(;parsed_args...)
end


if abspath(PROGRAM_FILE) == @__FILE__
    real_main()
end