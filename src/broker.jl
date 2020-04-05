using ZMQ
using ZeroMQ_jll
using ArgParse


function start_broker(worker_url::String, client_url::String)

    # ctx = Context()
    router = Socket(ROUTER)
    dealer = Socket(DEALER)

    ZMQ.bind(router, client_url)
    ZMQ.bind(dealer, worker_url)

    # missing: ZMQ.proxy(router, dealer)
    try
        try
            @info "worker=$worker_url client=$client_url"
            rc = ccall((:zmq_proxy, libzmq), Cint,  (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}), router, dealer, C_NULL)
            @info "done zmq_proxy $(rc)"
        catch
            @info "exception!"
        end
    finally
        # control never comes here... clean up anyway.
        ZMQ.close(router)
        ZMQ.close(dealer)
        # ZMQ.close(ctx)
    end
end

function broker_args()
    broker_args = ArgParseSettings(prog = "broker", autofix_names = true)  # turn "-" into "_" for arg names.

    @add_arg_table! broker_args begin
        
        "--worker"
        arg_type = String
        metavar = "URL"
        default = "tcp://127.0.0.1:9467"
        help = "ZMQ DEALER address to connect to"

        "--client"
        arg_type = String
        metavar = "URL"
        default = "ipc:///tmp/chloe-client"
        help = "ZMQ ROUTER address to connect to"
    end

broker_args = parse_args(ARGS, broker_args; as_symbols = true)
    broker_args
end
function broker_main()
    args = broker_args()
    Sys.set_process_title("chloe-broker")
    start_broker(args[:worker], args[:client])
end

if abspath(PROGRAM_FILE) == @__FILE__
    broker_main()
end