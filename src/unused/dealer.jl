import ZMQ
import ZeroMQ_jll
using JuliaWebAPI

# doesn't seem to work!
# trying to use inproc://worker ZMQ to create a DEALER/ROUTER
# *within* the chloe server process.

# But julia ZMQ is missing the proxy function and my attempts
# to monkeypatch it doesn't seem to work... so
# I need a second process to run the DEALER/ROUTER ... see bin/broker.py

function ping()
    return "OK $(current_task())"
end
function worker(async)
    process(
        JuliaWebAPI.create_responder([
                (ping, false)

            ], "inproc://workers", false,  "chloe"); async=async
        )
    @info "end worker"
end

function worker2(name)
    api = APIResponder(InProcTransport(Symbol(name)), DictMsgFormat(), "chloe", false)
    register(api, ping)
    process(
        api;async=true
    )
end
function invoker(name)
    APIInvoker(InProcTransport(Symbol(name)), DictMsgFormat())
end

function multi()
    channels = Array{Tuple{APIInvoker{InProcTransport,DictMsgFormat},APIResponder{InProcTransport,DictMsgFormat}}}(undef, 0)
    for i in 1:5
        name = "worker$i"
        push!(channels, (invoker(name), worker2(name)))
    end

end

function start_broker(router_url::String, dealer_url::String)

    # ctx = Context()
    router = ZMQ.Socket(ZMQ.ROUTER)
    dealer = ZMQ.Socket(ZMQ.DEALER)

    ZMQ.bind(router, router_url)
    ZMQ.bind(dealer, dealer_url)

    # missing: ZMQ.proxy(router, dealer)
    rc = ccall((:zmq_proxy, ZeroMQ_jll.libzmq), Cint,  (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}), router, dealer, C_NULL)
    @info "done proxy $(rc)"

    # control never comes here... clean up anyway.
    ZMQ.close(router)
    ZMQ.close(dealer)
    # ZMQ.close(ctx)
end

funcs = [worker, worker]
function workers()
    @sync for f in funcs
        @async f()
    end
end



start_broker("ipc:///tmp/chloe-client", "tcp://127.0.0.1:9467")