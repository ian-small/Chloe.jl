using ZMQ
using ZeroMQ_jll
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
function worker()
    process(
        JuliaWebAPI.create_responder([
                (ping, false)

            ], "inproc://workers", false,  "chloe"); async = false
        )
    @info "end worker"
end


function start_broker(url::String)

    # ctx = Context()
    router = Socket(ROUTER)
    dealer = Socket(DEALER)

    ZMQ.bind(router, url)
    ZMQ.bind(dealer, "inproc://workers")

    @info "listening on $(url)"
    # missing: ZMQ.proxy(router, dealer)
    rc = ccall((:zmq_proxy, libzmq), Cint,  (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}), router, dealer, C_NULL)
    println("done $(rc)")

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
Threads.@threads for w in [workers, ()->start_broker("tcp://127.0.0.1:9998")]
    w()     
end