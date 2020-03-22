using ZMQ
using ZeroMQ_jll
using JuliaWebAPI

# doesn't seem to work!

# function default_endpoint(f::Function)
#     endpt = string(f)
#     # separate the module (more natural URL, assumes 'using Module')
#     if '.' in endpt
#         endpt = rsplit(endpt, '.', limit = 2)[2]
#     end
#     endpt
# end
# function _add_spec(spec::Tuple, api::APIResponder)
#     fn = spec[1]
#     resp_json = (length(spec) > 1) ? spec[2] : false
#     resp_headers = (length(spec) > 2) ? spec[3] : Dict{String,String}()
#     api_name = (length(spec) > 3) ? spec[4] : default_endpoint(fn)
#     register(api, fn, resp_json = resp_json, resp_headers = resp_headers, endpt = api_name)
# end

# function create_responder(apispecs::Array, addr, bind, nid, ctx::Context, open = false)
#     api = APIResponder(ZMQTransport(addr, REP, bind, ctx), JSONMsgFormat(), nid, open)
#     for spec in apispecs
#         _add_spec(spec, api)
#     end
#     api
# end
function ping()
    return "OK $(Threads.threadid())"
end
function worker()
    bind = false
    process(
        JuliaWebAPI.create_responder([
                (ping, false)

            ], "inproc://workers", bind,  "chloe"); async = false
        )
    println("end worker")
end


function start_broker(url::String)

    ctx = Context()
    router = Socket(ctx, ROUTER)
    dealer = Socket(ctx, DEALER)

    ZMQ.bind(router, url)
    ZMQ.bind(dealer, "inproc://workers")
    # ZMQ.bind(dealer, "tcp://127.0.0.1:9998")

    @info "listening on $(url)"
    rc = ccall((:zmq_proxy, libzmq), Cint,  (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}), router, dealer, C_NULL)
    # missing: ZMQ.proxy(router, dealer)
    println("done $(rc)")

    # control never comes here
    ZMQ.close(router)
    ZMQ.close(dealer)
    ZMQ.close(ctx)
end

funcs = [worker, worker, ()->start_broker("tcp://127.0.0.1:9999") ]

Threads.@threads for f in funcs
    f()
end