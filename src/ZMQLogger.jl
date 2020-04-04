using Logging
using ZMQ

struct ZMQLogger <: AbstractLogger
    min_level::Logging.LogLevel
    lock::ReentrantLock
    topic::String
    message_limits::Dict{Any,Int}
    socket::Socket

    function ZMQLogger(endpoint::String, min_level::Logging.LogLevel = Logging.Warn;topic::String = "", message_limits = nothing)
        if message_limits == nothing
            message_limits = Dict{Any,Int}()
        end
        socket = ZMQ.Socket(ZMQ.PUB)
        ZMQ.connect(socket, endpoint)
        new(min_level, ReentrantLock(), topic, message_limits, socket)

    end
end

function Logging.handle_message(logger::ZMQLogger, level, message, _module, group, id,
    filepath, line; maxlog = nothing, kwargs...)
    
    if maxlog !== nothing && maxlog isa Integer
        remaining = get!(logger.message_limits, id, maxlog)
        logger.message_limits[id] = remaining - 1
        remaining > 0 || return
    end
    msglines = split(chomp(string(message)), '\n')
    buf = IOBuffer()
    io = IOContext(buf, stderr)
    println(io, msglines[1])
    for i in 2:length(msglines)
        println(io, "\t> ", msglines[i])
    end
    for (key, val) in kwargs
        println(io, "\t> ", key, " = ", val)
    end
    msg = String(take!(buf))
    prefix = (level == Logging.Warn ? "WARNING" : uppercase(string(level))) 
    topic = length(logger.topic) > 0 ? "$(logger.topic).$prefix" : prefix
    lock(logger.lock) do
        ZMQ.send(logger.socket, Message(topic); more = true)
        ZMQ.send(logger.socket, Message(msg); more = false)
    end
    nothing
end
function Logging.shouldlog(logger::ZMQLogger, level, _module, group, id)
    get(logger.message_limits, id, 1) > 0
end

function Logging.min_enabled_level(logger::ZMQLogger)
    logger.min_level
end

function Logging.catch_exceptions(logger::ZMQLogger)
    false
end

const LEVELS = Dict("info" => Logging.Info, "debug" => Logging.Debug, 
                    "warn" => Logging.Warn, "error" => Logging.Error)
MayBeString = Union{Nothing,String}

function set_global_logger(logfile::MayBeString, level::String = "warn"; quiet::Bool = true, topic = "")
    
    # don't add any line number guff even for debugging
    function quiet_metafmt(level, _module, group, id, file, line)
        color = Logging.default_logcolor(level)
        prefix = (level == Logging.Warn ? "Warning" : string(level)) * ':'
        return color, prefix, ""

    end
    llevel = get(LEVELS, level, Logging.Warn)

    if logfile === nothing
        logger = ConsoleLogger(stderr, llevel, meta_formatter = quiet ? quiet_metafmt : Logging.default_metafmt)
    else
        logger = ZMQLogger(logfile::String, llevel; topic = topic)
    end
    global_logger(logger) 
end