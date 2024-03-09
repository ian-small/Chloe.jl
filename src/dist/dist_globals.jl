import Logging
const LOGLEVELS = Dict("info" => Logging.Info, "debug" => Logging.Debug, "warn" => Logging.Warn,
    "error" => Logging.Error)

const ZMQ_WORKER = "tcp://127.0.0.1:9458"
const ZMQ_ENDPOINT = "ipc:///tmp/chloe5-client"
const ZMQ_BACKEND = "ipc:///tmp/chloe5-backend"
# change this if you change the API!
const VERSION = "5.0"
