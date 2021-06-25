import Logging
const LOGLEVELS = Dict("info" => Logging.Info, "debug" => Logging.Debug, "warn" => Logging.Warn, 
"error" => Logging.Error)

const ZMQ_WORKER = "tcp://127.0.0.1:9467"
const ZMQ_CLIENT = "ipc:///tmp/chloe-client"
const ZMQ_BACKEND = "ipc:///tmp/chloe-backend"
# change this if you change the API!
const VERSION = "1.3"
