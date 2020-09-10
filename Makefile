# run this if there is a DEALER/ROUTER frontend running
# see run-broker or run-chloe-broker:
run-chloe:
	JULIA_NUM_THREADS=8 julia --color=yes src/chloe_distributed.jl -l info --workers=4 --address=tcp://127.0.0.1:9467 --broker=@ipc:///tmp/chloe-client

# start up chloe with a broker in the background
run-chloe-broker:
	JULIA_NUM_THREADS=8 julia --color=yes src/chloe_distributed.jl -l info --workers=4 --broker=ipc:///tmp/chloe-client

# just run the broker
run-broker:
	julia -q --startup-file=no src/broker.jl --worker=tcp://127.0.0.1:9467 --client=ipc:///tmp/chloe-client

run-chloe-logger:
	JULIA_NUM_THREADS=8 julia --color=yes src/chloe_distributed.jl -l info --workers=4 --address=tcp://127.0.0.1:9467 --broker=ipc:///tmp/chloe-client --backend=ipc:///tmp/chloe-logger

run-chloe-backend:
	JULIA_NUM_THREADS=8 julia --color=yes src/chloe_distributed.jl -l info --workers=4 --address=tcp://127.0.0.1:9467 --backend=ipc:///tmp/chloe-backend --broker=@ipc:///tmp/chloe-client

.PHONY: run-chloe run-chloe-broker run-broker run-chloe-logger run-chloe-backend 
