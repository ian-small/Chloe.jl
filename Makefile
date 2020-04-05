# run this if there is a DEALER/ROUTER frontend running
# see run-broker:
run-chloe:
	JULIA_NUM_THREADS=8 julia --color=yes src/chloe_distributed.jl -l info --nprocs=4 --address=ipc:///tmp/chloe-worker

run-chloe-9467:
	#JULIA_NUM_THREADS=8 julia src/chloe_svr.jl -l debug --connect --nconn=2 --address=tcp://127.0.0.1:9467
	JULIA_NUM_THREADS=8 julia --color=yes src/chloe_distributed.jl -l info --nprocs=4 --address=tcp://127.0.0.1:9467

run-broker:
	python bin/broker.py broker --worker=ipc:///tmp/chloe-worker --client=ipc:///tmp/chloe-client

run-stiletto:
	JULIA_NUM_THREADS=32 julia src/chloe_svr.jl -l info --nthreads=3 --address=tcp://127.0.0.1:9467 --connect
