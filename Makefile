run-chloe:
	JULIA_NUM_THREADS=8 julia src/chloe_svr.jl -l info --connect --nthreads=3 --address=ipc:///tmp/chloe-worker
