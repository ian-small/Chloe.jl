#!/usr/bin/env python
import zmq
import click

# url_worker = "tcp://127.0.0.1:9999"
URL_WORKER = "ipc:///tmp/chloe-worker"
# url_client = "inproc://worker"
# url_client = "tcp://127.0.0.1:9998"
URL_CLIENT = "ipc:///tmp/chloe-client"
context = zmq.Context()


def proxy(url_worker, url_client):
    """Server routine"""
    # pylint: disable=no-member

    # Socket to talk to clients
    clients = context.socket(zmq.ROUTER)
    clients.bind(url_client)

    # Socket to talk to workers
    workers = context.socket(zmq.DEALER)
    workers.bind(url_worker)

    zmq.proxy(clients, workers)

    # We never get here but clean up anyhow
    clients.close()
    workers.close()
    context.term()


@click.group()
def cli():
    pass


@cli.command()
@click.option(
    "-a",
    "--address",
    default=URL_WORKER,
    help="address to connect to",
    show_default=True,
)
@click.option("--router", default=URL_CLIENT, help="run a broker endpoint too")
@click.option("-n", "--nprocs", default=4, help="number of processes to use")
@click.option(
    "-l",
    "--level",
    default="info",
    help="log level for server",
    show_default=True,
    type=click.Choice(["info", "warn", "debug", "error"]),
)
def distributed(address, nprocs, level, router):
    """Run the distributed server."""
    from invoke import Context
    from multiprocessing import cpu_count

    if router:
        from threading import Thread

        t = Thread(target=proxy, args=[address, router])
        t.daemon = True
        t.start()
    c = Context()

    def run():

        args = f"""-level={level} --nprocs={nprocs} --address={address}"""
        cmd = (
            f"""JULIA_NUM_THREADS={cpu_count()} julia src/chloe_distributed.jl {args}"""
        )
        # need a pty to send interrupt
        c.run(cmd, pty=True)
        click.secho("chloe terminated", fg="green", bold=True)

    run()


@cli.command()
@click.option(
    "--client",
    default=URL_CLIENT,
    help="network address to connect to chloe server",
    show_default=True,
)
@click.option(
    "--worker",
    default=URL_WORKER,
    help="network address where chloe workers connect to (chloe needs to be run with --connect flag)",
    show_default=True,
)
def broker(worker, client):
    """Run a DEALER/ROUTER ZMQ broker."""
    click.secho(f"workers connect to {worker}, clients connect to {client}", fg="green")
    proxy(worker, client)


if __name__ == "__main__":
    cli()  # pylint: disable=no-value-for-parameter
