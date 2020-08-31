#!/usr/bin/env python
import zmq
from zmq.devices import ThreadDevice
import click

URL_WORKER = "tcp://127.0.0.1:9467"
# URL_WORKER = "ipc:///tmp/chloe-worker"
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


def device(url_worker, url_client):
    # pylint: disable=no-member
    td = ThreadDevice(zmq.QUEUE, in_type=zmq.ROUTER, out_type=zmq.DEALER)
    td.bind_in(url_client)
    td.bind_out(url_worker)
    td.setsockopt_in(zmq.IDENTITY, 'ROUTER')
    td.setsockopt_out(zmq.IDENTITY, 'DEALER')
    td.start()


@click.group()
def cli():
    pass


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
