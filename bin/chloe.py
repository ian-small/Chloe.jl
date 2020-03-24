#!/usr/bin/env python
import zmq
import json
import click

# ADDRESS = "tcp://127.0.0.1:9999"
ADDRESS = "ipc:///tmp/chloe-client"
context = zmq.Context()


class Socket:
    def __init__(self, address, timeout=None):
        self.socket, self.poller = setup_zmq(address)
        self.timeout = timeout

    def poll(self):
        if self.timeout:
            events = self.poller.poll(self.timeout)
            if not events:
                self.socket.close(linger=0)
                self.poller.unregister(self.socket)
                raise RuntimeError(f"no response after {self.timeout} milliseconds")

    def msg(self, **kwargs):
        self.socket.send_json(kwargs)
        self.poll()
        resp = self.socket.recv_json()
        return resp["code"], resp["data"]


def setup_zmq(address):
    #  Socket to talk to server
    socket = context.socket(zmq.REQ)  # pylint: disable=no-member
    socket.connect(address)
    poller = zmq.Poller()
    poller.register(socket, zmq.POLLIN)
    return socket, poller


@click.group()
def cli():
    pass


def address(f):
    f = click.option(
        "-a",
        "--address",
        default=ADDRESS,
        help="network address to connect to julia server",
        show_default=True,
    )(f)
    f = click.option(
        "-t",
        "--timeout",
        type=float,
        help="wait timeout in milliseconds [default: wait forever]",
    )(f)
    return f


@cli.command()
@address
@click.option(
    "-o",
    "--output",
    required=True,
    help="output .sff filename (relative to the server)",
)
@click.argument("fasta")
def annotate(timeout, address, fasta, output):
    """Annotate a fasta file."""
    socket = Socket(address, timeout)
    code, data = socket.msg(cmd="chloe", args=[fasta, output])

    click.secho(
        str(data) if code == 200 else f"No Server at {address}",
        fg="green" if code == 200 else "red",
        bold=True,
    )
    print("returned data", data)


def num_threads(socket):
    _, data = socket.msg(cmd="threads")
    return data


@cli.command()
@click.option("-n", "--nthreads", default=0)
@address
def terminate(timeout, address, nthreads):
    """Shutdown the server."""
    socket = Socket(address, timeout)
    for _ in range(nthreads or num_threads(socket)):
        code, _ = socket.msg(cmd=":terminate")
        click.secho(
            "OK" if code == 200 else f"No Server at {address}",
            fg="green" if code == 200 else "red",
            bold=True,
        )


@cli.command()
@address
def ping(timeout, address):
    """Ping the server."""
    socket = Socket(address, timeout)
    code, data = socket.msg(cmd="ping")
    click.secho(
        str(data) if code == 200 else f"No Server at {address}",
        fg="green" if code == 200 else "red",
        bold=True,
    )


@cli.command()
@address
def workers(timeout, address):
    """Number of service workers"""
    socket = Socket(address, timeout)
    code, data = socket.msg(cmd="threads")
    click.secho(
        str(data) if code == 200 else f"No Server at {address}",
        fg="green" if code == 200 else "red",
        bold=True,
    )


if __name__ == "__main__":
    cli()
