#!/usr/bin/env python
from io import StringIO
from concurrent.futures import ThreadPoolExecutor
import re
import zmq
import click

PORT = re.compile("^[0-9]+$")
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


def addresses(f):
    def callback(ctx, param, value):
        if PORT.match(value):
            return f"tcp://127.0.0.1:{value}"
        return value

    f = click.option(
        "-a",
        "--address",
        default=ADDRESS,
        help="network address to connect to julia server",
        show_default=True,
        callback=callback,
    )(f)
    f = click.option(
        "-t",
        "--timeout",
        type=float,
        help="wait timeout in milliseconds [default: wait forever]",
    )(f)
    return f


@cli.command()
@addresses
@click.option("--parallel", is_flag=True, help="send annotation requests in parallel")
@click.option(
    "-o",
    "--output",
    required=True,
    help="output .sff filename or directory (relative to the server)",
)
@click.argument("fastas", nargs=-1)
def annotate(timeout, address, fastas, output, parallel):
    """Annotate fasta files."""

    def do_one(fasta, socket=None):
        if socket is None:
            socket = Socket(address, timeout)

        code, data = socket.msg(cmd="chloe", args=[fasta, output])

        click.secho(
            f"{fasta}: {str(data)}"
            if code == 200
            else f"{fasta}: No Server at {address}",
            fg="green" if code == 200 else "red",
            bold=True,
        )

    if parallel:
        with ThreadPoolExecutor(max_workers=len(fastas)) as pool:
            for fasta in fastas:
                pool.submit(do_one, fasta)
    else:
        socket = Socket(address, timeout)
        for fasta in fastas:
            do_one(fasta, socket)


@cli.command()
@addresses
@click.argument("fastas", nargs=-1)
def annotate2(timeout, address, fastas):
    """Annotate fasta files (send and receive file content)."""
    socket = Socket(address, timeout)
    for fasta in fastas:
        with open(fasta) as fp:
            fasta = fp.read()
        code, data = socket.msg(cmd="annotate", args=[fasta])

        if code != 200:
            click.secho(data, fg="red", bold=True)
            return

        ncid, sff = data["ncid"], data["sff"]
        click.secho(ncid, fg="green")
        with StringIO(sff) as fp:
            for line in fp:
                print(line, end="")


def num_conn(socket):
    _, data = socket.msg(cmd="nconn")
    return data


@cli.command()
@click.option("-n", "--nconn", default=0)
@addresses
def terminate(timeout, address, nconn):
    """Shutdown the server."""
    socket = Socket(address, timeout)
    # terminate each thread.
    nconn = nconn or num_conn(socket)
    click.secho(f"terminating {nconn} server(s) @ {address}", fg="magenta")
    for _ in range(nconn):
        code, _ = socket.msg(cmd=":terminate")
        click.secho(
            "OK" if code == 200 else f"No Server at {address}",
            fg="green" if code == 200 else "red",
            bold=True,
        )


@cli.command()
@addresses
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
@addresses
def workers(timeout, address):
    """Number of service workers"""
    socket = Socket(address, timeout)
    code, data = socket.msg(cmd="nconn")
    click.secho(
        str(data) if code == 200 else f"No Server at {address}",
        fg="green" if code == 200 else "red",
        bold=True,
    )


if __name__ == "__main__":
    cli()
