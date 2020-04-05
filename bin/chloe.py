#!/usr/bin/env python
from io import StringIO, BytesIO
import os
import gzip
from concurrent.futures import ThreadPoolExecutor
import re

import zmq
import click

PORT = re.compile("^[0-9]+$")
# ADDRESS = "tcp://127.0.0.1:9467"
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


# http://api.zeromq.org/2-1:zmq-setsockopt
def proxy(url_worker, url_client, hwm=1000):
    """Server routine"""
    # pylint: disable=no-member

    # Socket to talk to clients
    clients = context.socket(zmq.ROUTER)
    # number of *messages* in queue
    clients.setsockopt(zmq.RCVHWM, hwm)
    # clients.setsockopt(zmq.ZMQ_SWAP, swap)

    clients.bind(url_client)

    # Socket to talk to workers
    worker = context.socket(zmq.DEALER)
    worker.setsockopt(zmq.SNDHWM, hwm)
    # workers.setsockopt(zmq.ZMQ_SWAP, swap)
    workers.bind(url_worker)

    zmq.proxy(clients, worker)

    # We never get here but clean up anyhow
    clients.close()
    workers.close()
    context.term()


@click.group()
def cli():
    pass


@cli.command()
@click.option(
    "-r", "--remote", type=int, help="remote port to connect to (same as local)"
)
@click.option(
    "-l", "--local", default=9467, help="local port to connect to", show_default=True
)
@click.option("--router", help=f"run a broker endpoint too (e.g. {ADDRESS})")
@click.option("-n", "--nprocs", default=8, help="number of processes to use")
@click.option("--sleep", default=0.0, help="sleep seconds before trying to connect")
@click.option(
    "--level",
    default="info",
    help="log level for server",
    show_default=True,
    type=click.Choice(["info", "warn", "debug", "error"]),
)
def stiletto_ssh(remote, local, nprocs, level, sleep, router):
    """start up a ssh tunnel to stiletto chloe-svr."""
    from threading import Thread
    from time import sleep as Sleep
    from warnings import filterwarnings
    from fabric import Connection
    from cryptography.utils import CryptographyDeprecationWarning

    filterwarnings("ignore", category=CryptographyDeprecationWarning)

    if remote is None:
        remote = local

    remote_dir = "chloe-svr"
    julia = "/home/ianc/julia-1.4.0/bin/julia"

    if router:
        address = f"tcp://127.0.0.1:{local}"
        t = Thread(target=proxy, args=[address, router])
        t.daemon = True
        t.start()

    if sleep:
        Sleep(sleep)
    c = Connection("chloe-stiletto")  # entry in ~/.ssh/config
    # bind on remote and connect on local
    # This means that the broker should be running first
    def run():
        with c.forward_remote(local_port=local, remote_port=remote):
            with c.cd(remote_dir):
                args = f"""-l {level} --nprocs={nprocs} --address=tcp://127.0.0.1:{remote}"""
                cmd = (
                    f"""JULIA_NUM_THREADS=96 {julia} src/chloe_distributed.jl {args}"""
                )
                # need a pty to send interrupt
                c.run(cmd, pty=True)
                click.secho("stiletto terminated", fg="green", bold=True)

    try:
        run()
    except Exception as e:
        print(type(e), e)
        raise


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
            f"{fasta}: {str(data)}", fg="green" if code == 200 else "red", bold=True
        )

    if parallel:
        with ThreadPoolExecutor(max_workers=len(fastas)) as pool:
            for fasta in fastas:
                pool.submit(do_one, fasta)
    else:
        socket = Socket(address, timeout)
        for fasta in fastas:
            do_one(fasta, socket)


def maybegz_open(fasta, mode="rt"):
    if fasta.endswith(".gz"):
        return gzip.open(fasta, mode)
    else:
        return open(fasta, mode)


def gzcompress(b):
    fasta = BytesIO()
    with gzip.GzipFile(fileobj=fasta, mode="wb") as fp:
        fp.write(b)
    return fasta.getvalue()


@cli.command()
@addresses
@click.option("-o", "--output", help="output .sff filename or directory")
@click.option("--binary", is_flag=True, help="don't decompress")
@click.argument("fastas", nargs=-1)
def annotate2(timeout, address, binary, fastas, output):
    """Annotate fasta files (send and receive file content)."""
    socket = Socket(address, timeout)
    for fasta in fastas:
        if binary:
            with open(fasta, "rb") as fp:
                b = fp.read()
            if not fasta.endswith(".gz"):
                b = gzcompress(b)
            fasta = b.decode("latin1")
        else:
            with maybegz_open(fasta) as fp:
                fasta = fp.read()
        code, data = socket.msg(cmd="annotate", args=[fasta])

        if code != 200:
            click.secho(str(data), fg="red", bold=True)
            return

        ncid, sff = data["ncid"], data["sff"]
        click.secho(ncid, fg="green")
        if not output:
            with StringIO(sff) as fp:
                for line in fp:
                    print(line, end="")
        else:
            if os.path.isdir(output):
                tgt = os.path.join(output, f"{ncid}.sff")
            else:
                tgt = output
            with open(tgt, "wt") as fp:
                fp.write(sff)


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
    click.secho(str(data), fg="green" if code == 200 else "red", bold=True)


@cli.command()
@addresses
def workers(timeout, address):
    """Number of service workers"""
    socket = Socket(address, timeout)
    code, data = socket.msg(cmd="nconn")
    click.secho(str(data), fg="green" if code == 200 else "red", bold=True)


if __name__ == "__main__":
    cli()
