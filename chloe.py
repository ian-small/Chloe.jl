import zmq
import json
import click

ADDRESS = "tcp://127.0.0.1:9999"
context = zmq.Context()


def setup_zmq(address):
    #  Socket to talk to server
    socket = context.socket(zmq.REQ)  # pylint: disable=no-member
    socket.connect(address)
    return socket


@click.group()
def cli():
    pass


def address(f):
    f = click.option(
        "--address",
        default=ADDRESS,
        help="network address to connect to julia server",
        show_default=True,
    )(f)
    f = click.option(
        "--timeout",
        type=float,
        help="wait timeout in milliseconds [default: wait forever]",
    )(f)
    return f


def poll(socket, timeout):
    poller = zmq.Poller()
    poller.register(socket, zmq.POLLIN)
    events = poller.poll(timeout)
    if not events:
        socket.close(linger=0)
        raise RuntimeError(f"no response after {timeout} milliseconds")


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
    socket = setup_zmq(address)
    msg = dict(cmd="chloe", args=[fasta, output])
    print("sending", msg)
    socket.send_json(msg)
    if timeout is not None:
        poll(socket, timeout)

    resp = socket.recv_json()
    click.secho(f"got {resp}")


@cli.command()
@address
def terminate(timeout, address):
    """Shutdown the server."""
    socket = setup_zmq(address)
    msg = dict(cmd=":terminate")
    socket.send_json(msg)
    if timeout:
        poll(socket, timeout)
    resp = socket.recv_json()
    code = resp["code"]
    click.secho(
        "OK" if code == 200 else f"No Server at {address}",
        fg="green" if code == 200 else "red",
        bold=True,
    )


@cli.command()
@address
def ping(timeout, address):
    """Ping the server."""
    socket = setup_zmq(address)
    msg = dict(cmd="ping")
    socket.send_json(msg)
    if timeout:
        poll(socket, timeout)
    resp = socket.recv_json()
    code = resp["code"]
    click.secho(
        "OK" if code == 200 else f"No Server at {address}",
        fg="green" if code == 200 else "red",
        bold=True,
    )


if __name__ == "__main__":
    cli()

