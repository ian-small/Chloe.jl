#!/usr/bin/env python
import os
import json
import click

REF = "reference_1116"


@click.command()
@click.option(
    "-r",
    "--refsdir",
    default=REF,
    type=click.Path(file_okay=False),
    help="reference directory",
)
@click.option(
    "-o", "--output", type=click.Path(file_okay=True), help="output json file",
)
def refset(refsdir, output):
    "generate a ReferenceOrganisms.json file"

    if output is None:
        output = os.path.join(refsdir, "ReferenceOrganisms.json")

    def getid(f):
        with open(os.path.join(refsdir, f)) as fp:
            return fp.readline().split()[0][1:]

    fa = [f for f in os.listdir(refsdir) if f.endswith(".fa")]
    ids = [getid(f) for f in fa]
    tgt = [f.split("_")[0] for f in fa]
    # unique names...
    assert len(tgt) == len(set(tgt))
    click.secho(f"found {len(tgt)} organisms", fg="magenta")
    d = dict(zip(ids, tgt))
    with open(output, "w") as fp:
        json.dump(d, fp, sort_keys=True, indent=4)


if __name__ == "__main__":
    refset()  # pylint: disable=no-value-for-parameter
