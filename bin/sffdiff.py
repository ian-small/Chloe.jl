#!/usr/bin/env python
import os
import click


def diff(fa1, fa2, depth, coverage):
    def todict(l):
        # accD/1/CDS/1	+	59145	1491	0	1.010	0.537	0.999
        name, strand, pos, length, i1, f1, f2, f3, com = l[:-1].split("\t")
        return dict(
            name=name,
            strand=strand,
            pos=int(pos),
            length=int(length),
            phase=int(i1),
            avg=float(f1),
            depth=float(f2),
            coverage=float(f3),
            comment=com,
        )

    lines1 = open(fa1).readlines()
    lines2 = open(fa2).readlines()
    if lines1 == lines2:
        return
    # assert len(lines1) == len(lines2)
    assert lines1[0] == lines2[0]

    d1 = {d["name"]: d for l in lines1[1:] for d in [todict(l)]}
    d2 = {d["name"]: d for l in lines2[1:] for d in [todict(l)]}
    if d1 == d2:
        return
    if d1.keys() != d2.keys():
        s = set(d1) - set(d2)
        if s:
            click.secho(f"new > old {s}", fg="red")
        s = set(d2) - set(d1)
        if s:
            click.secho(f"old > new {s}", fg="red")

    for k in set(d1) & set(d2):
        dd1 = d1[k]
        dd2 = d2[k]
        for n in ["name", "strand", "pos", "length", "phase", "comment"]:
            if not dd1[n] == dd2[n]:
                print(k, n, dd1[n], dd2[n])

        for n in ["avg", "coverage"]:
            f1, f2 = dd1[n], dd2[n]
            if abs(f1 - f2) > coverage:
                print(k, n, f1, f2)

        for n in ["depth"]:
            f1, f2 = dd1[n], dd2[n]
            if abs(f1 - f2) > depth:
                print(k, n, f1, f2)


@click.command()
@click.option("--depth", default=0.2)
@click.option("--coverage", default=0.05)
def run(depth, coverage):
    for sff in os.listdir("testo"):
        if not sff.endswith(".sff"):
            continue

        fa1 = os.path.join("testo", sff)
        fa2 = os.path.join("testfa", sff)
        click.secho(f"diffing {sff}", fg="yellow")
        diff(fa1, fa2, depth=depth, coverage=coverage)


if __name__ == "__main__":
    run()  # pylint: disable=no-value-for-parameter
