#!/usr/bin/env python
import os
from os.path import isfile
import click

yellow = lambda s: click.style(s, fg="yellow")
red = lambda s: click.style(s, fg="red", bold="true")


def diff(fa1, fa2, depth, coverage, skip_comments):
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
            comment=com.strip(),
        )

    def ps(key, n, f1, f2):
        return f'{prefix}{key} [{yellow(n)}]: "{f1}" {yellow("!=")} "{f2}" [ref]'

    def fixcom(s):
        s = s.replace("possible pseudogene", "")
        s = s.replace(",", "")
        s = s.split()
        return " ".join(s)

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
    prefix = "\t"
    if d1.keys() != d2.keys():
        s = set(d1) - set(d2)
        if s:
            for t in sorted(s):
                click.echo(red(f"{prefix}new: {t}"))
        s = set(d2) - set(d1)
        if s:
            for t in sorted(s):
                click.echo(red(f"{prefix}old: {t}"))

    fields = ["name", "strand", "pos", "length", "phase"]

    def check_pos(dd1, dd2):
        p1, l1 = dd1["pos"], dd1["length"]
        p2, l2 = dd2["pos"], dd2["length"]

        if p1 != p2 or l1 != l2:
            ll = []
            if p1 - p2 == l1 - l2:
                ll.append(f"gene shifted {p1 -p2}")
            if p1 + l1 == p2 + l2:
                ll.append(f"start moved {p1 - p2}")
            if p1 == p2 and l1 != l2:
                ll.append(f"end moved  {l1 -l2}")
            if ll:
                click.secho(f"{prefix}{k} {', '.join(ll)}")

    for k in set(d1) & set(d2):
        dd1 = d1[k]
        dd2 = d2[k]

        check_pos(dd1, dd2)

        for n in fields:
            f1, f2 = dd1[n], dd2[n]
            if not f1 == f2:
                click.secho(ps(k, n, f1, f2))

        if not skip_comments:
            f1, f2 = dd1["comment"], dd2["comment"]
            if not fixcom(f1) == fixcom(f2):
                click.secho(ps(k, "comment", f1, f2))

        for n in ["avg", "coverage"]:
            f1, f2 = dd1[n], dd2[n]
            if abs(f1 - f2) > coverage:
                click.secho(ps(k, n, f1, f2))

        for n in ["depth"]:
            f1, f2 = dd1[n], dd2[n]
            if abs(f1 - f2) > depth:
                click.secho(ps(k, n, f1, f2))


@click.command()
@click.option("--depth", default=0.2, show_default=True, help='depth "closeness"')
@click.option(
    "--coverage", default=0.05, show_default=True, help='coverage "closeness"'
)
@click.option(
    "--src",
    default="testfa",
    type=click.Path(file_okay=False, dir_okay=True),
    help="source directory for .sff files",
)
@click.option("-c", "--skip-comments", is_flag=True, help="ignore comment differences")
def run(depth, coverage, skip_comments, src):

    for sff in os.listdir("testo"):
        if not sff.endswith(".sff"):
            continue

        fa1 = os.path.join("testo", sff)
        fa2 = os.path.join(src, sff)
        if not isfile(fa2):
            click.secho(f"can't find ref file: {fa2}", fg="red")
            continue
        click.echo(yellow(f"diffing {sff}"))
        diff(fa1, fa2, depth=depth, coverage=coverage, skip_comments=skip_comments)


if __name__ == "__main__":
    run()  # pylint: disable=no-value-for-parameter
