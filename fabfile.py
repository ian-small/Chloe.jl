# import sys
import os
from fabric import task
from click import style as color, echo, secho


local_dir = os.path.dirname(__file__) or "."
remote_dir = "/var/www/websites3/annotator"

HOSTS = ["ianc@croppal"]


def git_uptodate(res):
    # Already up-to-date or Already up to date
    # is there a better way?
    return res.stdout.lower().startswith("already up")


@task(hosts=HOSTS)
def update(c):
    """Update code from github and restart (if any changes)."""
    with c.cd(remote_dir):
        pwd = c.run("pwd", hide=True)
        echo("pulling from: " + color(pwd.stdout.strip(), fg="yellow"))
        result = c.run("uname -a", hide=True)
        # ,result.failed,result.return_code,result.succeeded
        echo(color(result.stdout.strip(), fg="green"))
        res = c.run("git pull", warn=True)
        if not res.failed and not git_uptodate(res):
            secho('restarting service use: "fab ping" to check restart', fg="magenta")
            c.run("python bin/chloe.py terminate")


@task(hosts=HOSTS)
def restart(c):
    """Restart annotator service."""
    with c.cd(remote_dir):
        c.run("python bin/chloe.py terminate")


@task(hosts=HOSTS)
def ping(c):
    """Ping annotator service."""
    with c.cd(remote_dir):
        c.run("python bin/chloe.py ping")
