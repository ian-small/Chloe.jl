import sys
from datetime import datetime
from julia import Julia
jl = Julia('/home/ianc/src/julia-1.3.1/bin/julia', compiled_modules=False)
from julia import Main

S = """
include("main.jl")
main(fasta_files=["../../chloe/xxx.fa"])
"""

def run():
    s = datetime.now()
    Main.eval(S)
    e = datetime.now()
    print(e - s, file=sys.stderr)

run()
run()
run()

