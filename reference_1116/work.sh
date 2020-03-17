ls *.gbff|awk '{print "perl gb2sff.pl "$1" >"$1".to.sff"}'|sh
ls *.gbff|awk '{print "perl gb2sff_intron.pl "$1" >"$1".to.sff.intron"}' |sh
ls *.fa|awk '{print "julia iansSuffixArrays.jl "$1}'|sh

ls *.gbff|awk '{print "perl gb_to_fa.pl "$1" >"$1".fa"}'|sh
ls *.gbff|awk '{print "perl gb2_new_format_v2.pl "$1" >"$1".to.sff.tmp"}'|sh
ls *.gbff|awk '{print "perl insert_intron.pl "$1".to.sff.tmp "$1".to.sff.intron > "$1".to.sff"}'|SH
