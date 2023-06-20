BEGIN {
    FS="\t"
}

{
    print "chr"$0
}