BEGIN{
    FS="\t";
    print "Sample" FS "Reads_in_Peak" FS "Mapped_Reads" FS "FRiP"
}

{
    print $1 FS $2 FS $4 FS $2 / $4
}