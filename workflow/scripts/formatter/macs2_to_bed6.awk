BEGIN {
    FS="\t"
}

NR > 1 {
    print $1 FS $2 FS $3 FS $4 FS $5 FS $6
}