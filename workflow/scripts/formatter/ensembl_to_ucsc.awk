BEGIN {
    FS="\t"
}

{
    {
        if ($3!="*") $3="chr"$3
    }

    {
        if ($7!="=") $7="chr"$7
    }

    print $0
}
