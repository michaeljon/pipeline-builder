BEGIN {
    FS = "\t"
    xcov = 0;
    ycov = 0;
}

$1 ~ /NC_000023.11/ {
    xcov = $3 / $2
}

$1 ~ /^NC_000024.10/ {
    ycov = $3 / $2
}

END {
    if (ycov > 0) {
        ratio = xcov / ycov
        sex = ratio < 2.0 ? "F" : "M"

        printf "%s\t%s\t%.4f\t%s\n", SEQUENCE, SAMPLE, ratio, sex
    } else {
        printf "%s\t%s\tN/C\tN/C\n", SEQUENCE, SAMPLE
    }
}