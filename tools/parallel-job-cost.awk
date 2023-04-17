BEGIN { 
    OFMT="%.2f"

    minstart = 2000000000; 
    maxstart = 0;

    minend = 2000000000; 
    maxend = 0;

    maxcost = 00000; 
    mincost = 99999; 
} 

$1 ~ /[0-9]+/ { 
    if ($3 < minstart) minstart = $3; 
    if ($3 > maxstart) maxstart = $3; 
    
    cost = $4; 
    endtime = $3 + $4; 

    if (endtime < minend) minend = endtime; 
    if (endtime > maxend) maxend = endtime; 

    if (cost > maxcost) maxcost = cost; 
    if (cost < mincost) mincost = cost; 
} 

END { 
    printf "Min cost : %.2f (s)\n", mincost
    printf "Max cost : %.2f (s)\n", maxcost
    printf "Job start: %s\n", strftime("%Y-%m-%d %H:%M:%S", minstart)
    # printf "Max start: %s\n", strftime("%Y-%m-%d %H:%M:%S", maxstart)
    # printf "Min end  : %s\n", strftime("%Y-%m-%d %H:%M:%S", minend)
    printf "Job end  : %s\n", strftime("%Y-%m-%d %H:%M:%S", maxend)
}