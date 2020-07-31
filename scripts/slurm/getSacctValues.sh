#!/bin/bash

export SACCT_FORMAT="JobID%20,JobName,User,Partition,NodeList,Elapsed,ElapsedRaw,State,ExitCode,MaxRSS,AllocTRES%32"
files=`find ~/TFM/resultados/ -name dot.*.out`

for i in $files; do
        dir=`dirname $i`
        jobid=`echo $i | sed -e 's/.*dot\.\(.*\).out/\1/'`
        if [ ! -f "$dir/sacct.txt" ]; then
                echo $dir
                sacct -j $jobid > $dir/sacct.txt
        fi
done

