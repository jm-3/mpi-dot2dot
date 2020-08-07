#!/bin/bash

SIZE=10240

for i in 1 2 4 8 16 32; do
        echo "${i}MB"
        mkdir ${i}MB
        lfs setstripe -S ${i}M -i -1 -c -1 ${i}MB
        cd ${i}MB
        echo "Writing"
        dd if=/dev/zero of=stripefile bs=${i}M count=$((SIZE/$i)) oflag=direct status=progress
        lfs getstripe stripefile
        echo "Reading"
        dd if=stripefile of=/dev/null bs=${i}M iflag=direct status=progress
        cd ..
done
