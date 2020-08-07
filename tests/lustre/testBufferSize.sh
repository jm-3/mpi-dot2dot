#!/bin/bash

SIZE=10240

mkdir S16MB
lfs setstripe -S 16M -i -1 -c -1 S16MB
cd S16MB

for i in 1 2 4 8 16 32; do
        echo "${i}MB"
        echo "Writing"
        dd if=/dev/zero of=stripefile bs=${i}M count=$((SIZE/$i)) oflag=direct status=progress
        lfs getstripe stripefile
        echo "Reading"
        dd if=stripefile of=/dev/null bs=${i}M iflag=direct status=progress
        rm stripefile
done
