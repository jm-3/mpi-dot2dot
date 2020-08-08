#!/bin/bash

./fio --ioengine=libaio --gtod_reduce=1 --bs=1024K  --size=1G --readwrite=write --thread --numjobs=24 --name=test
