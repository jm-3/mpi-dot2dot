#!/bin/bash

if [ ! $# -eq 2 ]; then
        echo "Forma de uso: $0 dataset algoritmo"
        exit 1
fi

DATASET=$1
ALGNAME=$2

ORIGFILE="$HOME/TFM/resultados/$DATASET/Dot-1.0.p3/threads-1/test-1/output.dot"


for i in `ls -1 ~/TFM/resultados/$DATASET/$ALGNAME/`; do
        for j in `ls -1 ~/TFM/resultados/$DATASET/$ALGNAME/$i`; do
                CHECKFILE="$HOME/TFM/resultados/$DATASET/$ALGNAME/$i/$j/$k/output.dot"
                OKFILE="$HOME/TFM/resultados/$DATASET/$ALGNAME/$i/$j/$k/ok-output.dot"
                if [ ! -f $OKFILE ]; then
                        diff -q $ORIGFILE $CHECKFILE
                        if [ "$?" -eq "0" ]; then
                                touch $OKFILE
                                rm -f $CHECKFILE
                        fi
                fi
        done
done
