#!/bin/bash

thre=0.95
rm yuu-${thre}-summary
for f in *yuu-${thre}.*stats.o
do 
    gawk -f parse.awk $f >> yuu-${thre}-summary
done
exit
rm mahito-summary
for f in *mahito.*stats.o
do 
    gawk -f parse.awk $f >> mahito-summary
done
