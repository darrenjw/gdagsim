#!/bin/sh

SPACE=20
TIME=30
ITERS=5

ulimit -d
ulimit -s

echo $SPACE $SPACE $TIME $ITERS
date
cd /home/ndjw1/src/gdag/Examples
sptimeprior $SPACE $SPACE $TIME
date
time sptime $SPACE $SPACE $TIME $ITERS  > sptime.err 2>&1
date

