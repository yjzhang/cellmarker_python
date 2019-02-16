#!/bin/bash

# kills all gunicorn processes with 8889

pids=`ps ax | grep gunicorn | grep "8889" | awk '{split($0,a," "); print a[1]}'`
# TODO
for pid in $pids; do
    kill -9 $pid
    echo "killed gunicorn process $pid"
done

