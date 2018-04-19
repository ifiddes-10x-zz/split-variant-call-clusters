#!/usr/bin/env bash
set -e
deck=$1
d=$PWD
pipeline=`basename $d`
cd /mnt/deck${deck}/ian
mkdir -p ${pipeline}
cd ${pipeline}
ln -sf ${d}/* ./