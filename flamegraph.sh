#!/bin/sh

FLAMEGRAPH_FILENAME=$1.svg
FLAMEGRAPH_DIR=./flamegraph/graphs
STACKS_FILENAME=$1.stacks
STACKS_DIR=./flamegraph/stacks

if [ $# -eq 0 ]
  then
    echo "Please supply bench name as argument."
fi

mkdir -p $FLAMEGRAPH_DIR
mkdir -p $STACKS_DIR

echo "Calling cargo bench in background to generate bench fn executable."
nohup cargo bench --bench $1 >/dev/null 2>&1 &
bench_exec=$(ls -t ./target/release/deps/$1-* | head -1)
echo "Calling dtrace to run this executable and collect its stack frames."
sudo nohup dtrace -x ustackframes=100 -c './'$bench_exec -o $STACKS_DIR/$STACKS_FILENAME -n 'profile-997 { @[ustack()] = count(); }' >/dev/null 2>&1
echo "Generating flamegraph."
./flamegraph/stackcollapse.pl $STACKS_DIR/$STACKS_FILENAME | grep closure | ./flamegraph/flamegraph.pl > $FLAMEGRAPH_DIR/$FLAMEGRAPH_FILENAME
echo "Done. Your flamegraph is now at: "$FLAMEGRAPH_DIR/$FLAMEGRAPH_FILENAME
