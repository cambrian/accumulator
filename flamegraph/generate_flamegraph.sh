if [ $# -eq 0 ]
  then
      echo "Please supply bench name as argument"
fi
echo "Calling cargo bench in background to generate bench fn executable."
nohup cargo bench --bench $1 >/dev/null 2>&1 &
bench_exec=$(ls -t ../target/release/deps/$1-* | head -1)
echo "Calling dtrace to run this executable and collect its stack frames."
sudo nohup dtrace -x ustackframes=100 -c './'$bench_exec -o out.stacks -n 'profile-997 { @[ustack()] = count(); }' >/dev/null 2>&1
echo "Generating flamegraph."
./stackcollapse.pl out.stacks | grep closure | ./flamegraph.pl > $1.svg
echo "Done. Your flamegraph is now at flamegraph/"$1".svg"
