#!/bin/sh

#!/usr/bin/env bash

if [ "$(uname)" == "Darwin" ]; then
  curl -O http://www.flintlib.org/flint-2.5.2.zip
  unzip
elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
  sudo apt-get install libgmp3-dev -y
  sudo apt-get install libflint-dev -y
  g++ -O3 vdf.cpp -lgmpxx -lgmp -lflint -lmpfr -lpthread
fi