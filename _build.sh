#!/bin/bash

if [ `uname` = "Linux" ]
then
  PROC_NUM=`grep processor /proc/cpuinfo | wc -l`
  if [ $PROC_NUM -gt 4 ]
  then
    PROC_NUM=$(($PROC_NUM/2))
  fi
  OPT="-j $PROC_NUM"
  echo "PROC_NUM = $PROC_NUM"
else
  OPT="WIN32=1"
fi

if ! [ -d ./build/ ]
then
  echo 'no build directory; make dir and run cmake'
  mkdir ./build && cd ./build && cmake -DCMAKE_BUILD_TYPE=Release .. && cd ..
fi

cd ./build && make $OPT && cd ..

#if [ $? -eq 0 ]
#then
#  test `which upx` && upx -q -9 rcrx rctx
#fi

./rfilt_test 2>q

