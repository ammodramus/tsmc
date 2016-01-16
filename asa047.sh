#!/bin/bash
#
cp asa047.h /$HOME/include
#
gcc -c asa047.c
if [ $? -ne 0 ]; then
  echo "Errors compiling asa047.c."
  exit
fi
#
mv asa047.o ~/libc/$ARCH/asa047.o
#
echo "Library installed as ~/libc/$ARCH/asa047.o"
