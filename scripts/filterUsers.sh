#!/bin/sh
# first argument = filePath
# 2nd argument = useId
cat $1| grep -e "^$2\t" | tee $1_User$2_data
