#!/bin/bash

if [ $(uname -r | grep "fc14") ]; then
  # Fedora 14
  system_release=$(uname -r | awk -F. '{print $5}')
elif [ $(uname -r | grep "fc") ]; then
  # Fedora
  system_release=$(uname -r | awk -F. '{print $4}')
elif [ $(uname -r | grep "el") ]; then
  # RedHat
  #system_release=$(uname -r | awk -F. '{print $4}')
  system_revision=$(lsb_release -r | awk -F' ' '{print $2}')
  #system_release="${system_release}_${system_revision}"
  system_release="el${system_revision}"
else
  # Other
  system_release="noarch"
fi

arch=$(uname -m)
system_name="${system_release}.${arch}"

#echo ${system_release}
echo ${system_name}

