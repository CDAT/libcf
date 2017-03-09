#!/usr/bin/env bash
cmd="ls"
echo $cmd
$cmd

cmd="pwd"
echo $cmd
$cmd

cmd="export PATH=${HOME}/miniconda/bin:${PATH}"
echo $cmd
$cmd

cmd="conda install -q  -c conda-forge hdf5 libnetcdf lapack clapack curl ossuuid gcc
echo $cmd
$cmd
