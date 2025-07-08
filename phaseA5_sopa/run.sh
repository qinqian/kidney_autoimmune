#!/bin/bash -ex
export -f conda
export -f __conda_activate
export -f __conda_reactivate
export -f __add_sys_prefix_to_path
export -f __conda_hashr
export -f __conda_exe
snakemake -j 12 --use-conda

