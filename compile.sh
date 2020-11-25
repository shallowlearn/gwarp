# This work was supported by the Intelligence Advanced
# Research Projects Activity (IARPA) via Department of
# Interior / Interior Business Center (DOI/IBC) contract
# number D17PC00280. The U.S. Government is authorized to
# reproduce and distribute reprints for Governmental purposes
# notwithstanding any copyright annotation
# thereon. Disclaimer: The views and conclusions contained
# herein are those of the authors and should not be
# interpreted as necessarily representing the official
# policies or endorsements, either expressed or implied, of
# IARPA, DOI/IBC, or the U.S. Government.

# Author: Bharath Comandur, cjrbharath@gmail.com
# Date: 11/24/2020

set -e

LIB_DIR=$CONDA_PREFIX/lib/

INCLUDE_DIR=$CONDA_PREFIX/include/

if [ ! -e "$LIB_DIR" ]
then
    printf "\nERROR: "$LIB_DIR" not found\n$"
    exit 1
fi

if [ ! -e "$INCLUDE_DIR" ]
then
    printf "\nERROR: "$INCLUDE_DIR" not found\n$"
    exit 1
fi

printf "\nCompiling\n"

command="g++ -std=c++11 -O3 -fopenmp -L $LIB_DIR -I $INCLUDE_DIR c++/gwarp_worker_cpp.cpp -o c++/gwarp++ -lgdal -lboost_program_options -Wl,-rpath=$LIB_DIR"

printf "\n$command\n"

$command

printf "\nDone\n\n"

set +e
