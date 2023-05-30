#!/usr/bin/env bash

rm -f ./pyJulia/*

RHO="0.5"

# spencers stuff to create the x csv
python3 ./pyFMFM/main.py -x -l pyJulia/

julia --project=. ./julia_code/create_rhs.jl ./pyJulia/x.csv $RHO ./pyJulia/ic.csv

python3 ./pyFMFM/main.py -r -l pyJulia/

julia --project=. ./julia_code/create_matrix.jl ./pyJulia/x.csv $RHO ./pyJulia/sbar.csv ./pyJulia/L.csv

