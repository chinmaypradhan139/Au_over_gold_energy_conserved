#!/usr/bin/env bash
rm a.out
ifort gaussquad.f90
./a.out
rm -r running
mkdir running

ifort mod_atom_interaction.f90 Tulli13_1990.f90 main.f90 -ipo -O3 -no-prec-div -static-intel -fp-model fast=2 -xHost -I${MKLROOT}/include/intel64/lp64 -I${MKLROOT}/include${MKLROOT}/lib/intel64/libmkl_blas95_lp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl

cp a.out running
cp fort.23 running
cp parallel_script_v4.py running
cp raw_x.txt running
cp raw_w.txt running
cp job.sh running
cp 528atom.dat running
cp gold_neigh_pos.dat running
cp gold_positions.inp running
cp gold_velocities.inp running
cp q_number.dat running
cd running

python3 parallel_script_v4.py 



