#!/bin/sh

# Set the number of nodes and processes per node
#PBS -l nodes=1:ppn=1

# Set max wallclock time
#PBS -l walltime=48:00:00

# Set maximum memory
#PBS -l mem=16gb

# Set name of job
#PBS -N Job1

# Use submission environment
#PBS -V
cd ~/PlasmoAlgorithms/examples/SIPLIB
julia sslp_example.jl
