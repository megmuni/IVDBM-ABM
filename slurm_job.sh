#!/bin/bash
#SBATCH --account=def-nicoleli
#SBATCH --time=0-0:15:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --gpus-per-node=2
#SBATCH --mem=32000M
#SBATCH --mail-user=meghana.munipalle@mail.mcgill.ca
#SBATCH --mail-type=ALL

module load cuda

chmod a+rwx config_scaffold.txt

#./bin/testRun --numticks 200 --inputfile config_scaffold.txt --wxw 3 --wyw 3 --wzw 3
gdb -ex=r --args ./bin/testRun --numticks 200 --inputfile config_scaffold.txt --wxw 3 --wyw 3 --wzw 3 #debugging with GDB
#valgrind --tool=memcheck --leak-check=yes --show-reachable=yes ./bin/testRun --numticks 100 --inputfile config_scaffold.txt --wxw 0.3 --wyw 0.3 --wzw 0.3