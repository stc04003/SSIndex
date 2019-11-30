#!/bin/bash

#SBATCH -J GSM
#SBATCH -N NODE                              # number of nodes
#SBATCH --ntasks-per-node CORE               # number of cores to evolve for each node
#SBATCH -p normal                            # Partition
#SBATCH -t 38:58:00
#SBATCH --mail-user=syhan.chiou@gmail.com
#SBATCH --mail-type=end

work=`pwd`

module restore gnuR
module load launcher

export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
#export LAUNCHER_WORKDIR=$work
export LAUNCHER_SCHED=interleaved
export LAUNCHER_JOB_FILE=SIMJOB
$LAUNCHER_DIR/paramrun
