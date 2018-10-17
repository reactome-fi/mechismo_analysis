#!/usr/bin/bash
#
# set max wallclock time hh:mm:ss
#SBATCH --time=24:00:00
#
# set number of nodes
#SBATCH --nodes=12
#
# set number of tasks
#SBATCH --ntasks=12
#
# set number of cores per node
#SBATCH --cpus-per-task=32
#
# set output filename
#SBATCH -o slurm-%j.out-%N
#
# mail all alerts (start, end and abortion)
#SBATCH --mail-type=ALL
#
# mailto
#SBATCH --mail-user=burkhajo@ohsu.edu
#
# set job name
#SBATCH --job-name=pinterfacewise_de

# execute on cluster
srun Rscript pInterfacewiseDEAnalyzer.R

# combine results
echo -e "Cancer.Type\tInterface\tDE.Gene\tDE.Gene.Wilcox.p\tNum.Interface.Samples\tNum.NoInterface.Samples\tInterface.NoInterface.Diff" > /home/exacloud/lustre1/WongLab/tmp/pinterfacewise_de.tsv
cat /home/exacloud/lustre1/WongLab/tmp/pslurm_proc_*_de.tsv >> /home/exacloud/lustre1/WongLab/tmp/pinterfacewise_de.tsv
rm /home/exacloud/lustre1/WongLab/tmp/pslurm_proc_*_de.tsv

# post-processing
Rscript pInterfacewiseDEPostprocessor.R

# compress results
zip /home/exacloud/lustre1/WongLab/tmp/pinterfacewise_de.zip /home/exacloud/lustre1/WongLab/tmp/pinterfacewise_de.tsv
rm /home/exacloud/lustre1/WongLab/tmp/pinterfacewise_de.tsv
