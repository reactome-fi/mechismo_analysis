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
#SBATCH --job-name=interfacewise_de

# execute on cluster
srun Rscript InterfacewiseDEAnalyzer.R

# combine results
echo -e "Cancer.Type\tInterface\tDE.Gene\tDE.Gene.Wilcox.p\tNum.Interface.Samples\tNum.NoInterface.Samples\tInterface.NoInterface.Diff" > /home/users/burkhajo/results/interfacewise_de.tsv
cat /home/users/burkhajo/results/slurm_proc_*_de.tsv >> /home/users/burkhajo/results/interfacewise_de.tsv
rm /home/users/burkhajo/results/slurm_proc_*_de.tsv

# post-processing
Rscript InterfacewiseDEPostprocessor.R

# compress results
zip /home/users/burkhajo/results/interfacewise_de.zip /home/users/burkhajo/results/interfacewise_de.tsv
rm /home/users/burkhajo/results/interfacewise_de.tsv
