#!/usr/bin/bash
#
# set max wallclock time hh:mm:ss
#SBATCH --time=12:00:00
#
# set number of nodes
#SBATCH --nodes=12
#
# set number of tasks
#SBATCH --ntasks=12
#
# set number of cores per node
#SBATCH --cpus-per-task=24
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
echo -e "Interface\tDE.Gene\tDE.Gene.Wilcox.p" > interfacewise_de.tsv
cat slurm_proc_*_de.tsv >> interfacewise_de.tsv
rm slurm_proc_*_de.tsv

# post-processing
Rscript InterfacewiseDEPostprocessor.R

# compress results
zip interfacewise_de.zip interfacewise_de.tsv
rm interfacewise_de.tsv
