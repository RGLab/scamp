#!/bin/bash
sbatch <<EOT
#!/bin/bash
#SBATCH --partition=partitionName
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=100:00:00
#SBATCH --mem=100000
#SBATCH -o ./logs/sepMixture_scenario_1_part_1_setting_$(($1))$(($2))$(($3))$(($4))
#SBATCH -J mix$(($1))$(($2))$(($3))$(($4))
#SBATCH --threads-per-core=1
#SBATCH --dependency=singleton

echo "Start of program at `date`" 
Rscript --no-save --no-restore /path/to/1_separatedMixtureSimulation.R -a $1 -b $2 -c $3 -d $4
echo "End of program at `date`"
EOT
