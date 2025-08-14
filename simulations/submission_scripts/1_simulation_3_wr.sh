#!/bin/bash
cd /home/mengbing/research/GSRA_walter/simulations/submission_scripts
#\$SLURM_ARRAY_TASK_ID
seeds=1000
write_slurm() {
    echo "#!/bin/bash
#SBATCH --job-name=wr_$1_$2_N$3
#SBATCH --time=00:10:00
#SBATCH --mail-type=END,FAIL,BEGIN
#SBATCH --mem=3g
#SBATCH --cpus-per-task=1
#SBATCH --array=501-${seeds}
#SBATCH -o ./reports/%x_%A_%a.out
##SBATCH --constraint=E5-2650v4

cd /home/mengbing/research/GSRA_walter/simulations/submission_scripts
module load R

Rscript --verbose ../1_simulation_3_wr.R \$SLURM_ARRAY_TASK_ID $1 $2 $3" > sim_wr_$1_R$2_N$3.slurm
if $4
then
    sbatch sim_wr_$1_R$2_N$3.slurm
fi
}

run=true
Ns=(100 400)

pt_setting="ADependOnZ1Z2"
responder_setting="Hdependent"
for N in "${Ns[@]}"; do
    write_slurm ${pt_setting} ${responder_setting} ${N} ${run}
done

pt_setting="Aconstant"
responder_setting="Z1dependent"
for N in "${Ns[@]}"; do
    write_slurm ${pt_setting} ${responder_setting} ${N} ${run}
done


# for pt_setting in "${pt_settings[@]}"; do
#     for responder_setting in "${responder_settings[@]}"; do
#         for N in "${Ns[@]}"; do
#             write_slurm ${pt_setting} ${responder_setting} ${N} ${run}
#         done
#     done
# done



