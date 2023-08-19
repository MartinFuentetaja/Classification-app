#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=KCNQ3
#SBATCH --cpus-per-task=4
#SBATCH --mem=50gb
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=./tmp/%x-%j.out
#SBATCH --error=./tmp/%x-%j.err
#SBATCH --mail-user=mfuentetaja002@ikasle.ehu.eus
#SBATCH --mail-type=ALL

module load GCCcore/11.3.0
module load Anaconda3/2023.03

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "You can take a coffe while Job is running!"

python gene_app.py -g $SLURM_JOB_NAME -l 872 -o $SLURM_JOB_NAME -cpu $SLURM_CPUS_PER_TASK -config_path "/home/einstein/martin/git_enviroment/Classification_app/config.txt" -column_names "My_Label,pLDDT,MTI_31" -n_splits_grid 14

#My_Label,Initial_A,Initial_C,Initial_E,Initial_F,Initial_G,Initial_H,Initial_K,Initial_L,Initial_M,Initial_N,Initial_P,Initial_Q,Initial_R,Initial_S,Initial_T,Initial_V,Initial_W,Initial_Y,Final_A,Final_C,Final_D,Final_E,Final_F,Final_G,Final_H,Final_I,Final_K,Final_L,Final_M,Final_N,Final_P,Final_Q,Final_R,Final_S,Final_T,Final_V,Final_W,Final_Y,coil,helix,neg_to_neg,neg_to_neu,neg_to_pos,neu_to_neg,neu_to_neu,neu_to_pos,pos_to_neg,pos_to_neu,pos_to_pos,np_to_np,np_to_p,p_to_np,p_to_p,a_to_a,a_to_na,na_to_a,na_to_na,pLDDT,MTI_31,d_size,d_vol,d_pol_e,d_ip_e,d_hf_e,d_msa,residue_conserv
#My_Label,pLDDT,MTI_31,d_size,d_vol,d_pol_e,d_ip_e,d_hf_e,d_msa,residue_conserv
