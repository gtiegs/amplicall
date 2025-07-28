#!/bin/bash
#SBATCH --job-name=clair3
#SBATCH --output=clair3_%A_%a.out
#SBATCH --error=clair3_%A_%a.err
#SBATCH --time=00:30:00
#SBATCH --mem=5G
#SBATCH --cpus-per-task=4        # 4 threads
#SBATCH --array=0-9              # Array for 9 BAM files
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=georgiatiegs@uvic.ca

# Define directories
INPUT_DIR="/home/gtiegs/projects/def-ld172/gtiegs/glcnac/GVLZ7G-IS/align_output_64269857/bam_files/"
OUTPUT_DIR="/home/gtiegs/projects/def-ld172/gtiegs/glcnac/GVLZ7G-IS/clair3_output/"
MODEL_PATH="basecalling_models/rerio/clair3_models/r1041_e82_400bps_sup_v430" 
REF="${INPUT_DIR}/WTamuc1412.fa"

# Define array of BAM files
BAM_FILES=(${INPUT_DIR}/*.sorted.bam)
BAM_FILE=${BAM_FILES[$SLURM_ARRAY_TASK_ID]}
BASENAME=$(basename "$BAM_FILE" .sorted.bam)    # Extracts the base name
OUTPUT_PATH="${OUTPUT_DIR}/${BASENAME}" # Defines the output path

# Run Clair3
singularity exec \
    -B ${INPUT_DIR},${OUTPUT_DIR} \ 
    clair3_latest.sif \
    /opt/bin/run_clair3.sh \
    --bam_fn=${BAM_FILE} \
    --ref_fn=${REF} \
    --threads=${SLURM_CPUS_PER_TASK} \
    --platform="ont" \
    --model_path=${MODEL_PATH} \
    --output=${OUTPUT_PATH} \
    --snp_min_af=0.03 \
    --indel_min_af=0.03 \
    --min_mq=20 \
    --no_phasing_for_fa \
    --include_all_ctgs \
    --haploid_sensitive \
    --var_pct_full=1.0 \
    --var_pct_phasing=1.0 \
    --ref_pct_full=1.0


