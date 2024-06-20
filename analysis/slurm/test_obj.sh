#!/bin/bash

#SBATCH --job-name=coco_obj_analysis
#SBATCH --output=coco_obj_analysis_%A_%a.out
#SBATCH --error=coco_obj_analysis_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --cpus-per-task=10
#SBATCH --array=1-10

QUEUE_FILE="/home/mescalin/miksch/Documents/cocopaths/analysis/slurm/queue_file.txt"
LOCK_FILE="/home/mescalin/miksch/Documents/cocopaths/analysis/slurm/queue_file.lock"
FILE_PATH="/home/mescalin/miksch/Documents/cocopaths/analysis/1_run/3_steps_out.tsv"
OUTPUT_FILE="/home/mescalin/miksch/Documents/cocopaths/analysis/slurm/output/test_out.txt"
CONDA_ENV="bioinf"


CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate $CONDA_ENV



while true; do
    # Lock the queue file
    exec 200>$LOCK_FILE
    flock -n 200 || exit 1

    # Read the next index from the queue file
    NEXT_INDEX=$(head -n 1 $QUEUE_FILE)

    # Check if the queue is empty
    if [ -z "$NEXT_INDEX" ]; then
        echo "Queue is empty. Exiting."
        flock -u 200
        exit 0
    fi

    # Remove the index from the queue file
    tail -n +2 "$QUEUE_FILE" > "$QUEUE_FILE.tmp" && mv "$QUEUE_FILE.tmp" "$QUEUE_FILE"

    # Unlock the queue file
    flock -u 200

    # Run the analysis
    python slurm_obj_test.py --file_path $FILE_PATH --index $NEXT_INDEX --output_file $OUTPUT_FILE

done
