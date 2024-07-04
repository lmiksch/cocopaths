#!/bin/bash

#SBATCH --job-name=coco_obj_analysis
#SBATCH --output=logs/coco_obj_analysis_%A_%a.out
#SBATCH --error=logs/coco_obj_analysis_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --cpus-per-task=20
#SBATCH --array=1-40
#SBATCH --time=48:00:00

echo "Starting job $SLURM_ARRAY_TASK_ID"  # Debug statement

QUEUE_FILE="/home/mescalin/miksch/Documents/cocopaths/analysis/slurm/queue_file.txt"
LOCK_FILE="/home/mescalin/miksch/Documents/cocopaths/analysis/slurm/queue_file.lock"
FILE_PATH="/home/mescalin/miksch/Documents/cocopaths/analysis/1_run/6_steps_out.tsv"
OUTPUT_FILE="/home/mescalin/miksch/Documents/cocopaths/analysis/slurm/output/6_steps_out.txt"
CONDA_ENV="bioinf"

CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate $CONDA_ENV

echo "Conda environment activated"  # Debug statement

while true; do
    echo "Acquiring lock on queue file"  # Debug statement

    # Lock the queue file
    exec 200>$LOCK_FILE
    flock -n 200 || exit 1

    echo "Lock acquired"  # Debug statement

    # Read the next index from the queue file
    NEXT_INDEX=$(head -n 1 $QUEUE_FILE)
    echo "Next index: $NEXT_INDEX"  # Debug statement

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
    echo "Lock released"  # Debug statement

    # Run the analysis
    echo "Running analysis for index $NEXT_INDEX"  # Debug statement
    python slurm_obj_test.py --file_path $FILE_PATH --index $NEXT_INDEX --output_file $OUTPUT_FILE

    echo "Analysis complete for index $NEXT_INDEX"  # Debug statement

done

echo "Job completed"  # Debug statement
