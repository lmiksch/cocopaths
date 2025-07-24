#!/bin/bash

# Define the output queue file
QUEUE_FILE="./queue_file.txt"

# Create or overwrite the queue file with indices from 0 to 421
for i in $(seq 0 598); do
    echo $i >> $QUEUE_FILE
done

echo "Queue file created at $QUEUE_FILE with indices from 0 to 540."