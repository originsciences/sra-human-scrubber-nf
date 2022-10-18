#!/bin/bash

set -e

for INPUT in input/*; do

    # Write the output to the current working directory
    OUTPUT="\${INPUT#input/}"

    # If the input if gzipped
    if [ "\${INPUT: -3}" == ".gz" ]; then

        echo "Decompressing \$INPUT"

        # Decompress the input, compress the output
        gunzip -c "\$INPUT" \
            | scrub.sh -p ${task.cpus} \
            | gzip -c \
            > "\$OUTPUT"

    else

        echo "Processing \$INPUT"

        # Uncompressed input and output
        cat "\$INPUT" \
            | scrub.sh -p ${task.cpus} \
            > "\$OUTPUT"

    fi

    echo "Done processing \$INPUT"

done