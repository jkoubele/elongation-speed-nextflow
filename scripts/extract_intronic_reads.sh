#!/bin/bash

#SBATCH --job-name=extract_intronic_reads
#SBATCH --ntasks=3

# Function to display usage information
usage() {
    echo "Usage: $0 -i <input_folder> -o <output_folder> -d <docker_image_path> -g <genome_folder>"
    echo "[-s <script_folder>]"
    exit 1
}

# Variables to hold arguments
input_folder=""
output_folder=""
docker_image_path=""
genome_folder=""
script_folder=""

# Parse command line arguments
while getopts ":i:o:d:g:s:" opt; do
    case ${opt} in
        i )
            input_folder=$OPTARG
            ;;
        o )
            output_folder=$OPTARG
            ;;
        d )
            docker_image_path=$OPTARG
            ;;
        g )
            genome_folder=$OPTARG
            ;;
        s )
            script_folder=$OPTARG
            ;;
        \? )
            echo "Invalid option: $OPTARG" 1>&2
            usage
            ;;
        : )
            echo "Invalid option: $OPTARG requires an argument" 1>&2
            usage
            ;;
    esac
done

# Check if mandatory arguments are provided
if [ -z "$input_folder" ] || [ -z "$output_folder" ] || [ -z "$docker_image_path" ] || \
[ -z "$genome_folder" ] || [ -z "$script_folder" ]; then
    echo "Error: Missing mandatory arguments"
    usage
fi

# Check if the docker image is available, and load it from disk if it's not
if ! docker images --format "{{.Repository}}" | grep -q "^bioinfo_tools$"; then
    docker load -i "$docker_image_path"
fi

echo "SHELL DOCKER WRAPPER DEBUG: Output folder=$output_folder"$
echo "SHELL DOCKER WRAPPER DEBUG: Input folder=$input_folder"$
ls "$output_folder"
ls "$input_folder"
# Create output folder if it doesn't exist
mkdir "$output_folder" -p

# Run docker with script extracting the intronic reads
echo "DEBUG: Running docker"
docker run --rm \
-v "$(realpath "$input_folder")":/input_folder \
-v "$output_folder":/output_folder \
-v "$genome_folder":/genome_folder \
-v "$script_folder":/script_folder \
--security-opt seccomp=unconfined \
bioinfo_tools /bin/sh -c "ls /; ls /input_folder; python3 /script_folder/extract_intronic_reads.py \
--input_folder /input_folder \
--output_folder /output_folder \
--genome_folder /genome_folder; \
chmod 777 -R /output_folder"