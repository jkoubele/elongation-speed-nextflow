#!/bin/bash

#SBATCH --job-name=compute_intronic_coverage
#SBATCH --ntasks=5
#SBATCH --mem=15G

# Function to display usage information
usage() {
    echo "Usage: $0 -i <input_folder> -o <output_folder> -d <docker_image_path>"
    echo "-s <script_folder> -g <genome_folder> -f <fai_file_name>"
    exit 1
}

# Variables to hold arguments
input_folder=""
output_folder=""
docker_image_path=""
genome_folder=""
script_folder=""
fai_file_name=""

# Parse command line arguments
while getopts ":i:o:d:g:s:f:" opt; do
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
        f )
            fai_file_name=$OPTARG
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
[ -z "$genome_folder" ] || [ -z "$script_folder" ] || [ -z "$fai_file_name" ]; then
    echo "Error: Missing mandatory arguments"
    usage
fi

# Check if the docker image is available, and load it from disk if it's not
if ! docker images --format "{{.Repository}}" | grep -q "^bioinfo_tools$"; then
    docker load -i "$docker_image_path"
fi

# Create output folder if it doesn't exist
mkdir "$output_folder" -p

# Run docker with script extracting the intronic reads

docker run --rm \
-v "$input_folder":/input_folder \
-v "$output_folder":/output_folder \
-v "$genome_folder":/genome_folder \
-v "$script_folder":/script_folder \
--security-opt seccomp=unconfined \
bioinfo_tools /bin/sh -c "sh /script_folder/compute_coverage.sh \
-i /input_folder \
-o /output_folder \
-f /genome_folder/$fai_file_name; \
chmod 777 -R /output_folder"