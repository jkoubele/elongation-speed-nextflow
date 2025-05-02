#!/bin/bash

# Function to display usage information
usage() {
    echo "Usage: $0 -i <input_folder> -o <output_folder> -g <genome_folder> -f <fai_file_name>"
    exit 1
}

script_directory="$(cd "$(dirname "$0")" && pwd)"
repository_path="$(dirname "$script_directory")"

echo $repository_path

# Variables to hold arguments
input_folder=""
output_folder=""
genome_folder=""
fai_file_name=""
docker_image_path="$repository_path"/docker_images/bioinfo_tools.tar
slurm_log_folder="$repository_path"/slurm_logs

# Parse command line arguments
while getopts ":i:o:g:f:" opt; do
    case ${opt} in
        i )
            input_folder=$OPTARG
            ;;
        o )
            output_folder=$OPTARG
            ;;
        g )
            genome_folder=$OPTARG
            ;;
        f )
            fai_file_name=$OPTARG
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
if [ -z "$input_folder" ] || [ -z "$output_folder" ] || [ -z "$genome_folder" ] || [ -z "$fai_file_name" ]; then
    echo "Error: Missing mandatory arguments"
    usage
fi

# Create output folder if it doesn't exist
mkdir "$output_folder" -p

for sub_folder in "$input_folder"/*; do
  sample_name=$(basename "$sub_folder")
  echo "Submitting sample $sample_name"
  sbatch --output="$slurm_log_folder"/%j_%x.log --error="$slurm_log_folder"/%j_%x.err \
  "$repository_path"/scripts/compute_coverage_docker_wrapper.sh \
  -i "$sub_folder" \
  -o "$output_folder"/"$sample_name" \
  -d "$docker_image_path" \
  -g "$genome_folder" \
  -s "$repository_path"/scripts \
  -f "$fai_file_name"
done