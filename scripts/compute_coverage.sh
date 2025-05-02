usage() {
    echo "Usage: $0 -b <bed_file> -f <fai_file> -o <output_folder>"
    exit 1
}

# Variables to hold argument
bed_file=""
fai_file=""
output_folder=""

# Parse command line argument
while getopts ":f:" opt; do
    case ${opt} in
        f )
            fai_file_name=$OPTARG
            ;;
        f )
            fai_file_name=$OPTARG
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

bed_file="/home/jakub/Desktop/elongation-speed-nextflow/data/intronic_reads/K002000093_54873/intronic_reads_plus_strand.bed"
fai_file="/home/jakub/Desktop/elongation-speed-nextflow/reference_genomes/WBcel235/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.fai"
output_folder="/home/jakub/Desktop/elongation-speed-nextflow/data/intronic_reads/K002000093_54873"

filename=$(basename "$bed_file")
filename="${filename%.*}"
echo "$filename"

# Check if mandatory argument is provided
if [ -z "$fai_file" ]; then
    echo "Error: Missing mandatory argument"
    usage
fi


sort -k 1,1 -k 2,2n $bed_file > $output_folder/"${filename}"_sorted.bed
bedtools genomecov -bga -split -i $output_folder/"${filename}"_sorted.bed -g $fai_file > \
/$output_folder/coverage_"${filename}".bedGraph
gzip $output_folder/"${filename}"_sorted.bed
