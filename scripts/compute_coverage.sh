usage() {
    echo "Usage: $0 -i <input_folder> -o <output_folder> -f <fai_file>"
    exit 1
}

# Variables to hold argument
input_folder=""
fai_file=""
output_folder=""

# Parse command line argument
while getopts "i:o:f:" opt; do
    case ${opt} in
        i )
            input_folder=$OPTARG
            ;;
        o )
            output_folder=$OPTARG
            ;;
        f )
            fai_file=$OPTARG
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

#input_folder="/cellfile/datapublic/jkoubele/elongation-speed-nextflow/data/intronic_reads/K002000093_54873"
#fai_file="/cellfile/datapublic/jkoubele/reference_genomes/WBcel235/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.fai"
#output_folder="/cellfile/datapublic/jkoubele/elongation-speed-nextflow/data/intronic_coverage/K002000093_54873"

if [ -z "$input_folder" ] || [ -z "$output_folder" ] || [ -z "$fai_file" ]; then
    echo "Error: Missing mandatory arguments"
    usage
fi

# Create output folder if it doesn't exist
mkdir "$output_folder" -p

sort -k 1,1 -k 2,2n $input_folder/intronic_reads_plus_strand.bed > $output_folder/intronic_reads_plus_strand_sorted.bed
bedtools genomecov -bga -split -i $output_folder/intronic_reads_plus_strand_sorted.bed -g $fai_file > \
$output_folder/coverage_plus_strand.bedGraph
gzip -f $output_folder/coverage_plus_strand.bedGraph
gzip -f $output_folder/intronic_reads_plus_strand_sorted.bed

sort -k 1,1 -k 2,2n $input_folder/intronic_reads_minus_strand.bed > $output_folder/intronic_reads_minus_strand_sorted.bed
bedtools genomecov -bga -split -i $output_folder/intronic_reads_minus_strand_sorted.bed -g $fai_file > \
$output_folder/coverage_minus_strand.bedGraph
gzip -f $output_folder/coverage_minus_strand.bedGraph
gzip -f $output_folder/intronic_reads_minus_strand_sorted.bed