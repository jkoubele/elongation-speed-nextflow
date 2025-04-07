import argparse
import bisect
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import NamedTuple, Optional

import pandas as pd
import pysam
from interval import interval as py_interval
from tqdm import tqdm


@dataclass
class Alignment:
    chromosome: str
    strand: str
    read_1: pysam.AlignedSegment
    read_2: Optional[pysam.AlignedSegment] = None


class Intron(NamedTuple):
    chromosome: str
    strand: str
    start: int
    end: int
    interval: py_interval


class ChromAndStrand(NamedTuple):
    chromosome: str
    strand: str


def ignore_read(read: pysam.AlignedSegment, mapq_threshold=255) -> bool:
    return True if ((read.mapping_quality < mapq_threshold) or
                    read.is_secondary or
                    read.is_supplementary or
                    read.is_unmapped) else False


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_folder',
                        default='/cellfile/datapublic/jkoubele/celegans_mutants/intronic_reads/K002000093_54883')
    parser.add_argument('--output_folder',
                        default='/cellfile/datapublic/jkoubele/celegans_mutants/intronic_reads/K002000093_54883')
    parser.add_argument('--genome_folder',
                        default='/cellfile/datapublic/jkoubele/reference_genomes/WBcel235')
    args = parser.parse_args()

    input_folder = Path(args.input_folder)
    output_folder = Path(args.output_folder)
    genome_folder = Path(args.genome_folder)

    introns_df = pd.read_csv(genome_folder / 'introns.bed', sep='\t',
                             names=['chromosome', 'start', 'end', 'name', 'score', 'strand'])

    print(f"{introns_df.head()=}")

    sample_name = input_folder.name
    output_folder.mkdir(exist_ok=True, parents=True)

    introns_by_chrom_and_strand: dict[ChromAndStrand, list[Intron]] = defaultdict(list)
    for chromosome, strand, start, end in zip(introns_df['chromosome'],
                                              introns_df['strand'],
                                              introns_df['start'],
                                              introns_df['end']):
        intron = Intron(chromosome=chromosome,
                        strand=strand,
                        start=start,
                        end=end,
                        interval=py_interval([start, end]))
        introns_by_chrom_and_strand[ChromAndStrand(chromosome=intron.chromosome,
                                                   strand=intron.strand)].append(intron)

    intron_starts_by_chrom_and_strand = {key: [intron.start for intron in value] for
                                         key, value in introns_by_chrom_and_strand.items()}
    intron_ends_by_chrom_and_strand = {key: [intron.end for intron in value] for
                                       key, value in introns_by_chrom_and_strand.items()}

    # input_bam_path = input_folder / 'Aligned.sortedByCoord.out.bam'
    input_bam_path = input_folder / 'downsampled.bam'
    output_unsorted_bam_path = output_folder / 'intronic_reads_unsorted.bam'
    output_sorted_bam_path = output_folder / 'intronic_reads_sorted.bam'

    print(f"{list(input_folder.iterdir())=}")

    bam_input = pysam.AlignmentFile(input_bam_path, "rb")
    paired_sequencing = True
    strandendess_type = "2"  # either '1' or '2', eligible for paired sequencing only
    overlap_bp_threshold = 5
    create_bam_output = True
    create_json_output = True

    assert strandendess_type in ['1', '2']

    if create_bam_output:
        bam_output = pysam.AlignmentFile(output_unsorted_bam_path, "wb", template=bam_input)

    if paired_sequencing:
        reads_1: dict[str, pysam.AlignedSegment] = {}
        reads_2: dict[str, pysam.AlignedSegment] = {}

        for read in tqdm(bam_input):
            if ignore_read(read):
                continue

            if read.is_read1 and read.query_name not in reads_2:
                reads_1[read.query_name] = read
                continue
            elif read.is_read2 and read.query_name not in reads_1:
                reads_2[read.query_name] = read
                continue

            if read.is_read1 and read.query_name in reads_2:
                read_1 = read
                read_2 = reads_2.pop(read.query_name)
            elif read.is_read2 and read.query_name in reads_1:
                read_1 = reads_1.pop(read.query_name)
                read_2 = read
            else:
                assert False, 'Inconsistent detection of paired reads'

            if read_1.reference_name != read_2.reference_name:
                print(f"{read_1=}")
                print(f"{read_2=}")
                assert False, 'Paired reads aligned to different chromosomes'

            if strandendess_type == '1':
                strand = '+' if read_1.is_forward else '-'
            elif strandendess_type == '2':
                strand = '+' if read_2.is_forward else '-'
            else:
                assert False, 'Invalid strandendess_type'

            alignment = Alignment(chromosome=read_1.reference_name,
                                  strand=strand,
                                  read_1=read_1,
                                  read_2=read_2)

            # Move to function
            intron_list = introns_by_chrom_and_strand.get(ChromAndStrand(alignment.chromosome, alignment.strand))
            if not intron_list:
                continue  # for a case that no introns are present on given contig, e.g. for MtDNA
            aligned_blocks = py_interval(*(
                alignment.read_1.get_blocks() if alignment.read_2 is None else alignment.read_1.get_blocks() + alignment.read_2.get_blocks()))
            intron_list = introns_by_chrom_and_strand[ChromAndStrand(alignment.chromosome, alignment.strand)]
            intron_starts = intron_starts_by_chrom_and_strand[ChromAndStrand(alignment.chromosome, alignment.strand)]
            intron_ends = intron_ends_by_chrom_and_strand[ChromAndStrand(alignment.chromosome, alignment.strand)]

            alignment_start = aligned_blocks[0][0]
            alignment_end = aligned_blocks[-1][1]

            intron_index_lower_bound = bisect.bisect_left(intron_ends, alignment_start)
            intron_index_upper_bound = bisect.bisect_right(intron_starts, alignment_end)

            aligned_to_intron = False
            for intron in intron_list[intron_index_lower_bound:intron_index_upper_bound]:
                interval_intersection = intron.interval & aligned_blocks
                overlap_length = sum([interval[1] - interval[0] for interval in interval_intersection])

                if overlap_length >= overlap_bp_threshold:
                    aligned_to_intron = True
            if aligned_to_intron:
                if create_bam_output:
                    bam_output.write(alignment.read_1)
                    if alignment.read_2 is not None:
                        bam_output.write(alignment.read_2)


    else:
        for read in bam_input:
            alignment = Alignment(chromosome=read.reference_name,
                                  strand='+' if read.is_forward else '-',
                                  read_1=read)

    bam_input.close()
    bam_output.close()

    pysam.sort("-o", str(output_sorted_bam_path), str(output_unsorted_bam_path), catch_stdout=False)
    output_unsorted_bam_path.unlink()
