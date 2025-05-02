## Celegans mutants
 * Extract intronic reads: ```sh batch_extract_intronic_reads.sh -i /home/jkoubele/celegans_mutants/BAM -o /home/jkoubele/celegans_mutants/intronic_reads -g /home/jkoubele/reference_genomes/WBcel235```
 * Compute coverage: ```sh batch_compute_coverage.sh -i /home/jkoubele/celegans_mutants/intronic_reads -o /home/jkoubele/celegans_mutants/intronic_coverage -g /home/jkoubele/reference_genomes/WBcel235 -f Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.fai```
## Drosophila mutants

 * Extract intronic reads: ```sh batch_extract_intronic_reads.sh -i /home/jkoubele/drosophila_mutants/BAM -o /home/jkoubele/drosophila_mutants/intronic_reads -g /home/jkoubele/reference_genomes/BDGP6.46```