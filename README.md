
# Ambimux

Software tool for demultiplexing single cell multiome data.

## Installation

The source code can be downloaded using
```
git clone --recurse-submodules git@github.com:marcalva/ambimux.git
```
The `--recurse-submodules` flag makes sure the `htslib` dependency is downloaded.
If not, you can also clone this manually.

After cloning, you can build htslib and then ambimux using the following
sequence of commands. `make hts` will download htslib and build the library,
while make will build ambimux
```bash
make hts
make
```

The only dependency is htslib, which is included as a subdirectory. Note that
ambimux requires htslib to be built successfully during the `make` step.
Thus all dependencies for htslib are also required for ambimux.

## Usage

An example usage of ambimux is below.
```
ambimux \
    --atac-bam "$atac_bam" \
    --rna-bam "$gex_bam" \
    --vcf "$vcf" \
    --gtf "$gtf" \
    --peaks "$peaks" \
    --samples "$samples" \
    --bc-wl "$bcsd" \
    --rna-mapq 30 \
    --atac-mapq 30 \
    --tx-basic \
    --eps 1-e6 \
    --max-iter 100 \
    --seed 123 \
    --threads 8 \
    --verbose \
    --out path/to/output/dir/ambimux

    --exclude $excl \
    --out-min 100 \
```
where the variables are replaced by your own files.

For a detailed explanation of each parameter, see below.

## Alignment data

Ambimux can use data from single cell ATAC, RNA, or both. The ATAC and RNA 
BAM files are specified with the `--atac-bam` and `--rna-bam` options.
These files need to be indexed, for example with [samtools](https://github.com/samtools/samtools).

The alignment records in the BAM files need the barcode tag field to provide
the droplet of origin. These are given by the `--rna-bc-tag` and
`--atac-bc-tag`. By default these are set to `CB` as provided by the 10X
mapping pipeline.

Similarly, the RNA BAM file requires the alignments to have the UMI tags to map 
the alignment to the corresponding molecule. The tag is specified 
with the `--rna-umi-tag` option (set to `UB` by default).

## Mapping filtering parameters

The `--max-nh` filter specifies the maximum number of hits for a query read. Setting
this option to `1` means that only uniquely mapped reads are considered while,
for example, setting this to `10` means that a query with up to 10 alignments is
used for demultiplexing. Setting this to `0` ignores any filtering, and all reads
will be used. The field in the record that the number of hits is found in is
tyipcally `NH` but this can be changed with `--nh-tag`.

An important point is that if reads are multimapped, the primary alignment given by
the aligner is used. All alignments that are tagged as `secondary` or 
`supplementary` in the record's flag [flag](http://www.htslib.org/doc/samtools-flags.html)
are skipped.

The `--min-phredq` parameter specifies the phred-based quality score of the
base calls in the reads. Low-quality base calls from the sequencer can lead to
reading an incorrect base and consequently reading the wrong allele at a
variant. The default minimum is set to 30, meaning a 1 in 1,000 chance of a
sequencing error. This can be set to 20, for example, if a looser but
reasonable criteria is desired.

The `--rna-mapq` and `--atac-mapq` specify the mapping quality score threshold of
the RNA and ATAC BAM files, respectively. This gives
the probability that the alignment position is incorrect. By default, the
threshold is set at 30, meaning a 1 in 1,000 chance of an error in mapping
position.

The `--region` parameter can be used to specify a region of analysis, although
using the entire genome is recommended for normal analyses.

## VCF file

Ambimux takes variant and genotype data from a VCF file. Currently, only the
genotype calls in the `GT` field are used. Ambimux will only use bi-allelic
SNPs and ignore indels. However, VCF files can have multiallelic SNPs present
in two or more records. Ambimux will not ignore these?

To avoid multiallic SNPs, it is recommended to run `bcftools norm`
with the `--multiallelics +snps` option first to merge these multiallelic SNPs.
Then, we can include only biallelic SNPs by filtering with `bcftools view` with
the `-m2 -M2 -v snps` options set. An example filtering strategy is below:

```bash
bcftools view \
    --apply-filters .,PASS \
    -S $sfile \
    -m 2 -M 2 \
    --types snps \
    "$vcf_in" \
| bcftools plugin fill-tags \
    - \
    -- -t AC,AN,AF,MAF \
| bcftools view \
    --include 'INFO/AF > 0 & INFO/AF < 1' \
| bcftools norm \
    --check-ref ws \
    --multiallelics +both \
    --fasta-ref "$fa" \
| bcftools view \
    -m 2 -M 2  \
| bcftools annotate \
    --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' \
    -O z \
> "$vcf_out" \
```
This takes a VCF input file `$vcf_in`, a list of sample IDs to include in
the file `$sfile`, and a fasta file `$fa` and outputs a filtered, normalized
VCF file `$vcf_out`. This filters to include biallelic SNPs only samples
present in `$sfile`, updates the allele frequencies, removes monomorphic SNPs,
normalizes SNPs to a reference genome, merges multiallic SNPs, and again
removes multiallic SNPs. The last annotates SNP IDs by their position and
alleles, and this step can be skipped if original IDs are desired.

Note the chromosome names should match the BAM file.

## Peak and gene annotations

The peak and gene annotations are provided by the `--peaks` and `--gtf`
options. The peaks file should be in BED format while the gene annotations
should be in GTF format. Note that the chromosome names should match with
the BAM file.

These aren't strictly necessary for demultiplexing and will not affect the
singlet assignments. They are used for peak, gene, and allele read counting,
as well as for summary statistics.

## Exclusion regions

ATAC reads that align to problematic regions of the genome can be excluded. 
To do so, provide a BED file of the exclusion regions with the 
`--exclude` option.

A list of exclusion regions used by ENCODE can be found
[here](https://github.com/Boyle-Lab/Blacklist.git).

## Gene, peak, and allele counts

Ambimux generates variant allele counts for the model and these are 
output by default. Additionally, gene and peak counts are produced. Counts are in 
[Matrix Market exchange format](https://math.nist.gov/MatrixMarket/formats.html).
The complete list of counts are 

- SNP variant allele counts per barcode from ATAC-seq
    - `atac.ac.ref.mtx.gz` SNP reference allele counts from ATAC-seq
    - `atac.ac.alt.mtx.gz` SNP alternate allele counts from ATAC-seq
    - `atac.ac.alt.mtx.gz` SNP other allele counts from ATAC-seq
- SNP variant allele counts per barcode from RNA-seq
    - `rna.ac.ref.mtx.gz` SNP reference allele counts from RNA-seq
    - `rna.ac.alt.mtx.gz` SNP alternate allele counts from RNA-seq
    - `rna.ac.alt.mtx.gz` SNP other allele counts from RNA-seq
- Gene counts per barcode from RNA-seq
    - `gc.spl.mtx.gz` Counts from spliced genes from RNA-seq
    - `gc.uns.mtx.gz` Counts from unspliced genes from RNA-seq
    - `gc.amb.mtx.gz` Counts from ambiguous genes from RNA-seq
- Peak counts per barcode from ATAC-seq
    - `atac.peaks.mtx.gz` Counts per peak from ATAC-seq

The rows always correspond to the barcodes, and are listed in 
`barcodes.txt.gz`. The columns for the allele counts correspond to the 
variants and are listed in `var.txt.gz`. Similarly, the columns for the 
gene counts are the genes and are listed in `gene.txt.gz`. ATAC-seq 
peaks are given in the `peaks.txt.gz` file.
The variant 
file lists the first five columns of the VCF file, providing the chromosome, 
position, ID, and reference and alternate alleles.
The gene file includes some information from the GTF file, including position, 
strand, gene type, gene ID, and gene name.

It is also possible to only generate these counts and not run the likelihood 
model by including the `--counts-only` argument.

## EM options

The `--out-min` option specifies the count threshold to to fix droplets as
empty. By default, this option is set to `100`, meaning droplets with less
than 100 RNA UMIs and less than 100 ATAC fragments are fixed as empty. If a
droplet has at least 100 UMIs or 100 fragments, it is considered a
candidate droplet.

The parameters for convergence include `--eps`, which is set to `1e-5` by
default, and `--max-iter`, which is set to `100` by default. The eps parameter
specifies when to stop the EM based on the percent change in the log likelihood.
The max-iter specifies the maximum number of iterations to run, regardless of
convergence by `--eps`.

The `--seed` option can be used to specify the random seed for reproducible
runs.

## Output and filtering

The `summary.txt` output file contains the main demultiplexing results for the
test droplets, removing fixed empty droplets with low coverage. 
This file contains the following columns

1. **Barcode** The barcode of the droplet. For multiome runs, this is usually the RNA barcode.
2. **n_rna_molecules** Total number of RNA UMIs aligned to the genome passing filters.
3. **n_atac_molecules** Total number of ATAC deduplicated fragments aligned to the genome passing filters.
4. **n_features** Number of genes detected in the RNA modality.
5. **n_peaks** Number of peaks detected in the ATAC modality.
6. **n_rna_variants** Number of variants passing filters detected in the RNA modality.
7. **n_atac_variants** Number of variants passing filters detected in the ATAC modality.
8. **atac_pct_mt** Percent of ATAC fragments aligned to the mitochondrial genome.
9. **rna_pct_mt** Percent of RNA UMIs aligned to the mitochondrial genome.
10. **FRIP** Fraction of ATAC fragments inside peaks.
11. **best_type** Droplet type with the maximum likelihood among empty, singlet, and doublet.
  This does not assign ambiguous droplets and further filtering is recommended.
12. **best_sample** Among best droplet type, the sample assignment with the highest likelihood.
13. **best_rna_alpha** The RNA ambient fraction estimate of the best droplet and sample assignment.
14. **best_atac_alpha** The ATAC ambient fraction estimate of the best droplet and sample assignment.
15. **PP0** Posterior probability of an empty droplet type.
16. **PP1** Posterior probability of a singlet droplet type.
17. **PP2** Posterior probability of a doublet droplet type.
18. **LLK0** Log likelihood of an empty droplet type.
19. **LLK1** Log likelihood of a singlet droplet type.
20. **LLK2** Log likelihood of a doublet droplet type.

To call singlets, we recommend applying a filter on the singlet posterior probability.
Thresholds of 0.9, 0.95, or 0.99 are suitable, and ultimately depend on the desired
level of confidence in the call. Alternatively, thresholding the log likelihood
differences can provide another filtering strategy.

In addition to calling singlets, you can filter out highly contaminated droplets
with low coverage. The **alpha** columns contain the ambient fraction estimates.
For example, droplets with **alpha** values greater than 0.25 could be excluded
from analyses.

## Program options

```

ambimux v0.1.0: single-cell demultiplexing

Options:

Input options

  -a, --atac-bam      Indexed ATAC BAM file.
  -r, --rna-bam       Indexed RNA BAM file.
  -v, --vcf           Indexed VCF file.
  -g, --gtf           GTF file.
  -p, --peaks         BED file containing peaks.
  -e, --exclude       BED file containing regions to exclude.
  -s, --samples       File listing sample IDs in VCF file to de-multiplex.

Output options

  -o, --out           Output file prefix [ambimux]. File names are appended with '.' delimiter
  -C, --counts-only   Flag argument to produce counts only and not fit a demultiplexing model.

Alignment options

  -w, --bc-wl         Optional file containing a whitelist of barcodes.
  -u, --rna-umi-tag   RNA BAM tag for UMI [UB].
  -b, --rna-bc-tag    RNA BAM tag for cell barcode [CB].
  -c, --atac-bc-tag   ATAC BAM tag for cell barcode [CB].
  -H, --nh-tag [type] BAM tag for the number of alignments of a read, where type is one 
                      of 'rna' or 'atac'  [NH].

Mapping thresholds

  -m, --max-nh        Only process reads with a maximum of this many alignments [0]. 
                      Set to 0 to ignore (default).
  -P, --min-phredq    Minimum base phred quality score in read to consider for variant calling [30].
  -z, --rna-mapq      Minimum MAPQ (mapping quality) of an RNA alignment [30].
  -Z, --atac-mapq     Minimum MAPQ (mapping quality) of an ATAC alignment [30].
  -R, --region        Region (hts format), for example 'chr21,chr21:10-,chr21-10-20'.

EM options

  -x, --out-min       Calculate the likelihood and demultiplex barcodes that have at least this many 
                      RNA UMIs or ATAC fragments. If there both the UMIs and fragments are below this, skip. [100]
  -h, --eps           Convergence threshold, where the percent change in parameters
                      must be less than this value [1e-5].
  -q, --max-iter      Maximum number of iterations to perform for EM [100].
  -d, --seed          Optional random seed to initialize parameters for EM.
  -T, --threads       Optional number of threads to use [1].

GTF options

  -t, --tx-basic      Read only transcripts tagged as 'basic' in the GTF file.

  -V, --verbose       Write status on output.
      --help          Print this help screen.

```

