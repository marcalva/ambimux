
# Ambimux

Software tools to process aligned multiome data.

## Installation

The source code can be downloaded using
```
git clone --recurse-submodules git@github.com:marcalva/ambimux.git
```
The `--recurse-submodules` flag makes sure the `htslib` dependency is downloaded.

To create the binary, run `make`. After a successful build, the `ambimux` 
binary file will be present in the current directory.

The only dependency is htslib, which is 
included as a subdirectory. Note that ambimux requires htslib to be built 
successfully during `make`.

## Usage

An example usage of ambimux is
```
./ambimux \
    --rna-bam $rna_bam \
    --atac-bam $atac_bam \
    --vcf $vcf \
    --gtf $gtf \
    --peaks $peaks \
    --exclude $excl \
    --samples $samples \
    --rna-mapq 30 \
    --atac-mapq 30 \
    --out-min 100 \
    --bc-wl $wl_bcs \
    --threads 8 \
    --verbose \
    --eps 10 \
    --max-iter 10 \
    --out mulitome
```
where the variables are replaced by your own files.

## Alignment data

Ambimux can use data from single cell ATAC, RNA, or both. The ATAC and RNA 
BAM files are specified with the `--atac-bam` and `--rna-bam` options.
These files need to be indexed, for example with [samtools](https://github.com/samtools/samtools).

The alignments in the BAM files need the barcode tags to map to its 
corresponding droplet. These are given by the `--rna-bc-tag` and 
`--atac-bc-tag`. By default these are set to `CB` as provided by the 10X 
mapping pipeline.

Similarly, the RNA BAM file requires the alignments to have the UMI tags to map 
the alignment to the corresponding molecule. The tag is specified 
with the `--rna-umi-tag` option (set to `UB` by default).

## Mapping thresholds

There are two main thresholds to consider with ambimux. The first is the
phred-based quality score of the base calls in the reads. Low-quality base
calls from the sequencer can lead to reading an incorrect base and consequently
reading the wrong allele at a variant. The default minimum is set to 30,
meaning a 1 in 1,000 chance of a sequencing error.

The second threshold is the mapping quality score of the alignment.  This gives
the probability that the alignment position is incorrect.  By default, the
threshold is set at 30, meaning a 1 in 1,000 chance of an error in mapping
position.

## VCF file

Ambimux takes variant and genotype data from a VCF file. It uses either hard
genotype calls in the `GT` field or genotype probabilities in the `GP` field.
Ambimux currently only uses bi-allelic SNPs and ignores indels. However, VCF
files can have multiallelic SNPs split into two or more records. Ambimux will
not ignore these SNPs and will incorrectly treat them as independent, with a
warning.  To avoid multiallic SNPs, it is recommended to run `bcftools norm`
with the `--multiallelics +snps` option first to merge these multiallelic SNPs.
Then, we can include only biallelic SNPs by filtering with `bcftools view` with
the `-m2 -M2 -v snps` options set.

## Peak and gene annotations

Peak annotations are input as a BED file with the `--peaks` option.
Gene annotations are input as a GTF file with the `--gtf` option.

## Exclusion regions

Reads that align to problematic regions of the genome can be excluded. 
To do so, provide a BED file of the exclusion regions with the 
`--exclude` option.

## Counting reads

Ambimux generates gene and variant allele counts for the model and these are 
output by default. Counts are in 
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

## Output filtering

Droplets with very few counts are not useful for analyses and are typically 
empty. We avoid calculating the likelihood for these droplets to save time 
and resources. To specify the threshold at which to keep or remove droplets, 
add the `--out-min` option, which is set at `100` by default.
If a barcode has at least this many RNA UMIs or this many ATAC fragments, 
it is included in the likelihood calculation. If neither UMIs or fragments 
pass this threshold, the barcode is skipped.

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
  -x, --out-min       Calculate the likelihood and demultiplex barcodes that have at least this many 
                      RNA UMIs or ATAC fragments. If there both the UMIs and fragments are below this, 
skip. [100]

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

  -h, --eps           Convergence threshold, where the percent change in parameters
                      must be less than this value [1e-5].
  -q, --max-iter      Maximum number of iterations to perform for EM [20].
  -d, --seed          Optional random seed to initialize parameters for EM.
  -T, --threads       Optional number of threads to use [1].

GTF options

  -t, --tx-basic      Read only transcripts tagged as 'basic' in the GTF file.

  -V, --verbose       Write status on output.
      --help          Print this help screen.

```

