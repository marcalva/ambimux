
# Ambimux

Software tool for demultiplexing single cell multiome data.

## Installation

Platform: ambimux is developed for Linux and tested on CentOS 7.

Clone and build:
```
git clone git@github.com:marcalva/ambimux.git
cd ambimux
make hts   # fetches and builds a pinned htslib release (1.17)
make       # builds ambimux against the local htslib
```
This usually takes no more than a few minutes.

Notes:
- The `make hts` target removes any existing `htslib/` directory and clones a
  pinned release (1.17), then builds static libs.
- ambimux depends on [htslib](https://github.com/samtools/htslib) and therefore on
  htslib’s system prerequisites (autotools, zlib, bzip2, xz, libcurl, etc.). Ensure
  these are installed on your system before running `make hts`.

## Usage

An example usage of ambimux is:
```
ambimux \
    --atac-bam "$atac_bam" \
    --rna-bam "$gex_bam" \
    --vcf "$vcf" \
    --gtf "$gtf" \
    --peaks "$peaks" \
    --exclude "$exclude" \
    --samples "$samples" \
    --bc-wl "$bcsd" \
    --rna-mapq 30 \
    --atac-mapq 30 \
    --tx-basic \
    --eps 1e-6 \
    --max-iter 100 \
    --threads 8 \
    --verbose \
    --out path/to/output/dir/ambimux \
    --out-min 100 \
```
where the variables are replaced by your own files.

Because ambimux retains per‑barcode molecule/fragment data (for deduplication
and counting) and all variant‑overlapping reads used by the ambient model, it
can take longer and use more memory than tools that operate on prefiltered
barcode sets. In one run with seven donors and ~1 billion RNA+ATAC reads,
ambimux completed in roughly 4 hours and used about 50 GB RAM. Actual
performance varies with hardware, --threads, number of barcodes, and variant
density.

You can test ambimux on simulated data using
[ambisim](https://github.com/marcalva/ambisim).

For a detailed explanation of each parameter, see below.

## Input data

### Alignment data

Ambimux can use data from single cell ATAC only, RNA only, or RNA+ATAC
multiome. The ATAC and RNA BAM files are specified with the `--atac-bam` and
`--rna-bam` options. These BAM files need to be coordinate-sorted and indexed,
for example with [samtools](https://github.com/samtools/samtools).

The alignment records in the BAM files need the barcode tag field to provide
the droplet of origin. These are given by the `--rna-bc-tag` and
`--atac-bc-tag`. By default these are set to `CB` as provided by the 10X
mapping pipeline. When running ambimux on multiome data, the barcodes
should match.

Similarly, the RNA BAM file requires the alignments to have the UMI tags to map 
the alignment to the corresponding molecule. The tag is specified 
with the `--rna-umi-tag` option (set to `UB` by default).

You can include a whitelist of barcodes with `--bc-wl` so that counts are
included for all listed barcodes. This is useful when running ambimux across
multiple experiments. Any additional barcodes in the BAM file that are not
found in the given whitelist are still included in the analysis and output.

### Genotype data

Ambimux takes variant and genotype data from a VCF file. Currently, only the
genotype calls in the `GT` field are used. Ambimux will only use bi-allelic
SNPs and ignore indels. However, VCF files can have multiallelic SNPs present
in two or more records. Ambimux will not ignore these and so these should
be removed prior to running ambimux.

One way to remove multiallelic SNPs across multiple records is to run `bcftools norm`
with the `--multiallelics +snps` option to merge these multiallelic SNPs.
We can further include only biallelic SNPs by filtering with `bcftools view` with
the `-m2 -M2 -v snps` options set. An example filtering strategy
using [bcftools](https://samtools.github.io/bcftools/) is:

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
VCF file `$vcf_out`. This filtering removes multiallelic SNPs, includes only samples
present in `$sfile`, updates the allele frequencies, removes monomorphic SNPs,
normalizes SNPs to a reference genome (so that the reference allele matches the
reference sequence in the fasta file), merges multiallelic SNPs, and again
removes multiallelic SNPs. The last step annotates SNP IDs by their position and
alleles, and this step can be skipped if original IDs are desired. Be careful
with strand issues, as these are not fixed here.

The chromosome names must match those in the BAM file.

### Peak and gene annotations

The peak and gene annotations are provided by the `--peaks` and `--gtf`
options. The peaks file should be in BED format while the gene annotations
should be in GTF format. Note that the chromosome names should match exactly with
the BAM file. So a chromosome named `9` will not match a chromosome named `chr9`.

These annotations are used for filtering out reads that align outside of
annotated features (see **Filtering inter-genic/inter-peak reads**).
They are also used for generating peak and gene counts, as well as
for summary statistics.

If the GTF file contains a `basic` tag, you can include only these transcript
isoforms by including the `--tx-basic` command line option. This can be useful
to exclude rare transcripts which aren't representative of the main isoform.

### Excluding problematic regions for ATAC

ATAC reads that align to problematic regions of the genome can be excluded. 
To do so, provide a BED file of the exclusion regions with the 
`--exclude` option. Reads that overlap any listed region will be
discarded.

A list of exclusion regions used by ENCODE can be found
[here](https://github.com/Boyle-Lab/Blacklist.git).

## Mapping filtering parameters

The `--max-nh` filter specifies the maximum number of hits for a query read. Setting
this option to `1` means that only uniquely mapped reads are considered while,
for example, setting this to `10` means that a query with up to 10 alignments is
used for demultiplexing. Setting this to `0` ignores any filtering, and all reads
will be used. The number-of-alignments tag is read separately for RNA and ATAC
using `--rna-nh-tag` and `--atac-nh-tag` (both default to `NH`).

Only the primary alignment for each read is used. For reads mapping to multiple
locations in the genome, only the query alignment with the `primary` flag present
in the [FLAG](http://www.htslib.org/doc/samtools-flags.html)
field is used. One and only one alignment will have this. Reads with
a `secondary` and/or `supplementary` alignment flag are discarded.

The `--min-phredq` parameter specifies the phred-based quality score of the
base calls in the reads. Low-quality base calls from the sequencer can lead to
reading an incorrect base and consequently reading the wrong allele at a
variant. The default minimum is set to 30, meaning a 1 in 1,000 chance of a
sequencing error.

The `--rna-mapq` and `--atac-mapq` ambimux options specify the mapping quality
score threshold of the RNA and ATAC BAM alignments, respectively. This gives the
probability that the alignment position is incorrect. By default, the threshold
is set at 30, which usually means roughly a 1 in 1,000 chance of an error in
mapping position.

The `--region` parameter can be used to specify a region of analysis, although
we strongly recommend using the entire genome for normal analysis.

## Filtering inter-genic/inter-peak reads

Ambimux discards inter-genic and inter-peak reads for demultiplexing by default.
Most single-cell analyses will involve gene and peak read counts, so
these inter-feature reads are implicitly excluded downstream. While most
RNA reads fall within the annotated transcriptome in typical experiments,
60-70% of ATAC reads can fall outside of peaks. Including inter-peak reads
can inflate the ATAC ambient estimates as randomly fragmented exogenous DNA will
more likely fall outside of accessible sites by chance. It is possible
to include all reads using `--mdl-reads all` on the command line. We
also include modeling only inter-feature reads for completeness
using `--mdl-reads inter` instead. The option of including
intra-feature reads `--mdl-reads intra` is set as default.

## Thresholding counts for fixed empty barcodes

The `--out-min` option specifies the count threshold to fix droplets as
empty. By default, this option is set to `100`, meaning droplets with less
than 100 RNA UMIs **and** less than 100 ATAC fragments are fixed as empty. If a
droplet has at least 100 UMIs or 100 fragments, it is considered a
candidate droplet.

## EM parameters

The parameters for convergence include `--eps`, which is set to `1e-6` by
default, and `--max-iter`, which is set to `100` by default. The `eps` parameter
specifies when to stop the EM based on the percent change in the log likelihood.
The `max-iter` option specifies the maximum number of iterations to run, regardless of
convergence by `--eps`.

## Output

The `summary.txt` output file contains the main demultiplexing results for the
test droplets, removing fixed empty droplets with low coverage. 
This file contains the following columns

1. **Barcode** The barcode of the droplet. For multiome runs, this is usually the RNA barcode.
2. **n_rna_molecules** Total number of RNA UMIs aligned to the genome passing filters.
3. **n_atac_molecules** Total number of ATAC deduplicated fragments aligned to the genome passing filters.
4. **n_features** Number of genes detected in the RNA modality.
5. **n_peaks** Number of peaks detected in the ATAC modality.
6. **n_rna_info** Number of informative RNA UMIs overlapping variants.
  This is the number of UMIs that are used in the model for demultiplexing
  and ambient fraction estimation.
7. **n_atac_info** Number of informative ATAC fragments overlapping variants.
  This is the number of fragments that are used in the model for demultiplexing
  and ambient fraction estimation.
8. **n_rna_variants** Number of unique variants detected in the RNA modality.
9. **n_atac_variants** Number of unique variants detected in the ATAC modality.
10. **rna_pct_mt** Percent of RNA UMIs aligned to the mitochondrial genome.
11. **atac_pct_mt** Percent of ATAC fragments aligned to the mitochondrial genome.
12. **FRIG** Fraction of RNA UMIs inside genes.
13. **FRIP** Fraction of ATAC fragments inside peaks.
14. **best_type** Droplet type (`Empty`, `Singlet`, or `Doublet`) with the maximum likelihood.
  This does not assign ambiguous droplets and further filtering is recommended.
15. **best_sample** The sample ID assignment with the highest likelihood of the best type.
  If the `best_type` is `Singlet`, then this will give the singlet assignment with
  the best likelihood. If a doublet, the samples are delimited with a ":"
  character. Ambiguous samples are not given and further filtering is recommended.
16. **best_rna_ambient** The RNA ambient fraction estimated for `best_sample`. This
  will be NA if there are no RNA reads. A fraction of 0 means no contamination.
17. **rna_ambient_info** Effective informative RNA reads for `best_sample` (higher is better; assignment‑specific).
18. **best_atac_ambient** The ATAC ambient fraction estimated for `best_sample`. This
  will be NA if there are no ATAC reads. A fraction of 0 means no contamination.
19. **atac_ambient_info** Effective informative ATAC reads for `best_sample` (higher is better; assignment‑specific).
20. **best_singlet** The best singlet sample with the highest likelihood.
21. **best_doublet** The best doublet samples with the highest likelihood. The
  samples are delimited with a ":" character.
22. **PP0** Posterior probability of being empty.
23. **PP1** Posterior probability of being singlet.
24. **PP2** Posterior probability of being doublet.
25. **LLK0** Log likelihood of being empty.
26. **LLK1** Log likelihood of being singlet.
27. **LLK2** Log likelihood of being doublet.

To call singlets, we recommend applying a filter on the singlet posterior probability.
Thresholds of 0.9, 0.95, or 0.99 are suitable, and ultimately depend on the desired
level of confidence in the call. You can also filter by log likelihood differences.

### Ambient fraction estimation

The contamination estimates per droplet are found in the
**best_rna_ambient** and **best_atac_ambient** columns of the `summary.txt`
output file. These give the fraction (from 0 to 1) of ambient contamination
estimated from the data.
These can be used for downstream analyses as covariates for filtering criteria.

Each ambient estimate is accompanied by an assignment‑specific “info” column
(rna_ambient_info, atac_ambient_info) that reflects the effective number of
informative reads contributing to that estimate. Variants that better
differentiate donors contribute more. Higher values indicate more reliable
ambient estimates.

The ambient estimates in these columns are specific
to the **best_sample** sample and can't be used for other sample assignments.
Each ambient estimate is specific to each assignment. For example,
a singlet assignment to the wrong donor will seem to have higher contamination
than an assignment to the correct donor. The estimates for every assignment
are can be obtained from the `alpha_rna.txt.gz` and `alpha_atac.txt.gz`
files, with the `samples.txt` file giving which columns correspond to which
assignment. An `Empty` droplet will always have a contamination of `1`.

As a quick rule of thumb, `rna_ambient_info` or `atac_ambient_info` below ~10
often yield noisier estimates. The raw informative counts (**n_rna_info**,
**n_atac_info**) remain useful context but are not weighted by variant
discriminability.

### Sample proportions

Ambimux outputs parameter estimates for the sample proportions in the data.
Two parameters are estimated: 1) sample proportions in the singlets/doublets,
and 2) sample proportion in the ambient pool. These parameter estimates
can be found in the `pi.txt.gz` file. The first column contains the sample IDs,
the second column (`Pi`) contains the droplet sample proportions, and the 
third column (`Pi_amb`) contains the ambient sample proportions.

### Empty, singlet, and doublet proportions

The parameter experiment-wide estimates for empty, singlet, and doublet percents
found in the `lambda.txt.gz` file. These give the proportion of empty, singlet,
and doublet droplets in the entire experiment.

### Barcode likelihoods per droplet assignment

The likelihood of each empty, singlet, or doublet possibility is written to
`sample_llk.txt.gz`. Each row corresponds to a test droplet (excluding fixed
empty droplets) and each column corresponds to an assignment. The corresponding
columns can be found in the `samples.txt` file. Each row in `samples.txt`
gives the empty, singlet, and doublet assignments and correspond to each
column in `sample_llk.txt.gz`.

### Gene, peak, and allele counts

Ambimux generates variant allele counts, gene counts, and peak counts. These
can be excluded from the output by adding the `--no-counts-o` command line
option. Counts are in [Matrix Market exchange
format](https://math.nist.gov/MatrixMarket/formats.html). The complete list of
counts are 

- SNP variant allele counts per barcode from ATAC-seq
    - `atac.ac.ref.mtx.gz` SNP reference allele counts from ATAC-seq
    - `atac.ac.alt.mtx.gz` SNP alternate allele counts from ATAC-seq
    - `atac.ac.oth.mtx.gz` SNP other allele counts from ATAC-seq
- SNP variant allele counts per barcode from RNA-seq
    - `rna.ac.ref.mtx.gz` SNP reference allele counts from RNA-seq
    - `rna.ac.alt.mtx.gz` SNP alternate allele counts from RNA-seq
    - `rna.ac.oth.mtx.gz` SNP other allele counts from RNA-seq
- Gene counts per barcode from RNA-seq
    - `gc.spl.mtx.gz` Counts from spliced genes from RNA-seq
    - `gc.uns.mtx.gz` Counts from unspliced genes from RNA-seq
    - `gc.amb.mtx.gz` Counts from ambiguous genes from RNA-seq
- Peak counts per barcode from ATAC-seq
    - `atac.peaks.mtx.gz` Counts per peak from ATAC-seq

The columns always correspond to the barcodes, and are listed in 
`barcodes.txt.gz`. The rows for the allele counts correspond to the 
variants and are listed in `var.txt.gz`.
This variant 
file lists the first five columns of the VCF file, providing the chromosome, 
position, ID, and reference and alternate alleles.
Similarly, the rows (genes) for the 
gene counts are listed in `gene.txt.gz`.
The gene file includes some information from the GTF file, including position, 
strand, gene type, gene ID, and gene name.
ATAC-seq peaks are given in the `peaks.txt.gz` file, including the chromosome,
start, and end positions.

It is also possible to only generate these counts and not run demultiplexing
by including the `--counts-only` argument.

## Program options

```

ambimux v0.5.0: single-cell demultiplexing

Options:

Input options

  -a, --atac-bam      Indexed ATAC BAM file.
  -r, --rna-bam       Indexed RNA BAM file.
  -v, --vcf           Indexed VCF file.
  -g, --gtf           GTF file.
  -p, --peaks         BED file containing peaks.
  -e, --exclude       BED file containing ATAC regions to exclude.
  -s, --samples       File listing sample IDs in VCF file to de-multiplex.

Output options

  -o, --out           Output file prefix. File names are appended with '.' delimiter [ambimux].
  -C, --counts-only   Flag argument to produce counts only and not fit a demultiplexing model.
  -n, --no-counts-o   Do not output variant allele, gene, and peak counts. by default, ambimux
                      output UMI/fragment counts.

Alignment options

  -w, --bc-wl         Optional file containing a whitelist of barcodes. One barcode per line.
  -u, --rna-umi-tag   RNA BAM tag for UMI [UB].
  -b, --rna-bc-tag    RNA BAM tag for cell barcode [CB].
  -c, --atac-bc-tag   ATAC BAM tag for cell barcode [CB].
  -N, --rna-nh-tag   RNA BAM tag for number of alignments [NH].
  -A, --atac-nh-tag  ATAC BAM tag for number of alignments [NH].

Mapping thresholds

  -m, --max-nh        Only process reads with a maximum of this many alignments. 
                      Set to 0 to ignore (default) [0].
  -P, --min-phredq    Minimum base phred quality score in read to consider for variant calling [30].
  -z, --rna-mapq      Minimum MAPQ (mapping quality) of an RNA alignment [30].
  -Z, --atac-mapq     Minimum MAPQ (mapping quality) of an ATAC alignment [30].
  -R, --region        Region (hts format), for example 'chr21,chr21:10-,chr21-10-20'.

Model options

  -x, --out-min       Calculate the likelihood and demultiplex barcodes that have at least this many 
                      RNA UMIs or ATAC fragments. If both the UMIs and fragments are below this, skip [100].
  -h, --eps           Convergence threshold, where the percent change in parameters
                      must be less than this value [1e-6].
  -j, --amb-eps       Convergence threshold for droplet ambient contamination parameter alpha,
                      the percent change in value must be less than this argument [1e-6].
  -q, --max-iter      Maximum number of iterations to perform for EM [100].
  -d, --amb-prior-w   Weight of the ambient fraction prior. Effectively gives the number
                      of reads added to the likelihood as a prior. Must be > 0 [1e-8].
  -i, --mdl-reads     The type of reads to use for demultiplexing and ambient estimation.
                      One of 'intra', 'inter', or 'all'. 'intra' (default) specifies the model
                      should use only reads contained in peaks or genes.
                      'inter' specifies to use only inter-peak or inter-gene reads that
                      fall outside peaks or genes. 'all' specifies to use all reads including
                      inter- and intra- gene and peak reads [intra].
  -T, --threads       Optional number of threads to use [1].

GTF options

  -t, --tx-basic      Read only transcripts tagged as 'basic' in the GTF file.

  -V, --verbose       Write status on output.
      --help          Print this help screen.

```
