# visual_plodinf
Script used to verify ploidy visually from bam data

Quick implementation of an approach i tried for estimating the plodiy of Illumina sequencing samples. I tested with ddRAD-seq data but should work with any "random/pseudorandom" sequencing data. The idea is to plot the observed allele proportion of all heterozygous sites found in a bam file. BCFtools is used to do a "naive" variant calling by finding sites where more than 1 base is observed, and extracting the counts of those bases. The allele proprotion for a particular base at a site would be the count of reads of said base divided by the total site depth.

The shape of the distribution should give you an idea of the ploidy of the sample as depending on the ploidy you expect well defined. The smallest allele proportion you can theoretically observe at heterozygous sites for a particular ploidy "p" is 1/p, and the highest is (p-1)/p. So the distribution of allele proportions should have peaks near x/p with x being a natural number in in the range [1,p-1]. So for a diplod a single peak at 0.5 should be observed, a triploid would have 2 peaks at 0.33 and 0.66 and a tetraploid would have three peaks at 0.25, 0.5, and 0.75. 

In practice, sampling variablity, genome heterozygosity and sequencing error will make the signal noisy, an maybe not even possible to discern. 

# Dependencies

"environment.yml" specifies the needed dependencies for a conda environment:

  - matplotlib>=3.10.3
  - bcftools>=1.22

# Usage

```
bam_file="/path/to/sorted/bam_file.bam"
out_path="/out/folder/path"
minimum_ad="15"
minimum_depth="30"
maximum_depth="80"
minimum_baseQ="20"
minimum_mapQ="30"
histogram_bins="41"

bash main_script.sh $bam_file $out_path $minimum_ad $minimum_depth $maximum_depth $minimum_baseQ $minimum_mapQ $histogram_bins
```

The program will generate a tsv file with allele proportion data per heterozygous site. And then "histogram_generator.py" will plot the histogram from the data.  

## Parameter explanation

- bam_file: Path to a sorted bam file
- out_path: Path to output folder
- minimum_ad: Minimum number of required reads per allele to include a site into the analysis
- minimum_depth: Minimum depth of a site to include into the analysis
- maximum_depth: Maximum depth of a site to include into the analysis
- minimum_baseQ: Minimum quality of a read base to include into the count of alleles
- minimum_mapQ: Minimum mapping quality(as defined in the bam file, varies by mapper) to include into the count of alleles
- histogram_bins: Number of bins to create the histogram 31 to 41 is a reasonable range.


# Caveats

- This approach does a naive variant calling so tune the following parameters to reasonable levels to minimize possible errors:
    - minimum_ad
    - minimum_depth
    - maximum_depth
    - minimum_baseQ
    - minimum_mapQ
- The script should be run from the repository folder so that the python scripts are found by the bash script.
- The whole script runs in a single thread and depending on the size of the bam it may be slow, if you have multiple samples try running them in parallel. 

