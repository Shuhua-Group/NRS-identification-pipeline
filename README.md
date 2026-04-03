# NRS Identification Pipeline User Guide

## Pipeline Overview

The NRS (Non-Reference Sequence) identification pipeline is a two-stage analysis pipeline for identifying and analyzing novel sequence regions in the human genome. NRS represents sequence information missing from current linear reference genomes (GRCh37/GRCh38) and is a crucial component of the pangenome.

The human pangenome aims to capture the complete genetic variation within human populations, and NRS, as a key element, is of significant importance for understanding human genetic diversity, disease susceptibility, and drug responses. This pipeline systematically identifies and analyzes these novel sequence regions through advanced analytical methods, providing a powerful tool for pangenome research.

The pipeline consists of two main parts:

- **Part1**: Non-Reference Sequence of human pangenome analysis, including decontamination, contig region analysis, and end mapping
- **Part2**: Clustering analysis, merging results from different samples and performing clustering

## Installation Requirements

### Part1 Dependencies
- g++ 10.3.0 or higher
- samtools
- blast
- megahit or popins2
- nucmer
- bowtie2

### Part2 Dependencies
- Python 3.6+
- Biopython
- threading module

## Running Steps

### Input File Description
- **Part1 Input Files**: BAM files after duplication marking in the sample NGS pipeline
- **Part2 Input Files**: Part1 output result folder

### Part1: Human Pangenome Analysis

#### Method 1: Using the Executable (Recommended)
```bash
# Place the executable from the bin folder into your environment
pangenomeDebug5.Part.clean --indir /path/to/input/dir --outdir /path/to/outdir
```

#### Method 2: Using Sample Path File
```bash
pangenomeDebug5.Part.clean --inPath /path/to/inputPathFile --outdir /path/to/outdir
```

#### Method 3: Compile Yourself
1. Ensure g++ version ≥ 10.3.0
2. Enter the src directory, modify the GXX parameter in Makefile
3. Execute `make` to compile
4. Use the compiled executable

### Part2: Clustering Analysis

```bash
python pangenomePart2.py --indir /path/to/part1/output --outdir /path/to/part2/output
```

## Result Description

### Part1 Results
- `outdir/pan_<sampleID>/contigRegion/unmappedContig`
- `outdir/pan_<sampleID>/contigRegion/BothEndMapping`
- `outdir/pan_<sampleID>/contigRegion/leftEndmapping`
- `outdir/pan_<sampleID>/contigRegion/rightEndmapping`
- Decontaminated contig: `outdir/pan_<sampleID>/contigRegion/novelContig_<sampleID>.fa`

### Part2 Results
- `outdir/leftEnding/UpRepCluTot2.fa`
- `outdir/rightEnding/UpRepCluTot2.fa`
- `outdir/bothEnding/UpRepCluTot.fa`
- `outdir/unmapped/UpG38Unmap.fa`
- `MergeLefRig/MergeContain.position`
- `MergeLefRig/MergeOverlap.position`

## Notes

1. Do not include "-" in input file names, as Part2 uses names for identification
2. To modify parameters, copy the ConfigFile to your own directory and use `--config /path/to/your/config` parameter
3. The default decontamination method is fast, which only searches human or ape datasets
4. Part2 uses 10 threads for G38 filtering by default
5. MergeLefRig/MergeOverlap.position file may have no results due to popins call failure

## Configuration File Description

Part1 configuration file contains the following main parameters:
- `IniEnv`: Program running environment
- `buildThread`: Number of threads
- `suffix`: Input file suffix
- `vrmt`: Decontamination operation version
- `blastIdentity`: Blast alignment quality parameter
- `blastLength`: Blast alignment length parameter
- `Assembletool`: Assembly tool choice (megahit or popins)
- `PopinsPath`: popins2 path

## Example Commands

### Part1 Example
```bash
# Method 1
pangenomeDebug5.Part.clean --indir /path/to/sample/bam/files --outdir /path/to/output/directory

# Method 2
pangenomeDebug5.Part.clean --inPath /path/to/sample/path/file --outdir /path/to/output/directory
```

### Part2 Example
```bash
python pangenomePart2.py --indir /path/to/part1/output --outdir /path/to/part2/output
```

## Version Information

- Part1: human_pangenomeDebug5
- Part2: Part2.2.Mt.reWriteLeftRight
