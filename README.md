# Rieux Short Reads Mapping 🧬🚀

[![DOI](https://zenodo.org/badge/1194379824.svg)](https://doi.org/10.5281/zenodo.19281220)

A robust, one-click bioinformatics pipeline for highly specific short-read mapping and target quantification in High-Performance Computing (HPC) environments running SLURM. 

Designed specifically for low-biomass clinical samples or complex microbiomes, this tool features **Dynamic Resource Allocation**, **Strict Competitive Mapping**, and a **Center of Mass (Midpoint)** strategy for absolute quantitative precision without double-counting reads at genomic window boundaries.

## 🌟 Key Features

* **Auto-Scaling (AI Resource Allocation):** Automatically calculates required RAM, Threads, and SLURM Partitions based on the exact size of your input FASTQ files and Reference Genomes. Includes built-in safety limits to prevent cluster bans.
* **Array-Based Horizontal Scaling:** Processes dozens or hundreds of samples simultaneously in isolated scratch environments, merging everything automatically upon completion.
* **Strict Competitive Mapping:** Maps reads simultaneously against a target and a sister-species (e.g., *Treponema* vs. *Leptospira*). Crucially, it employs a **Strict Biological Filter** to extract and save *only* the true target reads into the final BAM, completely discarding reads that tied or mapped to the negative control (like conserved rRNAs).
* **Center of Mass Quantification:** Converts reads to 1bp midpoints before intersecting with genomic features, eliminating the "double-counting" artifact in sliding window profiling.
* **Smart Aggregation:** Automatically merges contiguous 500bp windows belonging to the same genomic feature into single blocks for accurate gene-level analysis.

## 🛠️ Dependencies

The pipeline assumes the following modules are available on your HPC:
* `fastqc` (Tested on 0.12.1)
* `fastp` (Tested on 0.23.4)
* `bwa-mem2` (Tested on 2.3)
* `samtools` (Tested on 1.22+)
* `bedtools` (Tested on 2.31+)
* `multiqc`

## 📦 Installation & Permanent Activation

To install the tool and make it permanently available as a native command in your HPC terminal session, follow these steps:

1. **Clone the repository** to your preferred directory (e.g., `~/tools`):
   ```bash
   cd ~
   git clone [https://github.com/thiagoparentefiocruz/rieux_short_reads_mapping.git](https://github.com/thiagoparentefiocruz/rieux_short_reads_mapping.git)
   ```

2. Add the tool to your .bashrc so it loads automatically on every login:
   ```bash
   echo "source ~/rieux_short_reads_mapping/rieux_short_reads_mapping.sh" >> ~/.bashrc
   ```

3. Reload your bash profile to activate it immediately in your current session:
   ```bash
   source ~/.bashrc
   ```

You can now use rieux_short_reads_mapping from any folder in your cluster!

🚀 Usage & Parameters
The pipeline is invoked with a single command.

Mandatory Arguments
-i : Directory containing raw paired-end .fastq.gz files (searches recursively).

-r : Target Reference Genome (.fasta or .fna).

-g : Target Genome Annotation (.gff3). Essential for feature profiling.

Optional Arguments
-h : Host Reference Genome (.fasta). Reads mapping to this host will be depleted before target analysis.

-c : Control Reference Genome (.fasta). Concatenated with the target to force competitive mapping.

-m : Minimum Mapping Quality (MAPQ). Excludes ambiguously mapped reads. (Default: 0. Recommended with -c: 20).

-o : Output Directory. (Default: ./resultados_pipeline).

Quick Start Example
```bash
rieux_short_reads_mapping \
  -i ./raw_data \
  -r ./ref/target_pathogen.fasta \
  -g ./ref/target_pathogen.gff3 \
  -h ./ref/human_host.fasta \
  -c ./ref/sister_species_control.fasta \
  -m 20 \
  -o ./final_results
  ```

📊 Outputs (The final_reports/ Folder)
To keep your workspace clean, the pipeline extracts the 4 most important, biologically relevant files and places them in a final_reports/ directory inside your output folder:

* rastreabilidade_reads.tsv: The ultimate tracking table showing survival rates from raw reads to pure target mapping, including Effective Depth (X) and Mean MAPQ.
   
* resumo_super_features.tsv: The high-level infection profile. Summarizes absolute read counts per broad genomic feature (e.g., CDS, rRNA, tRNA, intergenic) per sample.

* matriz_blocos_contiguos.tsv: The aggregated feature matrix. Contiguous 500bp windows sharing the same annotation are merged into single blocks to prevent splitting genes across rows.

* multiqc_report.html: Interactive HTML quality report for pre- and post-trimming read metrics.

(Note: Intermediate scratch files and raw unaggregated matrices are automatically cleaned up to save disk space, but pure target .bam files are kept in the 05_mapping_target folder for manual inspection in IGV if needed).

🏗️ Architecture
The tool dynamically creates and submits a 3-stage SLURM dependency chain:
* Job 0 (Index): Safely generates competitive indexes on a compute node.

* Job 1 (Array): Spawns isolated tasks for each sample. Throttled to a maximum of 15 concurrent jobs to respect typical HPC group memory limits and avoid AssocGrpMemLimit locks.

* Job 2 (Merge): Automatically triggers upon array success, harvesting outputs, running bedtools intersections, aggregating contiguous blocks, and extracting the golden final_reports/.

Developed by thiagoparentefiocruz
