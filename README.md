# README for Split Probe Design Pipeline
## Purpose
This pipeline is used to design "split probes" as per [HCR v3](http://doi.org/10.1242/dev.165753) and [split-FISH](https://doi.org/10.1038/s41592-020-0858-0).
The advantage here is that both probes are required to bind to the target site in order for fluorescence signal to appear.
This is meant to increase the probe specificity.
This probe design is based on brute-force design of all possible probes, then filtering based on specificity.

## Layout

The pipeline runs by executing the [Snakefile](https://github.com/benjamingrodner/FISH_split_probe_design/blob/main/Snakefile). The Snakefile refers to a series of [rules](https://github.com/benjamingrodner/FISH_split_probe_design/tree/main/rules), which directly run chunks of python code as well as specify the input and output files for the chunks. The code in the rules executes python functions loaded as modules from [scripts](https://github.com/benjamingrodner/FISH_split_probe_design/tree/main/scripts). The Snakefile, rules, and functions reference the [configuration file](https://github.com/benjamingrodner/FISH_split_probe_design/blob/main/config_example.yaml) to determine the specific filepaths and variable settings for the run. 

## Package requirements
[primer3](https://github.com/primer3-org/primer3) - Make sure to add the path to your primer3 executable to the configuration file

[blastn](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

snakemake=6.1.0

python=3.6

numpy=1.15.4

pandas=0.24.1

biopython=1.72

openpyxl=3.0.2


## Configuration
Use the config_example.yaml as a template and put your configuration file in the directory where you want to design probes. Make all of your filepaths in the configuration file relative to the directory where you want to design probes. Make sure to add the path to your primer3 executable to the configuration file. 

Create a blast database using
```
makeblastdb -in {DATABSE_FILE} -dbtype nucl -parse_seqids
```
Your inputs must be fasta files and there can be multiple targets in each input fasta.


## Execution
Execute the following in the command line


```
cd {PROBE_DESIGN_DIRECTORY}
snakemake --snakefile {PATH_TO_PIPELINE_DIRECTORY}/Snakefile \
          -j {NUMBER_OF_CORES} \
          --configfile {CONFIG_FILENAME}
```

### Run probe design on example data 

    snakemake --configfile config_example.yaml -j 4


## Workflow
Below is a rule graph of the Snakemake pipeline. Note that the rule 'prep_probe_files' is a checkpoint with the rules 'blast_probes' and 'filter_blasts' dependent on its output and thus not evaluated in the DAG until 'prep_probe_files' is complete. Thus 'blast_probes' and 'filter_blasts' fit between 'prep_probe_files' and 'evaluate_pairs'. The rule 'filter_blasts' also depends on the output of the rule 'target_alignments', which is also not shown in the rule graph.


[](./rulegraph.svg)
<img src="./rulegraph.svg">

## Output File Structure 

          output_dir
          ├── crosstalk
          │   └── TARGET_FASTA_BASENAME
          │       ├── evaluated_crosstalk
          │       └── filtered_crosstalk
          ├── final_outputs
          │   └── TARGET_FASTA_BASENAME
          │       └── selection_table
          ├── helper_selection
          │   └── TARGET_FASTA_BASENAME
          │       ├── primer3
          │       │   ├── reverse_complement
          │       │   │   └── FASTA_ENTRY_ID
          │       │   └── FASTA_ENTRY_ID
          │       └── selection
          ├── probe_generate
          │   └── TARGET_FASTA_BASENAME
          │       ├── all_pairs_index
          │       ├── blast
          │       │   └── FASTA_ENTRY_ID
          │       │       ├── probes
          │       │       └── results
          │       ├── primer3
          │       │   └── FASTA_ENTRY_ID
          │       └── target_files
          │           └── reverse_complement
          ├── probe_selection
          │   └── TARGET_FASTA_BASENAME
          │       ├── blast
          │       │   └── FASTA_ENTRY_ID
          │       │       ├── filtered
          │       │       └── measured
          │       ├── evaluated_pairs
          │       └── selected_pairs
          ├── spacer_selection
          │   └── TARGET_FASTA_BASENAME
          │       ├── blast
          │       │   └── FASTA_ENTRY_ID
          │       │       ├── filtered
          │       │       ├── inputs
          │       │       └── outputs
          │       └── spacer_selection
          └── target_alignment
              └── TARGET_FASTA_BASENAME
                  ├── blast
                  └── filtered
