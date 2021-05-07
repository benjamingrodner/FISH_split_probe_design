# README for Split Probe Design Pipeline
## Configuration
Use the config_example.yaml as a template and put your configuration file in the directory where you want to design probes. Make all of your filepaths in the configuration file relative to the directory where you want to design probes. 

Create a blast database using

    makeblastdb -in {DATABSE_FILE} -dbtype nucl -parse_seqids

Your inputs must be fasta files and there can be multiple targets in each input fasta.

## Package requirements
snakemake=6.1.0

python=3.6

numpy=1.15.4

pandas=0.24.1

biopython=1.72

openpyxl=3.0.2

# Execution
Execute the following in the command line

     cd {PROBE_DESIGN_DIR}
     snakemake --snakefile {PATH_TO_PIPELINE_DIR}/Snakefile \
               -j {NUMBER_OF_CORES} \
               --configfile {CONFIG_FILENAME}
