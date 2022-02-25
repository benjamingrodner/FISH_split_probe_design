from scripts.blast_functions import blastProbes

# Blast each probe
rule blast_probes:
    input:
        config['probe_generate_dir'] + '{in_file}/blast/{target}/probes/{pid}.fasta'
    output:
        config['probe_generate_dir'] + '{in_file}/blast/{target}/results/{pid}.fasta.blast.out'
    run:
        output_dir = os.path.split(output[0])[0]
        blastProbes(query=input[0], database=config['blast_database'],
                    output_dir=output_dir, strand=config['blast_strand'])
