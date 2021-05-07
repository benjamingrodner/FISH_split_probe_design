# Ensures blast_probes gets run
rule aggregate_probes:
    input:
        aggregate_filtered_probe_blasts
    output:
        config['output_dir'] + '/{in_file}/blast/aggregates/{target}.txt'
    shell:
        "cat {input} > {output}"