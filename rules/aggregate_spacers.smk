# Ensures blast_probes gets run
rule aggregate_spacers:
    input:
        aggregate_filtered_spacer_blasts
    output:
        config['output_dir'] + '/{in_file}/spacer_selection/{target}/aggregate.txt'
    shell:
        "cat {input} > {output}"