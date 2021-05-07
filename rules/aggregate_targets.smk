# Ensures the checkpoint gets run for prep_target_files
rule aggregate_targets:
    input:
        aggregate_targets
    output:
        config['output_dir'] + '/primer3/{in_file}/aggregate.txt'
    shell:
        "echo {input} > {output}"