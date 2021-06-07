rule intermediate:
    input:
        config['output_dir'] + '/{in_file}/spacer_selection/{target}/blast/inputs/{spacer}.fasta'
    output:
        config['output_dir'] + '/{in_file}/spacer_selection/{target}/blast/outputs/{spacer}.fasta'
    shell:
        "cp {input} {output}"