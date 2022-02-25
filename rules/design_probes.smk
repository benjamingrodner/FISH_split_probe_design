from scripts.primer3_adapter import Primer3

# Design probes for each target sequence using primer3
rule design_probes:
    input:
        config['probe_generate_dir'] + '/{in_file}/target_files/{target}.fasta'
    output:
        config['probe_generate_dir'] + '/{in_file}/primer3/{target}/{target}.fasta',
        config['probe_generate_dir'] + '/{in_file}/primer3/{target}/{target}.csv'
    run:
        min_size, opt_size, max_size = [config['probe_length'][s]
                                        for s in ['min_size', 'opt_size', 'max_size']]
        p3_object = Primer3(fasta_filename=input[0],
                            min_size=min_size, opt_size=opt_size, max_size=max_size,
                            )
        target_dir = os.path.split(output[0])[0]
        p3_dir = os.path.split(target_dir)[0]
        print('\n LOOK HERE! \n',p3_dir)
        p3_object.design_probes(output_dir=p3_dir)
