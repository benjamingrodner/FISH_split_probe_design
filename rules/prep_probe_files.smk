# Split probes into individual fasta files in preparation for blasting
checkpoint prep_probe_files:
    input:
        config['output_dir'] + '/{in_file}/primer3/{target}/{target}.fasta'
    output:
        directory(config['output_dir'] + '/{in_file}/blast/{target}/probes')
    run:
        records = SeqIO.parse(input[0], 'fasta')
        for r in records:
            probe_id = r.id
            if not os.path.exists(output[0]):
                os.makedirs(output[0])
            out_fn = output[0] + '/' + probe_id + '.fasta'
            SeqIO.write(r, out_fn, 'fasta')