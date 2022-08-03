from scripts.primer3_adapter import Primer3


########################################################################
# Design all possible helper probes with primer3
########################################################################

rule design_helpers:
    input:
        config['probe_generate_dir'] + '/{in_file}/target_files/{target}.fasta'
    output:
        config['helper_select_dir'] + '/{in_file}/primer3/{target}/{target}.fasta',
        config['helper_select_dir'] + '/{in_file}/primer3/{target}/{target}.csv',
        config['probe_generate_dir'] + '/{in_file}/target_files/reverse_complement_helper/{target}.fasta',
        config['helper_select_dir'] + '/{in_file}/primer3/reverse_complement/{target}/{target}.fasta',
        config['helper_select_dir'] + '/{in_file}/primer3/reverse_complement/{target}/{target}.csv'
    run:
        # set a large range of helper sizes
        min_size, opt_size, max_size = 15, 20, 30
        # Check if we want to design helper probes
        if config['get_helpers']['same_strand']:
            # and rund the scripts
            p3_object = Primer3(fasta_filename=input[0],
                                min_size=min_size, opt_size=opt_size, max_size=max_size,
                                )
            target_dir = os.path.split(output[0])[0]
            p3_dir = os.path.split(target_dir)[0]
            p3_object.design_probes(output_dir=p3_dir)
        else:
            # Otherwise write empty files
            open(output[0], 'a').close()
            open(output[1], 'a').close()

        # check if we want to design helpers for the other strand
        basename, ext = os.path.splitext(input[0])
        if config['get_helpers']['other_strand']:
            # if so, make a reverse complement fasta
            if ~os.path.exists(output[2]):
                records = SeqIO.parse(input[0], ext[1:])
                for r in records:
                    target = r.id
                    seq = r.seq.reverse_complement()
                    new_record = SeqRecord(seq, id=target)
                    SeqIO.write(new_record, output[2], 'fasta')
                    print('Wrote: ', output[2])
            # Then run the scripts on the new file
            p3_object = Primer3(fasta_filename=output[2],
                                min_size=min_size, opt_size=opt_size, max_size=max_size,
                                )
            target_dir = os.path.split(output[3])[0]
            p3_dir = os.path.split(target_dir)[0]
            p3_object.design_probes(output_dir=p3_dir)
        else:
            # otherwise write empty files
            open(output[2], 'a').close()
            open(output[3], 'a').close()
            open(output[4], 'a').close()
