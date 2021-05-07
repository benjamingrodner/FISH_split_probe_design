import scripts.blast_functions as fn

# remove blast result interactions that are not strong
rule filter_blasts:
    input:
        config['output_dir'] + '{in_file}/blast/{target}/results/{pid}.fasta.blast.out',
        config['input_dir'] + '{in_file}/on_target_ids.csv',

    output:
        config['output_dir'] + '{in_file}/blast/{target}/filtered/{pid}.csv'
    run:
        # Read in blast table
        blast = fn.read_blast_table(input[0])
        # calculate mch, tm and gc contents
        blast['MCH'] = blast.apply(lambda x: fn.max_continuous_homology(x.qseq, x.sseq), axis=1)
        blast['TM'] = blast.apply(lambda x: fn.melting_temperature(x.qseq, x.sseq), axis=1)
        blast['GC'] = blast.apply(lambda x: fn.gc_content(x.qseq, x.length), axis=1)
        # remove all except the strongest interactinos
        mch_bool = blast.MCH > config['mch_filter']
        tm_bool = blast.TM > config['tm_filter']
        gc_bool = blast.GC > config['gc_filter']
        blast_filtered = blast.loc[mch_bool & tm_bool & gc_bool, :]
#         # Remove on target blasts
#         on_target_ids = pd.read_csv(input[1]).values.astype(list)
#         sseqids = blast.loc[~blast.sseqid.isin(on_target_ids)]
        # save the resulting file
        blast_filtered.to_csv(output[0], index=False)