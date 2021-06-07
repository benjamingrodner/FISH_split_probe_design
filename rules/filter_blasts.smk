import scripts.blast_functions as fn
import numpy as np

# remove blast result interactions that are not strong
rule filter_blasts:
    input:
        config['output_dir'] + '/{in_file}/blast/{target}/results/{pid}.fasta.blast.out',
        config['output_dir'] + '/{in_file}/target_alignment/{target}.csv'        
    output:
        config['output_dir'] + '/{in_file}/blast/{target}/filtered/{pid}.csv'
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
        # Ignore on-target blasts
        target_alignments = pd.read_csv(input[1])
        on_target_bool = pd.Series(np.zeros(gc_bool.shape))
        for i, r in target_alignments.iterrows():
            a_start, a_end = np.sort([r.sstart, r.send])
            otb_id = (blast.sseqid == r.sseqid)
            otb_start = (blast.sstart >= a_start) & (blast.sstart <= a_end)
            otb_end = (blast.send >= a_start) & (blast.send <= a_end)
            on_target_bool += otb_id & (otb_start | otb_end)
        on_target_bool = np.invert(on_target_bool.astype(bool))
        # Apply boolean filters
        blast_filtered = blast.loc[mch_bool & tm_bool & gc_bool & on_target_bool, :]
        # save the resulting file
        blast_filtered.to_csv(output[0], index=False)