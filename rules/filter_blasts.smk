import scripts.blast_functions as fn
import numpy as np

# remove blast result interactions that are not strong
rule filter_blasts:
    input:
        config['probe_generate_dir'] + '/{in_file}/blast/{target}/results/{pid}.fasta.blast.out',
        config['target_alignment_dir'] + '/{in_file}/filtered/{target}.csv'
    output:
        config['probe_select_dir'] + '/{in_file}/blast/{target}/measured/{pid}.csv',
        config['probe_select_dir'] + '/{in_file}/blast/{target}/filtered/{pid}.csv'
    run:
        # Read in blast table
        blast = fn.read_blast_table(input[0])
        # calculate mch, tm and gc contents
        if blast.shape[0] > 0:
            target_alignments = pd.read_csv(input[1])
            blast_measured = fn.measure_blasts(blast, target_alignments)
            blast_filtered = fn.filter_blasts(blast_measured,
                                              mch_filter=config['mch_filter'],
                                              gc_filter=config['gc_filter'],
                                              tm_filter=config['tm_filter'])
            # blast['MCH'] = blast.apply(lambda x: fn.max_continuous_homology(x.qseq, x.sseq), axis=1)
            # blast['TM'] = blast.apply(lambda x: fn.melting_temperature(x.qseq, x.sseq), axis=1)
            # blast['GC'] = blast.apply(lambda x: fn.gc_content(x.qseq, x.length), axis=1)
            # # Sort blast start and end
            # blast_startend = blast.loc[:,['sstart','send']].values
            # blast_startend.sort(axis=1)
            # blast['sstart_sort'] = blast_startend[:,0]
            # blast['send_sort'] = blast_startend[:,1]
            # # Create filter to remove all except the strongest interactinos
            # mch_bool = blast.MCH >= config['mch_filter']
            # tm_bool = blast.TM >= config['tm_filter']
            # gc_bool = blast.GC >= config['gc_filter']
            # # Create filter to ignore on-target blasts
            # target_alignments = pd.read_csv(input[1])
            # on_target_bool = pd.Series(np.zeros(gc_bool.shape))
            # for i, r in target_alignments.iterrows():
            #     ot_sseqid = (blast.sseqid == r.sseqid)
            #     # only ignore on target blasts to the same strand as target
            #     ot_sstrand = (blast.sstrand == r.sstrand)
            #     if r.aln_sub: # if only part of the db molecule is the target
            #         a_start, a_end = np.sort([r.sstart, r.send])
            #         otb_start = (blast.sstart_sort >= a_start) & (blast.sstart_sort <= a_end)
            #         otb_end = (blast.send_sort >= a_start) & (blast.send_sort <= a_end)
            #         ot_startend = (otb_start | otb_end)
            #     on_target_bool += ot_sseqid & ot_sstrand & ot_startend
            # off_target_bool = np.invert(on_target_bool.astype(bool))
            # # Apply boolean filters
            # blast_filtered = blast.loc[mch_bool & tm_bool & gc_bool & off_target_bool, :]
        else:
            blast['MCH'] = ''
            blast['TM'] = ''
            blast['GC'] = ''
            blast_measured = blast
            blast_filtered = blast
        # Check the removal of on targets
        should_be_filtered = [blast_filtered[blast_filtered.sseqid == t].shape[0] for t in target_alignments.sseqid]
        print('LOOK HERE: ', sum(should_be_filtered))
        # save the resulting file
        blast_measured.to_csv(output[0], index=False)
        blast_filtered.to_csv(output[1], index=False)
