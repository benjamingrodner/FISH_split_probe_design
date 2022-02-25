import scripts.blast_functions as bf
import pandas as pd
import numpy as np
import os

########################################################################
# Blast the various spacers plus flanking regions and filter the results
########################################################################

rule blast_flanking_spacers:
    input:
        config['spacer_select_dir'] + '/{in_file}/blast/{target}/inputs/{spacer}.fasta',
        config['target_alignment_dir'] + '/{in_file}/filtered/{target}.csv'
    output:
        config['spacer_select_dir'] + '/{in_file}/blast/{target}/outputs/{spacer}.fasta.blast.out',
        config['spacer_select_dir'] + '/{in_file}/blast/{target}/filtered/{spacer}.csv'
    run:
        output_dir = os.path.split(output[0])[0]
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        bf.blastProbes(query=input[0], database=config['blast_database'], output_dir=output_dir)
        blast = bf.read_blast_table(output[0])
        if blast.shape[0] > 0:
            # calculate mch, tm and gc contents
            # calculate mch, tm and gc contents
            blast = bf.measure_blasts(blast)
            # blast['MCH'] = blast.apply(lambda x: bf.max_continuous_homology(x.qseq, x.sseq), axis=1)
            # blast['TM'] = blast.apply(lambda x: bf.melting_temperature(x.qseq, x.sseq), axis=1)
            # blast['GC'] = blast.apply(lambda x: bf.gc_content(x.qseq, x.length), axis=1)
            # Sort blast start and end
            target_alignments = pd.read_csv(input[1])
            blast_filtered = bf.filter_blasts(blast, target_alignments,
                                              mch_filter=config['mch_filter'],
                                              tm_filter=config['tm_filter'],
                                              gc_filter=config['gc_filter'])
            # blast_startend = blast.loc[:,['sstart','send']].values
            # blast_startend = blast_startend.sort(axis=1)
            # blast['sstart_sort'] = blast_startend[:,0]
            # blast['send_sort'] = blast_startend[:,1]
            # # remove all except the strongest interactinos
            # mch_bool = blast.MCH >= config['mch_filter']
            # tm_bool = blast.TM >= config['tm_filter']
            # gc_bool = blast.GC >= config['gc_filter']
            # # Ignore on-target blasts
            # target_alignments = pd.read_csv(input[1])
            # on_target_bool = pd.Series(np.zeros(gc_bool.shape))
            # for i, r in target_alignments.iterrows():
            #     on_target_bool = (blast.sseqid == r.sseqid)
            #     if r.aln_sub: # if only part of the db molecule is the target
            #         a_start, a_end = np.sort([r.sstart, r.send])
            #         p_start, p_end = [blast.sstart_sort, blast.send_sort]
            #         # p_start, p_end = [blast.sstart, blast.send]
            #         otb_start = (p_start >= a_start) & (p_start <= a_end)
            #         otb_end = (p_end >= a_start) & (p_end <= a_end)
            #         on_target_bool *= (otb_start | otb_end)
            #     # only ignore on target blasts to the same strand as target
            #     on_target_bool *= (blast.sstrand == r.sstrand)
            #
            # # for i, r in target_alignments.iterrows():
            # #     a_start, a_end = np.sort([r.sstart, r.send])
            # #     otb_id = (blast.sseqid == r.sseqid)
            # #     otb_start = (blast.sstart >= a_start) & (blast.sstart <= a_end)
            # #     otb_end = (blast.send >= a_start) & (blast.send <= a_end)
            # #     on_target_bool += otb_id & (otb_start | otb_end)
            #
            # off_target_bool = np.invert(on_target_bool.astype(bool))
            # # Apply boolean filters
            # blast_filtered = blast.loc[mch_bool & tm_bool & gc_bool & off_target_bool, :]
        else:
            blast_filtered=blast
        # save the resulting file
        output_dir = os.path.split(output[1])[0]
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        blast_filtered.to_csv(output[1], index=False)
