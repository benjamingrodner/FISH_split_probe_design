import re
import os
import pandas as pd

rule select_flanking_spacers:
    input:
        config['probe_select_dir'] + '/{in_file}/selected_pairs/{target}.csv',
        config['spacer_select_dir'] + '/{in_file}/blast/{target}/inputs',
        aggregate_filtered_spacer_blasts,
    output:
        config['spacer_select_dir'] + '/{in_file}/spacer_selection/{target}.csv',
        config['final_outputs_dir'] + '/{in_file}/selection_table/{target}.csv'
    run:
        selected_pairs = pd.read_csv(input[0])
        # Count strong blasts and combine into a single df for each selected probe
        l_results_list = []
        r_results_list = []
        start_list = selected_pairs.start.values
        cols = ['spacer','seq_full','ot_count','max_ot_gc','max_ot_mch','max_ot_tm']
        for start in start_list:
            l_result_df = pd.DataFrame(columns=cols)
            r_result_df = pd.DataFrame(columns=cols)
            for filtered_fn in input[2:]:
                spacer = os.path.splitext(os.path.split(filtered_fn)[1])[0]
                fasta_fn = input[1] + '/' + spacer + '.fasta'
                if str(int(start)) in spacer:
                    blast_filtered = pd.read_csv(filtered_fn)
                    record = SeqIO.parse(fasta_fn, 'fasta')
                    seq = [r.seq for r in record][0]
                    ot_count = blast_filtered.shape[0]
                    if ot_count:
                        max_ot_gc = blast_filtered.GC.max()
                        max_ot_mch = blast_filtered.MCH.max()
                        max_ot_tm = blast_filtered.TM.max()
                    else:
                        max_ot_gc, max_ot_mch, max_ot_tm = 0,0,0
                    blast_summary = pd.Series([spacer, str(seq), ot_count, max_ot_gc, max_ot_mch, max_ot_tm], index=cols)
                    if 'L' in spacer:
                        l_result_df = l_result_df.append(blast_summary, ignore_index=True)
                    else:
                        r_result_df = r_result_df.append(blast_summary, ignore_index=True)
            l_results_list.append(l_result_df)
            r_results_list.append(r_result_df)
        # Pick the best spacer for each probe
        series = []
        for start, left_results, right_results in zip(start_list, l_results_list, r_results_list):
            sp_left = left_results.sort_values(by='ot_count', ascending=True).reset_index(drop=True).loc[0,:]
            sp_left.index = [n + '_l' for n in sp_left.index]
            sp_right = right_results.sort_values(by='ot_count', ascending=True).reset_index(drop=True).loc[0,:]
            sp_right.index = [n + '_r' for n in sp_right.index]
            combined = sp_left.append(sp_right)
            combined['start'] = start
            series.append(combined)
        spacer_selection = pd.DataFrame(series)
        selection_table_full = selected_pairs.merge(spacer_selection, on='start', how='left')
        starts = selection_table_full.start.astype(int).astype(str)
        flanking_id = config['target_flanking'][wildcards.target]
        target_name = config['target_naming'][wildcards.target]
        selection_table_full['probe_name_l'] = 'sp_' + flanking_id + '.' + target_name + '_' + starts + '.L'
        selection_table_full['probe_name_r'] = 'sp_' + flanking_id + '.' + target_name + '_' + starts + '.R'
        selection_table_full.to_csv(output[0])
        selection_table_full.to_csv(output[1])
