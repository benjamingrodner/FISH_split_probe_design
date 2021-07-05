from scripts.pair_functions import get_bad_pair_bool
import pandas as pd

# Evaluate the crosstalk between left probe of one pair and right probes of all other pairs
rule evaluate_crosstalk:
    input: 
        config['output_dir'] + '/{in_file}/spacer_selection/{target}/full_length_spacer_selection.csv'
    output:
        config['output_dir'] + '/{in_file}/evaluated_crosstalk/{target}/crosstalk_evaluation.csv' 
    run:
        pairs_df = pd.read_csv(input[0])
        right_probes = pairs_df['start'].tolist()
        right_sequences = pairs_df['seq_full_r'].tolist()
        crosstalk = [0] * len(right_probes)
        #Loop through left probes
        for index, row in pairs_df.iterrows():
            Lspacer = row['seq_full_l'][-21:-18]
            # Remove the duplicates
            rights = list(right_probes)
            rights.remove(row['start'])
            r_seqs = list(right_sequences)
            r_seqs.remove(row['seq_full_r'])
            #Loop through right probes
            for x in range(len(rights)):
                Rspacer = r_seqs[x][18:21]
                spacer_fns = [str(rights[x]) + lr + s for lr, s in zip(['_L_','_R_'],[Lspacer, Rspacer])]
                blast_filenames = ['{out_dir}/{in_file}/spacer_selection/{target}/blast/filtered/{spacer}.csv'.format(out_dir=config['output_dir'],
                                                            in_file=wildcards['in_file'],
                                                            target=wildcards['target'],
                                                            spacer=sfn) for sfn in spacer_fns]
#                 blast_filenames = [config['output_dir'] + '/' + wildcards['in_file'] + '/spacer_selection/' + str(row['start']) + '_L_' + Lspacer + '.csv', config['output_dir'] + '/' + wildcards['in_file'] + '/spacer_selection/' + str(rights[x]) + '_R_' + Rspacer + '.csv']
                blasts = [pd.read_csv(fn) for fn in blast_filenames]
                min_diff = config['minimum_bp_between_off_target_sites']
                #Count crosstalk events, left then right
                if get_bad_pair_bool(blasts, min_diff):
                    crosstalk[index] += 1
                    crosstalk[right_probes.index(rights[x])] += 1            
        crosstalk_df = pd.DataFrame(list(zip(right_probes,crosstalk)), columns = ['Start','Crosstalk Events'])
        crosstalk_df.to_csv(output[0])