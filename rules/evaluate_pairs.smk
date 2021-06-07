from scripts.pair_functions import get_bad_pair_bool
import pandas as pd


# Check whether or not probe pairs have off target interactions near each other
rule evaluate_pairs:
    input:
        config['output_dir'] + '/{in_file}/all_pairs_index/{target}.csv',
        aggregate_filtered_probe_blasts,
    output:
        config['output_dir'] + '/{in_file}/evaluated_pairs/{target}.csv',
    run:
        target = os.path.splitext(os.path.split(output[0])[1])[0]
        print('Evaluating probe pairs for target: ', target)
        pairs_df = pd.read_csv(input[0])
        probe_blast_dir = os.path.split(input[1])[0]
        bad_pair_bool = []
        left_ot_count = []
        right_ot_count = []
        combined_ot_count = []
        # iterate through pairs
        for index, (left, right) in pairs_df.iterrows():
            print('pair ', index , '\t\tof ', pairs_df.shape[0], end='\r')
            # get blast info on the  pair
            blast_filenames = [probe_blast_dir + '/' + str(int(pid)) + '.csv' 
                               for pid in [left, right]]
            blasts = [pd.read_csv(fn) for fn in blast_filenames]
            # decide if the pair is good or not
            min_diff = config['minimum_bp_between_off_target_sites']
            bad_pair_bool.append(get_bad_pair_bool(blasts, min_diff))
            # Count off target blasts for each and combined
            left_ot_count.append(blasts[0].shape[0])
            right_ot_count.append(blasts[1].shape[0])
            combined_ot_count.append(blasts[0].shape[0] + blasts[1].shape[0])
        pairs_df['bad_pair_bool'] = bad_pair_bool
        pairs_df['left_ot_count'] = left_ot_count
        pairs_df['right_ot_count'] = right_ot_count
        pairs_df['combined_ot_count'] = combined_ot_count
        pairs_df.to_csv(output[0], index=False)
            
                            
                
                
