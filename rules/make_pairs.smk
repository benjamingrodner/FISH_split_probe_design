from scripts.pair_functions import find_pairs
import pandas as pd

# Make probe pairs that are near each other
rule make_pairs:
    input:
        config['probe_generate_dir'] + '/{in_file}/primer3/{target}/{target}.csv'
    output:
        config['probe_generate_dir'] + '/{in_file}/all_pairs_index/{target}.csv'
    run:
        probe_df = pd.read_csv(input[0])
        starts = probe_df.start.values.astype(list)
        ends = probe_df.loc[:,['start','length']].sum(axis=1).values.astype(list)
        ids = probe_df.probe_num.values.astype(list)
        pairs_index = find_pairs(starts=starts, ends=ends, ids=ids, max_diff=config['maximum_pair_distance_bp'])
        pd.DataFrame(pairs_index, columns=['left','right']).to_csv(output[0], index=False)
