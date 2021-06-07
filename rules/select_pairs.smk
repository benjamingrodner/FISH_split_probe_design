from scripts.pair_functions import get_pair_info, get_probesets
import numpy as np
import pandas as pd

########################################################################
# Select probe pairs that don't have off target interactions with each other
    # and make sure they don't overlap
########################################################################

rule select_pairs:
    input:
        config['output_dir'] + '/{in_file}/evaluated_pairs/{target}.csv',
        config['output_dir'] + '/{in_file}/primer3/{target}/{target}.csv',
    output:
        config['output_dir'] + '/{in_file}/selected_pairs/{target}.csv',
        config['output_dir'] + '/{in_file}/selected_pairs/{target}_alternate.csv',
    run:
        # Get pair evaluation
        pairs = pd.read_csv(input[0])
        # Remove bad pairs
        pairs = pairs[pairs.bad_pair_bool == 0]
        # Get primer3 info
        probes = pd.read_csv(input[1])
        # Concatenate start and end info to pairs
        pairs = get_pair_info(pairs=pairs, primer3_table=probes)
        # repeat the whole selection process a fixed number of times
        repetitions = config['probe_selection']['repetitions']
        probesets = get_probesets(repetitions=repetitions, pairs=pairs)
        # select the largest probeset
        probeset_sizes = [probeset.shape[0] for probeset in probesets]
        selection_index = np.argmax(probeset_sizes)
        probeset_selection = probesets[selection_index]
        alternate_index = np.argsort(probeset_sizes)[-2]
        probeset_alternate = probesets[alternate_index]
#         # Safe the seleciton
#         probeset_selection.to_csv(output[0], index=False)
        # Save the seleciton
        for ps, out in zip([probeset_selection, probeset_alternate], output):
            ps.to_csv(out, index=False)