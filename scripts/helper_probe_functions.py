import pandas as pd

def get_helper_probe_ends(primer3_table):
    starts = []
    ends = []

    return pairs


def get_probesets(repetitions, pairs):
    probesets = []
    for second_sort in ['left_ot_count','right_ot_count']:
        p_set = pd.DataFrame([], columns=pairs.columns)
        pairs_ = pairs.sort_values(by=['combined_ot_count',second_sort], ascending=[True,True])
        while pairs_.shape[0] > 0:
            # pick the next pair
            selection = pairs_.iloc[0,:]
            p_set = p_set.append(selection, ignore_index=True)
            # Filter overlapping pairs
            bool_downstream = pairs_.start > selection.end
            bool_upstream = pairs_.end < selection.start
            pairs_ = pairs_[bool_upstream | bool_downstream].reset_index(drop=True)
        probesets.append(p_set)
    return probesets
