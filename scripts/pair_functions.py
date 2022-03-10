import numpy as np
import pandas as pd
import random

#############################################
# Functions for dealing with split probe pairs
#############################################

def find_pairs(starts, ends, ids, max_diff=2):
    split_probes = np.empty((0,2))
    max_diff = 2
    for x in range(len(starts)):
        for y in range(len(ends)):
            diff = starts[x] - ends[y]
            if (diff > 0) and (diff < max_diff):
                new_split = np.array([[ids[y], ids[x]]])
                split_probes = np.concatenate((split_probes, new_split), axis=0)
    return split_probes


def get_bad_pair_bool(blasts, min_diff):
    # get a list of seqids that both probes align to
    seqids = [bl.sseqid.values.astype(list) for bl in blasts]
    # compare the list for OTB overlaps
    intersection = list(set(seqids[0]) & set(seqids[1]))
#     print('intersection',len(intersection))
    # for the otb overlaps, get the sstarts/sends
    for inter in intersection:
        locs = [bl.loc[bl.sseqid == inter, ['sstart','send']]
                  for bl in blasts]
        # if you get one site too close together, then give them False bool
        for index, s in locs[0].sstart.items():
            for e in locs[1].loc[~locs[1].index.isin([s]), 'send'].values:
                diff = s - e
                if diff < min_diff:
                    return 1
    return 0


def get_pair_info(pairs, primer3_table):
    starts = []
    ends = []
    seqs_left = []
    seqs_right = []
    mids = []
    for index, row in pairs.iterrows():
        l_start, l_seq = primer3_table.loc[primer3_table.probe_num == row.left,
                                      ['start','seq']].values[0].astype(list)
        r_start, r_length, r_seq = primer3_table.loc[primer3_table.probe_num == row.right,
                                                 ['start', 'length', 'seq']]\
                                                .values[0].astype(list)
        r_end = r_start + r_length
        mid = np.mean([l_start, r_end])
        starts.append(l_start)
        ends.append(r_end)
        mids.append(mid)
        seqs_left.append(l_seq)
        seqs_right.append(r_seq)
    pairs['start'] = starts
    pairs['end'] = ends
    pairs['seq_right'] = seqs_right
    pairs['seq_left'] = seqs_left
    pairs['middle'] = mids
    return pairs


def get_probesets(pairs):
# def get_probesets(repetitions, pairs):
    probesets = []
#     # random method
#     for rep in range(repetitions):
#         p_set = pd.DataFrame([], columns=pairs.columns)
#         # Keep picking pairs and filtering until none are left
#         pairs_ = pairs.copy()
#         while pairs_.shape[0] > 0:
#             # Pick a pair at random
#             index = random.randint(0, pairs_.shape[0]-1)
#             selection = pairs_.iloc[index,:]    # Sort method
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


def get_pair_overlap(blasts):
    bl_ssids = [bl.sseqid.values for bl in blasts]
    return np.intersect1d(bl_ssids[1], bl_ssids[0])
