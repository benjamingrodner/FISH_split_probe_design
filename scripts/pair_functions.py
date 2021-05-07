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
    for index, (left, right, bad_pair_bool) in pairs.iterrows():
        s_left, sq_left = primer3_table.loc[primer3_table.probe_num == left, 
                                      ['start','seq']].values[0].astype(list)
        s_right, l_right, sq_right = primer3_table.loc[primer3_table.probe_num == right, 
                                                 ['start', 'length', 'seq']]\
                                                .values[0].astype(list)
        starts.append(s_left)
        ends.append(s_right + l_right)
        seqs_left.append(sq_left)
        seqs_right.append(sq_right)
    pairs['start'] = starts
    pairs['end'] = ends
    pairs['seq_right'] = seqs_right
    pairs['seq_left'] = seqs_left
    return pairs


def get_probesets(repetitions, pairs):
    probesets = []
    for rep in range(repetitions):
        p_set = pd.DataFrame([], columns=pairs.columns)
        # Keep picking pairs and filtering until none are left 
        pairs_ = pairs.copy()
        while pairs_.shape[0] > 0:
            # Pick a pair at random
            index = random.randint(0, pairs_.shape[0]-1)
            selection = pairs_.iloc[index,:]
            p_set = p_set.append(selection, ignore_index=True)
            # Filter overlapping pairs
            bool_left = pairs_.end < selection.start
            bool_right = pairs_.start > selection.end
            pairs_ = pairs_[bool_left | bool_right].reset_index(drop=True)
        probesets.append(p_set)
    return probesets
