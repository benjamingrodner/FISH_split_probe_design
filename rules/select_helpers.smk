import pandas as pd

rule select_helpers:
    input:
        config['crosstalk_dir'] + '/{in_file}/filtered_crosstalk/{target}.csv',
        config['helper_select_dir'] + '/{in_file}/primer3/{target}/{target}.csv',
        config['helper_select_dir'] + '/{in_file}/primer3/reverse_complement/{target}/{target}.csv'
    output:
        config['helper_select_dir'] + '/{in_file}/selection/{target}.csv'
    run:
        target_name = config['target_naming'][wildcards.target]
        h_set = pd.DataFrame([])
        # Check if we want any helper probes
        if config['get_helpers']['same_strand']:
            # Load the probe selection
            probes = pd.read_csv(input[0])
            # Load the primer3 file
            helpers = pd.read_csv(input[1])
            # get the end values
            helpers['end'] = helpers['start'] + helpers['length']
            # sort probes by length
            helpers_ = helpers.sort_values(by='length', ascending=False).copy()
            # filter helpers overlapping probes
            for index, row in probes.iterrows():
                bool_downstream = helpers_.start > row.end
                bool_upstream = helpers_.end < row.start
                helpers_ = helpers_[
                        bool_upstream | bool_downstream
                        ].reset_index(drop=True)
            # pick non-overlapping helpers
            h_set = pd.DataFrame([], columns=helpers.columns)
            while helpers_.shape[0] > 0:
                selection = helpers_.iloc[0,:]
                h_set = h_set.append(selection, ignore_index=True)
                bool_downstream = helpers_.start > selection.end
                bool_upstream = helpers_.end < selection.start
                helpers_ = helpers_[
                        bool_upstream | bool_downstream
                        ].reset_index(drop=True)
            starts = h_set['start'].astype(int).astype(str)
            h_set['helper_name'] = 'h.' + target_name + '.' + starts
        h_set_rc = pd.DataFrame([])
        # Check if we want helper probes on the other strand
        if config['get_helpers']['other_strand']:
            # load the primer3 file
            helpers_rc = pd.read_csv(input[2])
            # Get end values
            helpers_rc['end'] = helpers_rc['start'] + helpers_rc['length']
            # sort by length
            helpers_rc_ = helpers_rc.sort_values(
                    by='length', ascending=False
                    ).copy()
            # collect non-overlaping probes
            h_set_rc = pd.DataFrame([], columns=helpers.columns)
            while helpers_rc_.shape[0] > 0:
                selection = helpers_rc_.iloc[0,:]
                h_set_rc = h_set_rc.append(selection, ignore_index=True)
                bool_downstream = helpers_rc_.start > selection.end
                bool_upstream = helpers_rc_.end < selection.start
                helpers_rc_ = helpers_rc_[
                        bool_upstream | bool_downstream
                        ].reset_index(drop=True)
            starts = h_set_rc['start'].astype(int).astype(str)
            h_set_rc['helper_name'] = 'h_rc.' + target_name + '.' + starts

        # Check if any helpers were made
        hs = [h_set, h_set_rc]
        indices = np.where([_hs.shape[0] for _hs in hs])[0]
        if indices.shape[0]:
            # if so save them
            h_final = pd.DataFrame([], columns=hs[indices[0]].columns)
            for i in indices:
                h_final = h_final.append(hs[i])
            h_final.to_csv(output[0])
        else:
            # otherwise write an empty file
            open(output[0], 'a').close()
