import pandas as pd
# Filter the probes
rule filter_crosstalk:
    input:
        config['crosstalk_dir'] + '/{in_file}/evaluated_crosstalk/{target}.csv',
        config['spacer_select_dir'] + '/{in_file}/spacer_selection/{target}.csv'
    output:
        config['crosstalk_dir'] + '/{in_file}/filtered_crosstalk/{target}.csv'
    run:
        pairs_df = pd.read_csv(input[1])
        crosstalk_df = pd.read_csv(input[0])
        probe_starts = crosstalk_df['Start'].tolist()
        events = crosstalk_df['Crosstalk Events'].tolist()
        for x in range(len(probe_starts)):
            if events[x] > config['maximum_crosstalk_events']:
                pairs_df = pairs_df[pairs_df.start != probe_starts[x]]

        pairs_df.to_csv(output[0])
