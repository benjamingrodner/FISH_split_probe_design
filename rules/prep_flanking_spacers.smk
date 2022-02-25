import pandas as pd
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

########################################################################
# Add spacers and flanking regions to the selected probes and write
    # to individual fastas
########################################################################

checkpoint prep_flanking_spacers:
    input:
        config['probe_select_dir'] + '/{in_file}/selected_pairs/{target}.csv',
        config['probe_select_dir'] + '/{in_file}/selected_pairs/{target}_alternate.csv'
    output:
        directory(config['spacer_select_dir'] + '/{in_file}/blast/{target}/inputs')
    run:
        # Read in selected pairs
        pselect = 1 if config['select_alternate_design'] else 0
        selected_pairs = pd.read_csv(input[pselect])
        # Build spacers
        spacers = []
        bases = ['a','t','c','g']
        for i in bases:
            for j in bases:
                for k in bases:
                    spacers.append(i + j + k)
        # Split seqs into individual rows in df
        seq_left=pd.DataFrame(columns=['start','left', 'seq_left'])
        seq_right=pd.DataFrame(columns=['start','right', 'seq_right'])
        for index, row in selected_pairs.iterrows():
            sl = row[['start','left', 'seq_left']]
            seq_left = seq_left.append(sl)
            sr = row[['start','right', 'seq_right']]
            seq_right = seq_right.append(sr)
        seq_left['hand'] = 'L'
        seq_right['hand'] = 'R'
        cols = ['start','probe_id','seq','hand']
        seq_left.columns = cols
        seq_right.columns = cols
        seq_table = seq_left.append(seq_right)
        # Get flanking regions
        flanking_id = config['target_flanking'][wildcards.target]
        flanking_df = pd.read_csv(config['flanking_regions_filename'])
        flank_left, flank_right, flank_id= flanking_df.loc[flanking_df.ID == flanking_id,
                                                           ['left_flanking', 'right_flanking', 'ID']
                                                           ].values[0].astype(list)
        # Append flanking and write fastas
        if not os.path.exists(output[0]):
            os.makedirs(output[0])
        for start, pid, seq, hand in seq_table.values:
            for sp in spacers:
                if hand=='L':
                    seq_full = seq + sp + flank_left
                else:
                    seq_full = flank_right + sp + seq
                output_fn = output[0] + '/{start}_{hand}_{sp}.fasta'.format(start=int(start), hand=hand, sp=sp)
                id_full = '{start}_{hand}_{sp}'.format(start=int(start), hand=hand, sp=sp)
                record = SeqRecord(Seq(seq_full), id=id_full)
                SeqIO.write(record, output_fn, 'fasta')
