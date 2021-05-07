import pandas as pd
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

########################################################################
# Parse the selected probe pairs, add flanking regions and write an idt 
    # order sheet excel file
########################################################################

rule append_flanking:
    input:
        [['{out_dir}/{in_file}/selected_pairs/{target}.csv'\
          .format(out_dir=config['output_dir'], in_file=in_file, target=target) 
          for target in targets]
         for in_file, targets in zip(INPUT_BASENAMES, TARGETS)],
    output:
        config['output_dir'] + '/{in_file}/final_outputs/selection.fasta',
    run:
        ids = []
        sequences = []
        for fn in input:
            pairs = pd.read_csv(fn)
            # add flanking regions to sequences
            target = os.path.splitext(os.path.split(fn)[1])[0]
            flanking_id = config['target_flanking'][target]
            flanking_df = pd.read_csv(config['flanking_regions_filename'])
            flank_left, flank_right, flank_id= flanking_df.loc[flanking_df.ID == flanking_id, 
                                                               ['left_flanking', 'right_flanking', 'ID']
                                                               ].values[0].astype(list)
            pairs['full_seq_left'] = flank_left + config['spacer'] + pairs.seq_left
            pairs['full_seq_right'] = pairs.seq_right + config['spacer'] + flank_right
            # Add the sequences and their ids to the output lists
            target_name = config['target_naming'][target]
            starts = pairs.start.values
            l_seqs = pairs.full_seq_left.values
            r_seqs = pairs.full_seq_right.values
            for st, l_sq, r_sq in zip(starts, l_seqs, r_seqs):
                for side, sq in zip(['L','R'], [l_sq, r_sq]):
                    ids.append('sp_' +  flank_id + '.' + target_name + '.' + str(st) + '.' + side)
                    sequences.append(sq)
        # write fasta 
        records = [SeqRecord(Seq(seq),id=i) for seq, i in zip(sequences, ids)]
        SeqIO.write(records, output[0], 'fasta')
