import pandas as pd
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

rule merge_final_outputs:
    input:
        [['{out_dir}/{in_file}/spacer_selection/{target}/full_length_spacer_selection.csv'\
          .format(out_dir=config['output_dir'], in_file=in_file, target=target) 
          for target in targets]
         for in_file, targets in zip(INPUT_BASENAMES, TARGETS)]
    output:
        config['output_dir'] + '/final_outputs/selection.fasta',
    run:
        records = []      
        for fn in input:
            selection_table_full = pd.read_csv(fn)
            for index, row in selection_table_full.iterrows():
                records.append(SeqRecord(Seq(row.seq_full_l), id=row.probe_name_l))
                records.append(SeqRecord(Seq(row.seq_full_r), id=row.probe_name_r))
        SeqIO.write(records, output[0], 'fasta')