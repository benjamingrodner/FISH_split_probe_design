import pandas as pd
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

rule merge_final_outputs:
    input:
        probes = [['{out_dir}/{in_file}/filtered_crosstalk/{target}.csv'\
          .format(out_dir=config['crosstalk_dir'], in_file=in_file, target=target)
          for target in targets]
         for in_file, targets in zip(INPUT_BASENAMES, TARGETS)],
         helpers = [['{out_dir}/{in_file}/selection/{target}.csv'\
           .format(out_dir=config['helper_select_dir'], in_file=in_file, target=target)
           for target in targets]
          for in_file, targets in zip(INPUT_BASENAMES, TARGETS)]
    output:
        config['final_outputs_dir'] + '/selection.fasta',
    run:
        records = []
        for fn in input.probes:
            selection_table_full = pd.read_csv(fn)
            for index, row in selection_table_full.iterrows():
                records.append(SeqRecord(Seq(row.seq_full_l), id=row.probe_name_l))
                records.append(SeqRecord(Seq(row.seq_full_r), id=row.probe_name_r))
        # Check if we designed helpers
        bool_help_same = config['get_helpers']['same_strand']
        bool_help_other = config['get_helpers']['other_strand']
        if bool_help_same or bool_help_other:
            # if so, write to the fasta
            for fn in input.helpers:
                try:
                    helper_table = pd.read_csv(fn)
                    for index, row in helper_table.iterrows():
                        records.append(SeqRecord(Seq(row.seq), id=row.helper_name))
                except:
                    pass
        SeqIO.write(records, output[0], 'fasta')
