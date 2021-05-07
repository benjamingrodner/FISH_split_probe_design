import pandas as pd
from Bio import SeqIO
import argparse

def fasta_to_idt_excel(fasta_filename, excel_filename):
    records = SeqIO.parse(fasta_filename, "fasta")
    idt_df = pd.DataFrame()
    for r in records:
        _scale = '100nm' if len(r.seq) > 60 else '25nm'
        _row = pd.DataFrame({'Name': [r.id],
                             'Sequence': [r.seq],
                             'Scale': [_scale],
                             'Purification': ['STD']})
        idt_df = idt_df.append(_row)
    idt_df.to_excel(excel_filename, index=False)
    return 'Wrote: ' + excel_filename

def main():
    parser = argparse.ArgumentParser('')
    parser.add_argument('-in', '--input_fasta', dest = 'input_fasta', 
                        type = str)
    parser.add_argument('-out', '--output_excel', dest = 'output_excel', 
                        type = str)
    args = parser.parse_args()
    
    fasta_to_idt_excel(fasta_filename=args.input_fasta, 
                       excel_filename=args.output_excel)
    
    return

if __name__ == '__main__':
    main()
