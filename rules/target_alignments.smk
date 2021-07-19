from scripts.blast_functions import blastTargets
import pandas as pd


# Blast each probe
rule target_alignments:
    input:
        config['output_dir'] + '/{in_file}/target_files/{target}.fasta'
    output:
        config['output_dir'] + '/{in_file}/target_alignment/{target}.fasta.blast.out',
        config['output_dir'] + '/{in_file}/target_alignment/{target}.csv'
    run:
        output_dir = os.path.split(output[0])[0]
        blastTargets(query=input[0], database=config['blast_database'],
                    output_dir=output_dir, strand=config['blast_strand'])
        blast_output = pd.read_csv(output[0], sep='\t', header=None)
        blast_output.columns = ['qseqid', 'sseqid', 'pident', 'qcovhsp', 'length',
                                'mismatch', 'gapopen', 'qstart', 'qend', 'sstart',
                                'send', 'evalue', 'bitscore', 'staxids']
        target_alignments = blast_output.loc[blast_output.qcovhsp >=
                                             config['target_alignment_qcovhsp'],
                                         ['sseqid','sstart','send']]
        target_alignments.to_csv(output[1], index=False)
