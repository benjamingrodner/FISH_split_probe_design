from scripts.blast_functions import blastTargets, read_target_blast_table
import pandas as pd


# Blast each target sequence against the database to identify on target sequences
rule target_alignments:
    input:
        config['probe_generate_dir'] + '/{in_file}/target_files/{target}.fasta'
    output:
        config['target_alignment_dir'] + '/{in_file}/blast/{target}.fasta.blast.out',
        config['target_alignment_dir'] + '/{in_file}/filtered/{target}.csv'
    run:
        output_dir = os.path.split(output[0])[0]
        if config['target_alignment']['target_in_db']: # Necessary to run the blast?
            # If you did the alignment beforehand you can load your own preformatted results
            if config['on_target_id_filename']:
                otid_fn = config['input_dir'] + '/' + config['on_target_id_filename']
                target_alignments = pd.read_csv(otid_fn, sep='\t')
                with open(output[0],'w') as f:
                    pass
            else:  # Blast against the database
                blastTargets(query=input[0], database=config['blast_database'],
                            output_dir=output_dir, strand=config['blast_strand'])
                try:
                    blast_output = read_target_blast_table(output[0])
                    # Filter results by coverage and identity
                    min_qcovhsp = int(config['target_alignment']['min_qcovhsp'])
                    min_pident = int(config['target_alignment']['min_pident'])
                    bool_cov = blast_output.qcovhsp >= min_qcovhsp
                    bool_pid = blast_output.pident >= min_pident
                    target_alignments = blast_output.loc[bool_cov & bool_pid,
                                                         ['sseqid','sstart','send',
                                                         'sstrand']]
                    target_alignments['aln_sub'] = 1  # only set to zero for manual tar
                except:
                    target_alignments = pd.DataFrame([], columns=['sseqid','sstart',
                                                                  'send','sstrand',
                                                                  'aln_sub'])
        else:
            blast_output = pd.DataFrame([], columns=['qseqid', 'sseqid', 'pident',
                                    'qcovhsp', 'length', 'mismatch', 'gapopen',
                                    'qstart', 'qend', 'sstart','send', 'evalue',
                                    'bitscore', 'staxids'])
            blast_output.to_csv(output[0], index=False)

        target_alignments.to_csv(output[1], index=False)
