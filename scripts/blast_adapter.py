import subprocess
import os

def blastProbes(query, database, output_dir, strand='minus', blast_extension='.blast.out'):
    basename = os.path.split(query)[1]
    output = output_dir + '/' + basename + blast_extension
    out_format = out_format = '6 qseqid sseqid pident qcovhsp length mismatch \
                               gapopen qstart qend sstart send evalue bitscore \
                               staxids qseq sseq'
    command_line_blast = ['blastn', '-db', database, '-query', query, 
                          '-out', output, '-outfmt', out_format, '-task',
                          'blastn-short', '-max_hsps', '1', '-max_target_seqs',
                          '100000', '-strand', strand, '-evalue', '100', 
                          '-num_threads', '1']       
    subprocess.call(command_line_blast)
    return

 

