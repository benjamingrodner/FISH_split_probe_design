# functions for filtering blast results

from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import GC as gc 
import pandas as pd
import numpy as np
from Bio.Seq import Seq
import os
import subprocess

def blastProbes(query, database, output_dir, strand='both', blast_extension='.blast.out'):
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


def read_blast_table(filename):
    blast = pd.read_csv(filename, sep='\t')
    blast.columns = ['qseqid', 'sseqid', 'pident', 'qcovhsp', 'length', 'mismatch',
                     'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 
                     'bitscore', 'staxids', 'qseq', 'sseq']
    return blast
    

def max_continuous_homology(qseq, sseq):
    if qseq != sseq:
        snp_indices = np.where(np.array(list(qseq)) != np.array(list(sseq)))[0]
        diffs = np.diff(snp_indices)
        mch = np.max(np.append(diffs,[snp_indices[0], len(qseq) - 1 - 
                                               snp_indices[-1]]))
    else:
        mch = len(qseq)
    return mch


def melting_temperature(qseq, sseq, Na = 390, dnac1_oligo = 100):
    sseq_template = str(Seq(str(sseq)).reverse_complement()[::-1])
    try:
        tm = mt.Tm_NN(qseq, c_seq=sseq_template, Na = Na, saltcorr = 7, 
                      dnac1 = dnac1_oligo, dnac2 = 1)
    except:
        tm = mt.Tm_NN(qseq, Na = Na, saltcorr = 7, 
                      dnac1 = dnac1_oligo, dnac2 = 1)
    return tm


def gc_content(qseq, length):
    return gc(qseq)/100 * length


# def _check_each(reference_blast, sseqid, sstart, send, ):
#     return 

# def _check_on_target(sseqid, sstart, send, on_target_blast):
#     ont_slice = on_target_blast.loc[:,['sseqid','sstart','send']]
#     bool_list = []
#     for index, (ot_sseqid, ot_sstart, ot_send) in ont_slice.iterrows():
#         bool_sseqid = sseq_id == ot_sseqid
#         bool_sstart = (sstart > ot_sstart) & (sstart < ot_send)
#         bool_send =  (send > ot_sstart) & (send < ot_send)
#         if bool_sseqid & bool_sstart & bool_send:
#             bool_list.append(1)  
#     return any(bool_list)

# def filter_on_target_blasts(blast, on_target_blast):
#     blast['on_target'] = blast.apply(lambda x: _check_on_target(x.sseqid, x.sstart, x.send, 
#                                                             on_target_blast), axis=1)
#     return blast[blast.on_target == 0]
    
    
    