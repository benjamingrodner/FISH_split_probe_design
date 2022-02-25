# functions for filtering blast results

from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import GC as gc
import pandas as pd
import numpy as np
from Bio.Seq import Seq
import os
import subprocess


def get_blast_columns(include_seqs=True):
    blast_columns = ['qseqid', 'sseqid', 'pident', 'qcovhsp', 'length',
                     'mismatch','gapopen', 'qstart', 'qend', 'sstart', 'send',
                     'sstrand','evalue', 'bitscore', 'staxids']
    if include_seqs:
        blast_columns += ['qseq','sseq']
    return blast_columns


def blastProbes(query, database, output_dir=False, strand='both', blast_extension='.blast.out'):
    basename = os.path.split(query)[1]
    output_dir = os.path.split(query)[0] if not output_dir else output_dir
    output = output_dir + '/' + basename + blast_extension
    print(output)
    out_format = '6 ' + ' '.join(get_blast_columns())
    # out_format = '6 qseqid sseqid pident qcovhsp length mismatch \
    #                            gapopen qstart qend sstart send sstrand evalue bitscore \
    #                            staxids qseq sseq'
    command_line_blast = ['blastn', '-db', database, '-query', query,
                          '-out', output, '-outfmt', out_format, '-task',
                          'blastn-short', '-max_hsps', '1', '-max_target_seqs',
                          '100000', '-strand', strand, '-evalue', '100',
                          '-num_threads', '1']
    subprocess.call(command_line_blast)
    return


def read_blast_table(filename):
    blast_columns = get_blast_columns()
    # blast_columns = ['qseqid', 'sseqid', 'pident', 'qcovhsp', 'length', 'mismatch',
    #                 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'sstrand','evalue',
    #                 'bitscore', 'staxids', 'qseq', 'sseq']
    try:
        blast = pd.read_csv(filename, sep='\t', header=None)
        blast.columns = blast_columns
    except:
        blast = pd.DataFrame([],columns=blast_columns)
    return blast


def blastTargets(query, database, output_dir, strand='both', blast_extension='.blast.out',
                 task='megablast', threads='1'):
    basename = os.path.split(query)[1]
    output = output_dir + '/' + basename + blast_extension
    out_format = '6 ' + ' '.join(get_blast_columns(include_seqs=False))
    # out_format = out_format = '6 qseqid sseqid pident qcovhsp length mismatch \
    #                            gapopen qstart qend sstart send sstrand evalue bitscore \
    #                            staxids'
    command_line_blast = ['blastn', '-db', database, '-query', query,
                          '-out', output, '-outfmt', out_format, '-task',
                          task, '-max_hsps', '1', '-max_target_seqs',
                          '100000', '-strand', strand, '-evalue', '100',
                          '-num_threads', threads]
    subprocess.call(command_line_blast)
    return


def read_target_blast_table(filename):
    blast_columns = get_blast_columns(include_seqs=False)
    # blast_columns = ['qseqid', 'sseqid', 'pident', 'qcovhsp', 'length', 'mismatch',
    #                 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'sstrand','evalue',
    #                 'bitscore', 'staxids']
    try:
        blast = pd.read_csv(filename, sep='\t', header=None)
        blast.columns = blast_columns
    except:
        blast = pd.DataFrame([],columns=blast_columns)
        print('\n\nLOOK HERE: Target not in database.\n\n')
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


def _sort_blast_startend(blast):
    blast_startend = blast.loc[:,['sstart','send']].values
    blast_startend.sort(axis=1)
    blast['sstart_sort'] = blast_startend[:,0]
    blast['send_sort'] = blast_startend[:,1]
    return blast


def _get_on_target_bool(blast, target_alignments):
    blast = _sort_blast_startend(blast)
    on_target_bool = pd.Series(np.zeros(blast.shape[0]))
    for i, r in target_alignments.iterrows():
        ot_sseqid = (blast.sseqid == r.sseqid)
        # only ignore on target blasts to the same strand as target
        ot_sstrand = (blast.sstrand == r.sstrand)
        if r.aln_sub: # if only part of the db molecule is the target
            a_start, a_end = np.sort([r.sstart, r.send])
            otb_start = (blast.sstart_sort >= a_start) & (blast.sstart_sort <= a_end)
            otb_end = (blast.send_sort >= a_start) & (blast.send_sort <= a_end)
            ot_startend = (otb_start | otb_end)
        else: # if the whole db molecule is the target
            ot_startend = pd.Series(np.ones(blast.shape[0]))
        # To be considered on-target, blast must fulfill all conditions
        on_target_bool += (ot_sseqid & ot_sstrand & ot_startend)
    return np.invert(on_target_bool.astype(bool))

def measure_blasts(blast):
    blast['MCH'] = blast.apply(lambda x: max_continuous_homology(x.qseq, x.sseq), axis=1)
    blast['TM'] = blast.apply(lambda x: melting_temperature(x.qseq, x.sseq), axis=1)
    blast['GC'] = blast.apply(lambda x: gc_content(x.qseq, x.length), axis=1)
    return blast

def filter_blasts(blast,target_alignments, mch_filter, tm_filter, gc_filter):
    # Create filter to remove all except the strongest interactinos
    mch_bool = blast.MCH >= mch_filter
    tm_bool = blast.TM >= tm_filter
    gc_bool = blast.GC >= gc_filter
    # Create filter to ignore on-target blasts
    ot_bool = _get_on_target_bool(blast, target_alignments)
    # Apply boolean filters
    return blast.loc[mch_bool & tm_bool & gc_bool & ot_bool, :]

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
