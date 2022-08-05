import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from scripts import fncs
import re
import subprocess
import glob
import pandas as pd
from Bio.SeqUtils import MeltingTemp as mt

##############################################################################
# Primer 3 class
##############################################################################
'''
Input a fasta file with the target sequences to the Primer3 object
Use the design_probes function to generate a csv table with columns
probe_num', 'seq', 'start', 'length', 'N', 'GC', 'Tm', 'self_any_th', 'self_end_th', 'hairpin', 'quality'
'''

class Primer3Error(Exception):
    pass


class Primer3(object):
    def __init__(self, fasta_filename, primer3_path
                 min_tm=50, Na=390, dnac1_oligo=100, 
                 min_size=18, opt_size=20, max_size=22):
        self.fasta_filename = fasta_filename
        self.primer3_path = primer3_path
        self.targets = [record.id for record in SeqIO.parse(fasta_filename, 'fasta')]
        self.original_dir = os.getcwd()
        self.min_tm = min_tm
        self.Na = Na
        self.dnac1_oligo = dnac1_oligo
        self.min_size = min_size
        self.opt_size = opt_size
        self.max_size = max_size

    def _make_fasta_output_dir(self, primer3_dir):
        p3_dir = os.path.abspath(fncs.remove_slash(primer3_dir))
        basename = os.path.basename(self.fasta_filename)
        basename = re.sub('.fasta', '', basename)
        self.output_dir = '{}/{}'.format(p3_dir, basename)
        fncs.confirm_dir_exists(self.output_dir)

    def _write_primer3_settings_file(self):
        primer3_settings = ['Primer3 File - http://primer3.sourceforge.net',
                            'P3_FILE_TYPE=settings',
                            '',
                            'P3_FILE_ID=FISH probe design',
                            'P3_FILE_FLAG=1',
                            'PRIMER_FIRST_BASE_INDEX=1',
                            'PRIMER_TASK=generic',
                            'PRIMER_EXPLAIN_FLAG=1',
                            'PRIMER_NUM_RETURN=10000',
                            'PRIMER_PICK_LEFT_PRIMER=0',
                            'PRIMER_PICK_INTERNAL_OLIGO=1',
                            'PRIMER_PICK_RIGHT_PRIMER=0',
                            'PRIMER_INTERNAL_OPT_SIZE=' + str(self.opt_size),
                            'PRIMER_INTERNAL_MIN_SIZE='+ str(self.min_size),
                            'PRIMER_INTERNAL_MAX_SIZE=' + str(self.max_size),
                            'PRIMER_INTERNAL_MIN_TM=' + str(self.min_tm),
                            'PRIMER_INTERNAL_MAX_SELF_ANY_TH=1000.00',
                            'PRIMER_INTERNAL_MAX_HAIRPIN_TH=1000.0',
                            'PRIMER_INTERNAL_MAX_NS_ACCEPTED=0',
                            'PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=1',
                            'PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT=0',
                            'PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/programs/primer3-2.3.5/src/primer3_config/',
                            'PRIMER_LOWERCASE_MASKING=0',
                            'PRIMER_PICK_ANYWAY=1',
                            '=']
        if not os.path.exists(self.primer3_settings_filename):
            with open(self.primer3_settings_filename, 'a') as file:
                for line in primer3_settings:
                    file.write(line + '\n')
            print('Wrote: ', self.primer3_settings_filename)
        return

    def _write_input(self):
        open(self.primer3_input, 'w').close()
        for seq_record in SeqIO.parse(self.fasta_filename, "fasta"):
            with open(self.primer3_input, 'a') as f:
                f.write('PRIMER_EXPLAIN_FLAG=1\n')
                f.write('P3_FILE_FLAG=1\n')
                f.write('SEQUENCE_ID=' + seq_record.id + '\n')
                f.write('SEQUENCE_TEMPLATE=' + str(seq_record.seq) + '\n')
                f.write('=\n')
        print('Wrote primer3 input: ', self.primer3_input)

    def _run_primer3(self):
        subprocess.check_call([self.primer3_path, '-p3_settings_file', self.primer3_settings_filename, '-output', self.primer3_output, '-format_output', self.primer3_input])
        return

    def _get_probes_as_df(self, probe_filename):
        probes = pd.read_csv(probe_filename, skiprows=3, header=None, sep='\s+')
        probes.columns = ['probe_num', 'seq', 'start', 'length', 'N', 'GC', 'Tm', 'self_any_th', 'self_end_th', 'hairpin', 'quality']
        return probes

    def get_probe_dfs(self):
        return [pd.read_csv(fn) for fn in self.csv_filenames]


    def _Tm(self, sequence):
        return mt.Tm_NN(sequence, Na = self.Na, saltcorr = 7,
                        dnac1 = self.dnac1_oligo, dnac2 = 1)

    def _write_csv_dataframes(self):
        int_filenames = glob.glob(self.output_dir + '/*.int')
        self.csv_filenames = []
        for int_f in int_filenames:
            df = self._get_probes_as_df(int_f)
            # Generate more accurate melting temperature values
            df['Tm_NN'] = df['seq'].apply(self._Tm)
            csv_filename = re.sub('\.int','.csv', int_f)
            df.to_csv(csv_filename, index=False)
            self.csv_filenames.append(csv_filename)

    def _write_probe_fastas(self):
        self.probe_fasta_filenames = []
        for csv_f in self.csv_filenames:
            csv = pd.read_csv(csv_f)
            seqs = csv.loc[:,['seq','probe_num']].values
            records = [SeqRecord(Seq(seq[0]), id=str(seq[1])) for seq in seqs]
            probe_fasta_filename = re.sub('.csv','.fasta', csv_f)
            SeqIO.write(records, probe_fasta_filename, 'fasta')
            self.probe_fasta_filenames.append(probe_fasta_filename)

    def design_probes(self, output_dir):
        self._make_fasta_output_dir(output_dir)
        self.primer3_settings_filename = self.output_dir + '/primer3_settings.txt'
        self._write_primer3_settings_file()
        self.primer3_input = self.output_dir + '/primer3_input.txt'
        self._write_input()
        self.primer3_output = self.output_dir + '/primer3_output.txt'
        os.chdir(self.output_dir)
        self._run_primer3()
        os.chdir(self.original_dir)
        self._write_csv_dataframes()
        self._write_probe_fastas()
        print('Designed probes for: ', self.fasta_filename)
        print('Probes for each feature written in: ', self.output_dir)
