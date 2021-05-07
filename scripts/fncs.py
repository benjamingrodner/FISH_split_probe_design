import os
import pandas as pd
from Bio import SeqIO

# Execute primer3 with a shell command
# It is important that the settings filename is formatted with no space character in the third line
# It is also necessary that the input be a csv file of the format (for each line):
# 'SEQUENCE_ID={id}','SEQUENCE_TEMPLATE={sequence}','SEQUENCE_INCLUDED_REGION={starting base number e.g. 0},{ending base number}', 'P3_FILE_FLAG=1', 'PRIMER_EXPLAIN_FLAG=1','='


def remove_slash(dir):
    return dir[:-1] if dir[len(dir) - 1] == '/' else dir


def confirm_dir_exists(filepath):
    if not os.path.exists(filepath):
        os.makedirs(filepath)
        print('Made dir: ', filepath)

        



