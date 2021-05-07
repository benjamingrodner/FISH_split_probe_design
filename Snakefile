##################
# Split-probe design
# Ben Grodner 5/5/21
##################
from Bio import SeqIO


##################
# Functions
##################

def get_targets(input_filenames):
    targets = []
    for fn in input_filenames:
        ext = os.path.splitext(fn)[1][1:]
        records = SeqIO.parse(fn, ext)
        ids = [r.id for r in records]
        targets.append(ids)
    return targets


def aggregate_filtered_probe_blasts(wildcards):
    checkpoint_output = checkpoints.prep_probe_files.get(**wildcards).output[0]
    return expand(config['output_dir'] + 
                  "/{in_file}/blast/{target}/filtered/{pid}.csv",
                  in_file=wildcards.in_file, target=wildcards.target,
                  pid=glob_wildcards(os.path.join(checkpoint_output, 
                                                  "{pid}.fasta")).pid)


##################
# Parameters
##################

INPUT_FILENAMES_FULL = [config['input_dir'] + '/' +  in_file 
                        for in_file in config['input_filenames']]
INPUT_BASENAMES = [os.path.splitext(fn)[0] for fn in config['input_filenames']]
INPUT_EXT = os.path.splitext(INPUT_FILENAMES_FULL[0])[1]
TARGETS = get_targets(INPUT_FILENAMES_FULL)

##################
# Snake rules
##################
# ruleorder: aggregate > blast_probes > prep_probe_files
rule all:
    input:
#     # first inputs
#         INPUT_FILENAMES_FULL,
#     # prep_target_files
#         [['{}/{}/target_files/{}.fasta'.format(config['output_dir'], in_file, target) 
#           for target in targets] 
#          for in_file, targets in zip(INPUT_BASENAMES, TARGETS)],
#     # design probes
#         [['{out_dir}/{in_file}/primer3/{target}/{target}.fasta'\
#           .format(out_dir=config['output_dir'], in_file=in_file, target=target)
#           for target in targets]
#          for in_file, targets in zip(INPUT_BASENAMES, TARGETS)],
#    # make probe pairs
#         [['{}/{}/all_pairs_index/{}.csv'.format(config['output_dir'], in_file, target) 
#           for target in targets] 
#          for in_file, targets in zip(INPUT_BASENAMES, TARGETS)],
#    # prep probe files
#         [['{}/{}/blast/{}/probes'.format(config['output_dir'], in_file, target) 
#           for target in targets]
#          for in_file, targets in zip(INPUT_BASENAMES, TARGETS)],
#    # blast_probes
#         [['{}/{}/blast/aggregates/{}.txt'.format(config['output_dir'], in_file, target) 
#           for target in targets]
#          for in_file, targets in zip(INPUT_BASENAMES, TARGETS)],
#    # evaluate_probes
# #         'outputs/arcobacter_16s/evaluated_pairs/NZ_CP032099.1:c1898141-1896624.csv'
#         [['{}/{}/evaluated_pairs/{}.csv'.format(config['output_dir'], in_file, target) 
#           for target in targets]
#          for in_file, targets in zip(INPUT_BASENAMES, TARGETS)],
#     # pick_pairs
#         [['{out_dir}/{in_file}/selected_pairs/{target}.csv'.format(out_dir=config['output_dir'], 
#                                                                    in_file=in_file, target=target) 
#           for target in targets]
#          for in_file, targets in zip(INPUT_BASENAMES, TARGETS)],
    # get idt order sheet
        [config['output_dir'] + '/{in_file}/final_outputs/order_sheet.xlsx'\
         .format(in_file=in_file) for in_file in INPUT_BASENAMES],
        
include: 'rules/prep_target_files.smk'
include: 'rules/design_probes.smk'
include: 'rules/make_pairs.smk'
include: 'rules/prep_probe_files.smk'
include: 'rules/blast_probes.smk'
include: 'rules/filter_blasts.smk'
include: 'rules/aggregate_probes.smk'
# include: 'rules/evaluate_pairs.smk'
include: 'rules/select_pairs.smk'
include: 'rules/append_flanking.smk'
include: 'rules/idt_order_sheet.smk'
