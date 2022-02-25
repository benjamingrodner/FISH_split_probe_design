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
    return expand(config['probe_select_dir'] +
                  "/{in_file}/blast/{target}/filtered/{pid}.csv",
                  in_file=wildcards.in_file, target=wildcards.target,
                  pid=glob_wildcards(os.path.join(checkpoint_output,
                                                  "{pid}.fasta")).pid)


def aggregate_filtered_spacer_blasts(wildcards):
    checkpoint_output = checkpoints.prep_flanking_spacers.get(**wildcards).output[0]
    return expand(config['spacer_select_dir'] +
                  "/{in_file}/blast/{target}/filtered/{spacer}.csv",
                  in_file=wildcards.in_file, target=wildcards.target,
                  spacer=glob_wildcards(os.path.join(checkpoint_output,
                                                  "{spacer}.fasta")).spacer)

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
    # aggregate spacers
#         [[config['output_dir'] + '/{in_file}/spacer_selection/{target}/aggregate.txt'\
#          .format(in_file=in_file, target=target)
#          for target in targets]
#          for in_file, targets in zip(INPUT_BASENAMES, TARGETS)],
#     # select flanking spacers
#         [[config['output_dir'] + '/{in_file}/spacer_selection/{target}/full_length_spacer_selection.csv'\
#           .format(in_file=in_file, target=target)
#          for target in targets]
#          for in_file, targets in zip(INPUT_BASENAMES, TARGETS)],
    # # evaluate crosstalk
    #     [[config['output_dir'] + '/{in_file}/evaluated_crosstalk/{target}/crosstalk_evaluation.csv'.format(in_file=in_file, target=target) for target in targets] for in_file, targets in zip(INPUT_BASENAMES, TARGETS)]
        # [[config['output_dir'] + '/{in_file}/filtered_crosstalk/{target}/filtered_crosstalk.csv'.format(in_file=in_file, target=target) for target in targets] for in_file, targets in zip(INPUT_BASENAMES, TARGETS)]
    # # Get helpers
    #     [[config['output_dir']
    #             + '/{in_file}/helper/selection/{target}.csv'.format(
    #             in_file=in_file, target=target)
    #             for target in targets]
    #             for in_file, targets in zip(INPUT_BASENAMES, TARGETS)]
#     # merge final outputs
#         config['output_dir'] + '/final_outputs/selection.fasta'
    # get idt order sheet
        config['final_outputs_dir'] + '/idt_order_sheet.xlsx'

include: 'rules/prep_target_files.smk'
include: 'rules/target_alignments.smk'
include: 'rules/design_probes.smk'
include: 'rules/make_pairs.smk'
include: 'rules/prep_probe_files.smk'
include: 'rules/blast_probes.smk'
include: 'rules/filter_blasts.smk'
# include: 'rules/aggregate_probes.smk'
include: 'rules/evaluate_pairs.smk'
include: 'rules/select_pairs.smk'
include: 'rules/prep_flanking_spacers.smk'
# include: 'rules/intermediate.smk'
include: 'rules/blast_flanking_spacers.smk'
# include: 'rules/aggregate_spacers.smk'
include: 'rules/select_flanking_spacers.smk'
include: 'rules/evaluate_crosstalk.smk'
include: 'rules/filter_crosstalk.smk'
include: 'rules/design_helpers.smk'
include: 'rules/select_helpers.smk'
include: 'rules/merge_final_outputs.smk'
include: 'rules/idt_order_sheet.smk'
