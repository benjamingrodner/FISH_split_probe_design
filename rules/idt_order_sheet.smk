
########################################################################
# write an idt order sheet excel file
########################################################################

rule idt_order_sheet:
    input:
        config['output_dir'] + '/{in_file}/final_outputs/selection.fasta',
    output:
        config['output_dir'] + '/{in_file}/final_outputs/order_sheet.xlsx',
    params:
        scripts_path = config['pipeline_dir']
#     conda:
#         '../conda_envs/openpyxl.yaml'
    shell:
        "python {params.scripts_path}/scripts/write_idt_order_sheet.py "
        "-in {input} -out {output}"
