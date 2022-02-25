from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

################################################
# Split target files into individual fasta
################################################

rule prep_target_files:
    input:
        config['input_dir'] + '/{in_file}' + INPUT_EXT
    output:
        config['probe_generate_dir'] + '/{in_file}/target_files/{target}.fasta'
    run:
        out_dir = os.path.split(output[0])[0]
        name = os.path.split(input[0])[1]
        basename, ext = os.path.splitext(name)
        rev_comp = config['design_reverse_complement']
        records = SeqIO.parse(input[0], ext[1:])
        for r in records:
            target = r.id
            seq = r.seq.reverse_complement() if  rev_comp else r.seq
            new_record = SeqRecord(seq, id=target)
            out_fn = out_dir + '/' + target + '.fasta'
            SeqIO.write(new_record, out_fn, 'fasta')
#     shell:
#         """
#         awk '/^>/ {close(OUT); OUT={params.in_dir} "/" {wildcards.in_file} "/" substr($0,2) ".fasta"};OUT {print >OUT}' {input}
#         """
