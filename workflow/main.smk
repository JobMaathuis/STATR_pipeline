configfile: 'config/config.yaml'

rule all:
    input:
        'all.done',
        expand(config['wdir'] + config['result'] + '4.meta_analysis/{sample}.txt', sample=config['samples']),
        expand(config['wdir'] + config['result'] + '4.meta_analysis/{sample}_profile.gff', sample=config['samples'])

include: config['wdir'] + config['rules'] + 'find_degs' + config['smk']
include: config['wdir'] + config['rules'] + 'create_ribosome_profile' + config['smk']
include: config['wdir'] + config['rules'] + 'decompile_alignment' + config['smk']
include: config['wdir'] + config['rules'] + 'trim_reads' + config['smk']
include: config['wdir'] + config['rules'] + 'align_reads' + config['smk']


