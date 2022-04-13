rule parse_genome_annotation:
    """ Parses the genome annotation of the file to a readable and clean format (GFF3) """
    message: 'parsing genome annotation of ' + config['genome']
    input:
        config['wdir'] + config['files'] + config['genome'] + '.gff3'
    log:
        config['wdir'] + 'logs/4.meta_analysis/parse_genome.log'
    output:
        config['wdir'] + config['result'] + '4.meta_analysis/' + config['genome'] + '.gff'
    shell:
        'mkdir -p ' + config['wdir'] + config['result'] + '4.meta_analysis;'
        '(python3 ' + config['wdir'] + config['resources'] + config['py_scripts'] + 'ParseGenomeAnnotation.py -i {input} -o {output}) 2> {log}'


rule calculate_ribosome_density:
    """ Exploits the periodicity in a way that the average ribosome density over a coding DNA sequence can be calculated"""
    message: 'checking periodicity of {wildcards.sample}'
    input:
        bed = config['wdir'] + config['result'] + '3.decompiled/{sample}.bed',
        gff = config['wdir'] + config['result'] + '4.meta_analysis/' + config['genome'] + '.gff'
    log:
        config['wdir'] + 'logs/4.meta_analysis/{sample}_density.log'
    output:
        config['wdir'] + config['result'] + '4.meta_analysis/{sample}.txt'
    shell:
        '(python3 ' + config['wdir'] + config['resources'] + config['py_scripts'] + 'CheckPeriodicity.py -i {input.bed} -a {input.gff} -o {output}) 2> {log}'


rule generate_profile:
    """ Generates a genome-wide Ribo-Seq profile """
    message: 'generating profile of {wildcards.sample}'
    input:
        bed = config['wdir'] + config['result'] + '3.decompiled/{sample}.bed',
        gff = config['wdir'] + config['result'] + '4.meta_analysis/' + config['genome'] + '.gff'
    log:
        config['wdir'] + 'logs/4.meta_analysis/{sample}_profile.log'
    output:
        config['wdir'] + config['result'] + '4.meta_analysis/{sample}_profile.gff'
    shell:
        '(python3 ' + config['wdir'] + config['resources'] + config['py_scripts'] + 'GenerateProfile.py -i {input.bed} -e 5 -a {input.gff} -n 1000000 -o {output}) 2> {log}'