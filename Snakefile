# this module takes variant and region file
# output the distance of each variant to nearest region

rule prepare_junction:
    input:
        region = lambda wildcards: config[wildcards.taskname]['region']
    output:
        temp('junction/junction1_{taskname}.bed'),
        temp('junction/junction2_{taskname}.bed')
    shell:
        '''zcat {input.region} | awk -F"\\t" -v OFS="\\t" '{{print $1,$2,$2+1,$4,$5,$6,1}}' > {output[0]}; \
        zcat {input.region} | awk -F"\\t" -v OFS="\\t" '{{print $1,$3-1,$3,$4,$5,$6,2}}' > {output[1]}'''

rule junction_merge:
    input:
        'junction/junction1_{taskname}.bed',
        'junction/junction2_{taskname}.bed'
    output:
        'junction/merged_{taskname}.bed.gz'
    shell:
        'bedtools sort -i <(cat {input[0]} {input[1]}) | gzip > {output[0]}'

rule extend_junction:
    input:
        'junction/merged_{taskname}.bed.gz'
    params:
        lambda wildcards: config[wildcards.taskname]['params']['prefiltering_window_size'],
        lambda wildcards: config[wildcards.taskname]['params']['genome_size']
    output:
        'junction/window__{taskname}.bed.gz'
    shell:
        'bedtools slop -i {input[0]} -g {params[1]} -b {params[0]} | gzip > {output[0]}'

rule prefiltering:
    input:
        variant = lambda wildcards: config[wildcards.taskname]['variant'][wildcards.dataname],
        region = 'junction/window__{taskname}.bed.gz'
    output:
        temp('temp/intersect__{taskname}__{dataname}.bed.gz')
    shell:
        '''bedtools intersect -a <( \
            cat {input.variant}| awk -F"\\t" -v OFS="\\t" '{{print $1,$2,$3}}') \
            -b <(zcat {input.region} | awk -F"\\t" -v OFS="\\t" '{{print $1,$2,$3}}') \
            -wao | gzip > {output[0]}'''

def get_all_split_files(chrms):
    chrms = chrms.split(',')
    out = []
    for chrm in chrms:
        out.append('temp/intersect__{{taskname}}__{{dataname}}__{chrm}.bed.gz'.format(chrm = chrm))
    return out

rule split_by_chrm:
    input:
        'temp/intersect__{taskname}__{dataname}.bed.gz'
    params:
        lambda wildcards: config[wildcards.taskname]['params']['chromosomes']
    output:
        lambda wildcards: get_all_split_files(config[wildcards.taskname]['params']['chromosomes'])
    shell:
        'python scripts/split_by_chrm.py --input {input[0]} \
            --chrm {params[0]} \
            --prefix temp/intersect__{wildcards.taskname}__{wildcards.dataname}__ \
            --suffix .bed.gz'

rule get_distance:
    input:
        'temp/intersect__{taskname}__{dataname}__{chrm}.bed.gz'
    output:
        'output/distance__{taskname}__{dataname}__{chrm}.tab.gz'
    shell:
        'Rscript scripts/bed2distance.R --intersection {input[0]} --output {output[0]}'

rule merge_chrm:
    input:
        lambda wildcards: get_all_split_files(config[wildcards.taskname]['params']['chromosomes'])
    output:
        'output/distance__{taskname}__{dataname}.tab.gz'
    shell:
        'bedtools sort -i <(zcat {input}) | gzip > {output[0]}'
