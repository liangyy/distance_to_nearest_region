# this module takes variant and region file
# output the distance of each variant to nearest region

rule prepare_junction:
    input:
        region = lambda wildcards: config[wildcards.taskname]['region']
    output:
        temp('junction/junction1_{taskname}.bed'),
        temp('junction/junction2_{taskname}.bed')
    shell:
        'zcat {input.region} | awk -F"\\t" -v OFS="\\t" '{{print $1,$2,$2+1,$4,$5,$6,1}} > {output[0]}'
        'zcat {input.region} | awk -F"\\t" -v OFS="\\t" '{{print $1,$3-1,$3,$4,$5,$6,2}} > {output[1]}'

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
        lambda wildcards: config[wildcards.taskname]['prefiltering_window_size'],
        lambda wildcards: config[wildcards.taskname]['genome_size']
    output:
        'junction/window__{taskname}.bed.gz'
    shell:
        'bedtools slop -i {input[0]} -g {params[1]} -b {params[0]} | gzip > {output[0]}'

rule prefiltering:
    input:
        varaint = lambda wildcards: config[wildcards.taskname]['variant'][wildcards.dataname],
        region = 'junction/window__{taskname}.bed.gz'
    output:
        temp('temp/intersect__{taskname}__{dataname}.bed.gz')
    shell:
        'bedtools intersect -a <( \
            cat {input.variant}| awk -F"\\t" -v OPS="\\t" '{{print $1,$2,$3}}') \
            -b <(zcat {input.region} | awk -F"\\t" -v OPS="\\t" '{{print $1,$2,$3}}') \
            -wao | gzip > {output[0]}'

rule get_distance:
    input:
        'temp/intersect__{taskname}__{dataname}.bed.gz'
    output:
        'output/distance__{taskname}__{dataname}.tab.gz'
    shell:
        'Rscript bed2distance.R --intersection {input[0]} --output {output[0]}'
