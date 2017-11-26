# this module takes variant and region file
# output the distance of each variant to nearest region

rule all:
    input:
        lambda config: 'report/{taskname}.html'.format(taskname = list(config.keys()[0])

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

rule split_by_chrm:
    input:
        'temp/intersect__{taskname}__{dataname}.bed.gz'
    params:
        lambda wildcards: wildcards.chrm
    output:
        temp('temp/intersect__{taskname}__{dataname}__{chrm}.bed.gz')
    shell:
        'python scripts/split_by_chrm.py --input {input[0]} \
            --output {output[0]} --chrm {params[0]}'

rule get_distance:
    input:
        'temp/intersect__{taskname}__{dataname}__{chrm}.bed.gz'
    output:
        temp('temp/distance__{taskname}__{dataname}__{chrm}.tab.gz')
    shell:
        'Rscript scripts/bed2distance.R --intersection {input[0]} --output {output[0]}'

rule merge_chrm:
    input:
        lambda wildcards: get_all_split_files(config[wildcards.taskname]['params']['chromosomes'], wildcards.taskname, wildcards.dataname)
        # 'temp/intersect__{taskname}__{dataname}__chr1.bed.gz'
    output:
        'output/distance__{taskname}__{dataname}.tab.gz'
    shell:
        'bedtools sort -i <(zcat {input}) | gzip > {output[0]}'

def get_all_data(taskname, config):
    out = []
    for data in config[taskname]['variant']:
        out.append('output/distance__{taskname}__{dataname}.tab.gz'.format(taskname = taskname, dataname = data))
    return out

def get_all_data_str(taskname, config):
    out = get_all_data(taskname, config)
    return ','.join(out)

rule plot_rmd:
    params:
        lambda wildcards: get_all_data_str(wildcards.taskname, config),
        lambda wildcards: config[wildcards.taskname]['params']['prefiltering_window_size']
    output:
        'report/{taskname}.rmd'
    run:
        rmd = '''---
title: "Histogram of distance to splicing junction"
output:
    html_document:
        number_sections: true
        toc: true
        toc_depth: 3
        toc_float: true
author: Yanyu Liang
date: "`r format(Sys.time(), '%d %B, %Y')`"
---

```{{r, results='asis'}}
library(stringr)
library(pander)
files <- strsplit('{file_str}', ',')[[1]]
headers <- str_match(files, 'distance__(.+)__(.+).tab.gz')
for (i in 1 : nrow(e)){{
  file <- headers[i, 1]
  taskname <- headers[i, 2]
  data <- headers[i, 3]
  cat("#", paste(data, 'in', taskname), "\n")
  distance <- read.table(file, sep = '\t', header = F)
  count <- table(distance$V5)
  names(count) <- c('dist <= {thre_dist}', 'dist > {thre_dist}')
  pander(count)
  hist(distance$V4, main = 'distance to nearest splicing junction')
}}
```
'''.format(file_str = params[0], thre_dist = params[1])
        o = open(output[0], 'w')
        o.write(rmd)
        o.close()

rule plot_html:
    input:
        lambda wildcards: get_all_data(wildcards.taskname, config)
    output:
        'report/{taskname}.html'
    shell:
        '''Rscript -e "rmarkdown::render('{input[0]}')"'''
