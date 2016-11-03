import os
import time
from pathlib import Path
from snakemake.utils import listfiles, update_config, format
from pyfaidx import Fasta

default_config = {
    'fastqs': '/home/ecl/dev/ensemble/fastqs',
    'bin': '/home/ecl/dev/ensemble/bin',
    'insert_size': 260,
    'metavelvet_kmer': 31,
    'soap_kmer': 31,
    'abyss_kmer': 31,
    'contig_len1': 150,
    'contig_len2': 300,
    'threads': 8,
    'partition_size': 100000
}

if config:
    for key in default_config.keys():
        if key not in config:
            config[key] = default_config[key]
    update_config(default_config, config)
else:
    config = default_config

CWD=Path(os.getcwd())
WD=Path(config['fastqs'])

starttime = int(time.time())

def build_sample_list(data_fp, filename_fmt, excluded):
    """
    Build a list of samples from a data filepath and filename format.
    :param data_fp: a Path to data files
    :param filename_fmt: a string giving wildcards for {sample} and (optionally)
    {rp}, for read pair (e.g. R1 or R2).
    :param exclude: a list of sample names to exclude
    :returns: A dictionary of samples, with sample names as keys
    """
    files = list(listfiles(str(data_fp/filename_fmt)))
    Samples = {t[1]['sample']: {} for t in files if t[1]['sample'] not in excluded}
    for f in files:
        fpath = f[0]
        wcards = f[1]
        if wcards['sample'] in excluded:
            continue
        rp = wcards['rp'] if 'rp' in wcards.keys() else False
        if rp:
            Samples[wcards['sample']][rp] = fpath
            Samples[wcards['sample']]['paired'] = True
        else:
            Samples[wcards['sample']]['file'] = fpath
            Samples[wcards['sample']]['paired'] = False
    return Samples

def filter_fasta(infile, length, label):
    
    in_fa = Fasta(infile)
    n = 0
    if label:
        label = "_"+label
    for seq in in_fa:
        if len(seq) > length:
            outstr = ">Contig_{n}{label}\n{seq}\n".format(
                n=n, label=label, seq=seq)
            n += 1
            yield outstr
        
Samples = build_sample_list(WD, "{sample}_{rp}.fastq", {})

rule all:
    input:
        expand('final_contigs/{sample}-contigs.fa', sample=Samples.keys())

rule _metavelvet:
    output:
        'metavelvet/{sample}/meta-velvetg.contigs.fa'
    input:
        r1 = str(WD/'{sample}_R1.fastq'),
        r2 = str(WD/'{sample}_R2.fastq')
    params:
        out_fp = 'metavelvet/{sample}'
    log:
        expand('logs/metavelvet/{{sample}}-{starttime}.log', starttime=starttime)
    shell:
        """
        {config[bin]}/velveth {params.out_fp} {config[metavelvet_kmer]} \
        -shortPaired -fastq {input.r1} {input.r2} >& {log} && \
        {config[bin]}/velvetg {params.out_fp} -exp_cov auto -ins_length \
        {config[insert_size]} >>& {log} && \
        {config[bin]}/meta-velvetg {params.out_fp} -ins_length \
        {config[insert_size]} >>& {log}
        """

rule _soap_config:
    output:
        temp('soapdenovo/{sample}/soap_config.txt')
    input:
        r1 = str(WD/'{sample}_R1.fastq'),
        r2 = str(WD/'{sample}_R2.fastq')
    run:
        with open(output[0], 'w') as out:
            out.write((
                "max_rd_len=300\n"
                "[LIB]\n"
                "avg_ins={insert_size}\n"
                "nreverse_seq=0\nasm_flags=3\nrank=1\n"
                "q1={r1}\n"
                "q2={r2}\n").format(
                    insert_size=config['insert_size'],
                    r1=input['r1'],
                    r2=input['r2']))
        
rule _soapdenovo:
    output:
        'soapdenovo/{sample}/{sample}.contig.fa'
    input:
        'soapdenovo/{sample}/soap_config.txt'
    params:
        out_fp = 'soapdenovo/{sample}/{sample}'
    log:
        expand('logs/soapdenovo/{{sample}}-{starttime}.log', starttime=starttime)
    shell:
        """
        {config[bin]}/SOAPdenovo-63mer all -K {config[soap_kmer]} \
        -s {input} -R -o {params.out_fp} -p {config[threads]} >& {log} && \
        cp {params.out_fp}.contig {params.out_fp}.contig.fa 
        """

rule _abyss:
    output:
        'abyss/{sample}/{sample}-unitigs.fa'
    input:
        r1 = str(WD/'{sample}_R1.fastq'),
        r2 = str(WD/'{sample}_R2.fastq')
    params:
        out_fp = 'abyss/{sample}'
    log:
        expand(
            '{cwd}/logs/abyss/{{sample}}-{starttime}.log',
            cwd=str(CWD),
            starttime=starttime)
    shell:
        """
        cd {params.out_fp} && {config[bin]}/abyss-pe name={wildcards.sample} \
        k={config[abyss_kmer]} in='{input.r1} {input.r2}' >& {log}
        """

rule _abyss_partition:
    input:
        r1 = str(WD/'{sample}_R1.fastq'),
        r2 = str(WD/'{sample}_R2.fastq')
    output:
        'abyss_partition/{sample}/unitig-combined.fa'
    params:
        out_fp = 'abyss_partition/{sample}'
    log:
        expand(
            '{cwd}/logs/abyss_partition/{{sample}}-{starttime}.log',
            cwd=str(CWD),
            starttime=starttime)
    run:
        touch(log)
        n_records = int(sum(1 for line in open(input['r1']))/4)
        r1_chunks = expand(
            str(CWD/'abyss_partition/{sample}/chunk_{i}_R1.fastq'),
            sample=wildcards['sample'],
            i=range(n_records // config['partition_size'] + 1))
        r1_str = " ".join(r1_chunks)
        r2_chunks = expand(
            str(CWD/'abyss_partition/{sample}/chunk_{i}_R2.fastq'),
            sample=wildcards['sample'],
            i=range(n_records // config['partition_size'] + 1))
        r2_str = " ".join(r2_chunks)
        shell(
            """
            rbt fastq-split {r1_str} < {input.r1} >& /dev/null && \
            rbt fastq-split {r2_str} < {input.r2} >& /dev/null
            """)
        part_names = []
        for r1, r2 in zip(r1_chunks, r2_chunks):
            part_name = wildcards.sample + re.search(
                "chunk_(.+?)_R1.fastq", r1).group(1)
            part_names.append(part_name)
            shell(
                """
                cd {params.out_fp} && {config[bin]}/abyss-pe name={part_name} \
                k={config[abyss_kmer]} in='{r1} {r2}' >>& {log}
                """)
        part_files = " ".join(params['out_fp']+'/{p}-unitigs.fa'.format(p=p) for p in part_names)
        shell("cat {part_files} > {output[0]}")


rule _idba_ud:
    input:
        'paired/{sample}.fasta'
    output:
        'idba_ud/{sample}/contig.fa'
    params:
        out_fp = 'idba_ud/{sample}'
    log:
        expand('logs/idba_ud/{{sample}}-{starttime}.log', starttime=starttime)
    shell:
        """
        idba_ud -l {input} --num_threads {config[threads]} \
        --out {params.out_fp} &> {log} ||
        if [ ! -a {output} ]; then cp {params.out_fp}/contig-100.fa {output}; fi
        """

rule _pair:
    input:
        r1 = str(WD/'{sample}_R1.fastq'),
        r2 = str(WD/'{sample}_R2.fastq')
    output:
        aq=temp('paired/{sample}.assembled.fastq'),
        dq=temp('paired/{sample}.discarded.fastq'),
        uf=temp('paired/{sample}.unassembled.forward.fastq'),
        ur=temp('paired/{sample}.unassembled.reverse.fastq'),
        fa='paired/{sample}.fasta'
    params:
        out_fp = 'paired/{sample}'
    log:
        expand('logs/pair/{{sample}}-{starttime}', starttime=starttime)
    shell:
        """
        pear -f {input.r1} -r {input.r2} -o {params.out_fp} >& {log} &&\
        fq2fa paired/{wildcards.sample}.assembled.fastq {output[fa]}
        """
        
rule _filter_len:
    input:
        '{filename}.fa'
    output:
        '{filename}.fa.{len}filtered'
    run:
        shell(
            """
            vsearch --sortbylength {input} \
            --relabel {wildcards.filename} --minseqlength {wildcards.len} \
            --output {output} >& /dev/null
            """)
        
rule _combine_individual_assemblers:
    input:
        A='abyss/{sample}/{sample}-unitigs.fa.150filtered',
        S='soapdenovo/{sample}/{sample}.contig.fa.150filtered',
        a='abyss_partition/{sample}/unitig-combined.fa.150filtered',
        I='idba_ud/{sample}/contig.fa.150filtered'
#        V='metavelvet/{sample}/meta-velvetg.contigs.fa.150filtered'
    output:
        'combined/{sample}-combined.fa'
    shell:
        """cat {input.A} {input.S} {input.a} {input.I} > {output}"""

rule _CAP3:
    input:
        'combined/{sample}-combined.fa'
    output:
        'combined/{sample}-combined-cap.fa'
    log:
        expand('logs/cap3/{{sample}}-{starttime}.log', starttime=starttime)
    shell:
        """
        {config[bin]}/cap3 {input} >& {log} && \
        cat {input}.cap.singlets {input}.cap.contigs > {output}
        """

rule _final_filter:
    input:
        'combined/{sample}-combined-cap.fa.300filtered'
    output:
        'final_contigs/{sample}-contigs.fa'
    shell:
        "cp {input} {output}"    
        
rule clean:
    shell:
        """
        rm -rf abyss abyss_partition chunks combined idba_ud metavelvet soapdenovo paired &&\
        echo "Cleanup finished." &&
        if [ -e final_contigs ]; then
            timestamp=$(date +%s) &&\
            mv final_contigs "final_contigs-$timestamp" &&\
            echo "Final contigs moved to final_contigs-$timestamp"
        fi
        """
rule clean_logs:
    shell:
        "rm -rf logs && echo 'Logs deleted'"
