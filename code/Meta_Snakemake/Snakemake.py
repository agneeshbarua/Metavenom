#Wildcards-----------------------------------------------
SAMPLES = ["SRR5351966","SRR5351967","SRR5351968","SRR5351969","SRR5351970","SRR5351971"]

#Rules----------------------i------------------------------

rule all:
        input: expand("results/output/{samples}",samples = SAMPLES)

rule kallisto_quant:
        input: r1 = "../Salivary-gland/{samples}_1.fastq",
               r2 = "../Salivary-gland/{samples}_2.fastq",
               idx ='../Salivary-gland/kallisto/mus_musculus/transcriptome.idx'
        output:
                "results/output/{samples}"
        threads: 8
        shell: "kallisto quant -i {input.idx} -o {output} -b 100 {input.r1} {input.r2}"