
Get the high coverage bams
    wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA172/ERA172924/bam/NA12889_S1.bam
    wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA172/ERA172924/bam/NA12890_S1.bam
    wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA172/ERA172924/bam/NA12878_S1.bam

    samtools index NA12889_S1.bam
    samtools index NA12890_S1.bam
    samtools index NA12878_S1.bam

    ~/src/lumpy-sv/scripts/lumpy_smooth \
        NA12878_S1.bam NA12889_S1.bam NA12890_S1.bam \
    > NA12878.NA12889.NA12890.vcf

    svtyper \
        -i <(bcftools view -G ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz 2:184569079  | awk '{ if ($1 ~ /^#/) {print $0;} else {print "chr"$0;} }')
        -B NA12878_S1.bam,NA12890_S1.bam,NA12889_S1.bam \
    | bgzip -c > NA12878.NA12889.NA12890.svt.vcf.gz

    bcftools view \
        -i 'SVTYPE="INV" && QUAL>0' \
        NA12878.NA12889.NA12890.svt.vcf.gz \
        chr2:89161083-89185670 \
        chr12:47290448-47309758 \
        chr12:12544868-12546613 \
    | bcftools annotate -x INFO/PRPOS,INFO/PREND \
    | bgzip -c \
    > invs.vcf.gz

    
Get the low coverage bams
