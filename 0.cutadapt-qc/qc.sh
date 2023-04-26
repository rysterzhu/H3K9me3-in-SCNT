ddir=~/data/9.NT-ChIP/revision.20230314
wdir=~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links/zd.20230314
mkdir -p $wdir ${wdir/0.links/1.qc} ${wdir/0.links/3.cutadapt} ${wdir/0.links/4.qc_cutadapt}

for i in 1 2; do 
ln -fs $ddir/_*R$i*.fastq.gz $wdir/.....R$i.fastq.gz
done 

nohup fastqc -o ${wdir/0.links/1.qc} $wdir/*gz -t 64 &
for i in $wdir/*R1.fastq.gz; do o=${i/0.links/3.cutadapt}
nohup cutadapt -j 8 --trim-n -q 25,25 -m 30 -a ATCAACGCAGAGTAC -a AGATCGGAAGAGC -A ATCAACGCAGAGTAC -A AGATCGGAAGAGC -o $o -p ${o/R1/R2} $i ${i/R1/R2} > ${o/R1.fastq.gz/log} &
done
wait
nohup fastqc -o ${wdir/0.links/4.qc_cutadapt} ${wdir/0.links/3.cutadapt}/*.gz -t 32 &
wait
cd ${wdir/0.links/1.qc}
multiqc --no-data-dir -f -dd 10 ./ -n $wdir/raw.qc.html
cd ${wdir/0.links/3.cutadapt}
multiqc --no-data-dir -f -dd 10 ./ -n $wdir/cutadapt.html
cd ${wdir/0.links/4.qc_cutadapt}
multiqc --no-data-dir -f -dd 10 ./ -n $wdir/clean.qc.html
wait-m
