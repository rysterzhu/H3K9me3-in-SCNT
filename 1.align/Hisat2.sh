ddir=~/workspace/9.NT-ChIP/0.cutadapt-qc/3.cutadapt/o.20191021
wdir=~/workspace/9.NT-ChIP/1.align/9.RNA_hisat2
mkdir -p $wdir/logs
for i in $ddir/*R1.fastq.gz; do o=`basename $i .R1.fastq.gz`
nohup hisat2 --dta-cufflinks --no-discordant --no-mixed -p 4 --no-unal -x /home/share/hisat2_index/new_index/mm10_tran -1 $i -2 ${i/R1./R2.} -S $wdir/$o.sam > $wdir/logs/$o.hisat2.log 2>&1 &
done
wait-m hisat2
for i in $wdir/*sam; do
(samtools view -Shb $i -o ${i/sam/bam}
samtools sort ${i/sam/bam} -o ${i/sam/sorted.bam}
samtools index ${i/sam/sorted.bam}
) &
done

###########BLASTn
ddir=~/workspace/9.NT-ChIP/0.cutadapt-qc/3.cutadapt/o.20191021
wdir=~/workspace/9.NT-ChIP/1.align/9.RNA_hisat2/blastn.test
step=100
for i in $ddir/*fastq.gz;do o=$(basename $i .fastq.gz)
(zcat  $i| awk -v step=$step 'NR%(4*step)==1{print ">"$0} NR%(4*step)==2{print $0}' > $wdir/$o.pick.fa
nohup blastn -max_target_seqs 1 -num_threads 16 -query $wdir/$o.pick.fa -db /home/share/BlastDB/blastdb_new/nt_new/nt \
-outfmt "7 qacc sacc evalue pident gaps qcovus stitle ssciname" \
-out $wdir/$o.self -evalue 1e-30 -perc_identity 99.33
echo -e 'count\tratio\tspecies' > $o.statistic.tab
awk 'BEGIN {FS="[\t]+" ; OFS="\t"} $0!~/^#/{organism[$8]+=1 ; total+=1} END {for (item in organism) print organism[item],organism[item]/total*100,item }' $o.self | sort -k2nr,2 >> $o.statistic.tab  )&
done
#不是污染问题
wdir=~/workspace/9.NT-ChIP/1.align/b.links-RNA-reps
mkdir -p $wdir/b.picard
for i in $wdir/*sorted.bam;do o=$(basename $i .sorted.bam)
nohup picard CollectRnaSeqMetrics REF_FLAT=/home/share/ann/mm10.refFlat STRAND_SPECIFICITY=NONE INPUT=$i OUTPUT=$wdir/b.picard/${o}.collect_rna_metrics.txt CHART_OUTPUT=$wdir/b.picard/${o}.collect_rna_metrics_chart.pdf > $wdir/b.picard/${o}.collect_rna_metrics.log
done
#效果不好

#featurecount
ddir=~/workspace/9.NT-ChIP/1.align/9.RNA_hisat2
wdir=~/workspace/9.NT-ChIP/b.RNA/1.DESeq2
mkd $wdir/logs
nohup featureCounts -T 36 -p -B -C -M -t exon -g gene_id -a /home/share/gtf/mm10.gtf -o $wdir/NT.featureCount $ddir/*sorted.bam > $wdir/logs/NT.featureCount.log &

#cufflinks
ddir=~/workspace/9.NT-ChIP/1.align/9.RNA_hisat2
wdir=~/workspace/9.NT-ChIP/b.RNA/2.cufflinks
mkd $wdir/logs
for i in $ddir/*sorted.bam;do o=$(basename $i .sorted.bam)
nohup cufflinks2 -o $wdir/$o -u -p 4 -G /home/share/gtf/mm10.gtf -b /home/share/bowtie_index/mm10.fa $i > $wdir/logs/$o.log &
done

#############################新数据
ddir=~/workspace/9.NT-ChIP/0.cutadapt-qc/3.cutadapt/q.20200614
wdir=~/workspace/9.NT-ChIP/1.align/9.RNA_hisat2/new_batch
mkdir -p $wdir/logs
for i in $ddir/*RNA*R1.fastq.gz; do o=`basename $i .R1.fastq.gz`
nohup hisat2 --dta-cufflinks --no-discordant --no-mixed -p 4 --no-unal -x /home/share/hisat2_index/new_index/mm10_tran -1 $i -2 ${i/R1./R2.} -S $wdir/$o.sam > $wdir/logs/$o.hisat2.log 2>&1 &
done
wait-m hisat2
for i in $wdir/*sam; do
(samtools view -Shb $i -o ${i/sam/bam}
samtools sort ${i/sam/bam} -o ${i/sam/sorted.bam}
samtools index ${i/sam/sorted.bam}
) &
done

~/workspace/9.NT-ChIP/1.align/b.links-RNA-reps
ln -s ~/workspace/9.NT-ChIP/1.align/9.RNA_hisat2/*sorted.bam* ./
ln -s ~/workspace/9.NT-ChIP/1.align/9.RNA_hisat2/new_batch/*sorted.bam* ./

4cell_RNA_rep5量太少
nohup featureCounts -T 36 -p -B -C -M -t exon -g gene_id -a /home/share/gtf/mm10.gtf -o all.featureCount *sorted.bam > all.log &

wdir=~/workspace/9.NT-ChIP/1.align/9.RNA_hisat2/a.NT
mkdir -p $wdir/1.picard
mult 10
for i in $wdir/*sorted.bam;do read -u 9;{ o=$(basename $i .sorted.bam)
nohup picard CollectRnaSeqMetrics REF_FLAT=/home/share/ann/mm10.refFlat STRAND_SPECIFICITY=NONE INPUT=$wdir/$o.sorted.bam OUTPUT=$wdir/1.picard/${o}.collect_rna_metrics.txt CHART_OUTPUT=$wdir/1.picard/${o}.collect_rna_metrics_chart.pdf > $wdir/1.picard/${o}.collect_rna_metrics.log 
echo >&9; }&
done

###########################2020年8月17日 IL
ddir=~/workspace/9.NT-ChIP/0.cutadapt-qc/3.cutadapt/r.20200817
wdir=~/workspace/9.NT-ChIP/1.align/9.RNA_hisat2/IL
mkdir -p $wdir/logs
for i in $ddir/*RNA*R1.fastq.gz; do o=`basename $i .R1.fastq.gz`
nohup hisat2 --dta-cufflinks --no-discordant --no-mixed -p 4 --no-unal -x /home/share/hisat2_index/new_index/mm10_tran -1 $i -2 ${i/R1./R2.} -S $wdir/$o.sam > $wdir/logs/$o.hisat2.log 2>&1 &
done

for i in $wdir/*sam; do
(samtools view -Shb $i -o ${i/sam/bam}
samtools sort ${i/sam/bam} -o ${i/sam/sorted.bam}
samtools index ${i/sam/sorted.bam}
) &
done

ln -s ~/workspace/9.NT-ChIP/1.align/9.RNA_hisat2/IL/*sorted.bam* ~/workspace/9.NT-ChIP/1.align/b.links-RNA-reps/
key=20200817
nohup featureCounts -T 36 -p -B -C -M -t exon -g gene_id -a /home/share/gtf/mm10.gtf -o $key.featureCount *sorted.bam > $key.log &


###########################2020年9月20日 IL
ddir=~/workspace/9.NT-ChIP/0.cutadapt-qc/3.cutadapt/s.20200920
wdir=~/workspace/9.NT-ChIP/1.align/9.RNA_hisat2/IL
mkdir -p $wdir/logs
for i in $ddir/*RNA*R1.fastq.gz; do o=`basename $i .R1.fastq.gz`
nohup hisat2 --dta-cufflinks --no-discordant --no-mixed -p 4 --no-unal -x /home/share/hisat2_index/new_index/mm10_tran -1 $i -2 ${i/R1./R2.} -S $wdir/$o.sam > $wdir/logs/$o.hisat2.log 2>&1 &
done

for i in $wdir/*[456].sam; do
(samtools view -Shb $i -o ${i/sam/bam}
samtools sort ${i/sam/bam} -o ${i/sam/sorted.bam}
samtools index ${i/sam/sorted.bam}
) &
done

ln -s ~/workspace/9.NT-ChIP/1.align/9.RNA_hisat2/IL/*sorted.bam* ~/workspace/9.NT-ChIP/1.align/b.links-RNA-reps/
key=20200920
nohup featureCounts -T 36 -p -B -C -M -t exon -g gene_id -a /home/share/gtf/mm10.gtf -o $key.featureCount *sorted.bam > $key.log &

wdir=~/workspace/9.NT-ChIP/1.align/9.RNA_hisat2/IL
mkdir -p $wdir/1.picard
for i in $wdir/*sorted.bam;do o=$(basename $i .sorted.bam)
nohup picard CollectRnaSeqMetrics REF_FLAT=/home/share/ann/hg19.refFlat STRAND_SPECIFICITY=NONE INPUT=$i OUTPUT=$wdir/1.picard/${o}.collect_rna_metrics.txt CHART_OUTPUT=$wdir/1.picard/${o}.collect_rna_metrics_chart.pdf > $wdir/1.picard/${o}.collect_rna_metrics.log &
done



#############################siKat5
#v.
ddir=~/workspace/9.NT-ChIP/0.cutadapt-qc/3.cutadapt/w.20201227
wdir=~/workspace/9.NT-ChIP/1.align/9.RNA_hisat2/siKat5
for i in $ddir/*RNA*R1.fastq.gz; do o=`basename $i .R1.fastq.gz`;
nohup hisat2 --dta-cufflinks --no-discordant --no-mixed -p 4 --no-unal -x /home/share/hisat2_index/new_index/mm10_tran -1 $i -2 ${i/R1./R2.} -S $wdir/$o.sam > $wdir/logs/$o.hisat2.log 2>&1 &
done
wait
for i in $wdir/*.sam; do
(samtools view -Shb $i -o ${i/sam/bam}
samtools sort ${i/sam/bam} -o ${i/sam/sorted.bam}
samtools index ${i/sam/sorted.bam}
) &
done
##########删除了siKat5,重做所有的si在Morula ICM TE的smart-seq
ddir=~/workspace/9.NT-ChIP/0.cutadapt-qc/3.cutadapt/6.si-MIT-smart
wdir=~/workspace/9.NT-ChIP/1.align/9.RNA_hisat2/d.siMIT
mkd $wdir/logs
for i in $ddir/*RNA*R1.fastq.gz; do o=`basename $i .R1.fastq.gz`;
nohup hisat2 --dta-cufflinks --no-discordant --no-mixed -p 4 --no-unal -x /home/share/hisat2_index/new_index/mm10_tran -1 $i -2 ${i/R1./R2.} -S $wdir/$o.sam > $wdir/logs/$o.hisat2.log 2>&1 &
done
wait-m hisat2
for i in $wdir/*.sam; do
(samtools view -Shb $i -o ${i/sam/bam}
samtools sort ${i/sam/bam} -o ${i/sam/sorted.bam}
samtools index ${i/sam/sorted.bam}
samtools flagstat ${i/sam/sorted.bam} > ${i/sam/flagstat}
) &
done
wait-m samtools
 
wdir=~/workspace/9.NT-ChIP/1.align/9.RNA_hisat2/d.siMIT
mkdir -p $wdir/1.picard
for i in $wdir/*sorted.bam;do o=$(basename $i .sorted.bam)
nohup picard CollectRnaSeqMetrics REF_FLAT=/home/share/ann/mm10.refFlat STRAND_SPECIFICITY=NONE INPUT=$i OUTPUT=$wdir/1.picard/${o}.collect_rna_metrics.txt CHART_OUTPUT=$wdir/1.picard/${o}.collect_rna_metrics_chart.pdf > $wdir/1.picard/${o}.collect_rna_metrics.log &
done





########################siMIT NT 
ddir=~/workspace/9.NT-ChIP/0.cutadapt-qc/3.cutadapt/8.siMIT-NT/
ddir=~/workspace/9.NT-ChIP/0.cutadapt-qc/3.cutadapt/y.20210425/2.NT-RNA
wdir=~/workspace/9.NT-ChIP/1.align/9.RNA_hisat2/e.siMIT-NT
mkd $wdir/logs
for i in $ddir/*RNA*R1.fastq.gz; do o=`basename $i .R1.fastq.gz`;
nohup hisat2 --dta-cufflinks --no-discordant --no-mixed -p 4 --no-unal -x /home/share/hisat2_index/new_index/mm10_tran -1 $i -2 ${i/R1./R2.} -S $wdir/$o.sam > $wdir/logs/$o.hisat2.log 2>&1 &
done
wait-m hisat2
for i in $wdir/*.sam; do
(samtools view -Shb $i -o ${i/sam/bam}
samtools sort ${i/sam/bam} -o ${i/sam/sorted.bam}
samtools index ${i/sam/sorted.bam}
samtools flagstat ${i/sam/sorted.bam} > ${i/sam/flagstat}
) &
done
wait-m samtools

mkdir -p $wdir/1.picard
for i in $wdir/*sam;do o=$(basename $i .sam)
nohup picard CollectRnaSeqMetrics REF_FLAT=/home/share/ann/mm10.refFlat STRAND_SPECIFICITY=NONE INPUT=$wdir/$o.sorted.bam OUTPUT=$wdir/1.picard/${o}.collect_rna_metrics.txt CHART_OUTPUT=$wdir/1.picard/${o}.collect_rna_metrics_chart.pdf > $wdir/1.picard/${o}.collect_rna_metrics.log &
done


ddir=~/workspace/9.NT-ChIP/0.cutadapt-qc/3.cutadapt/9.siMIT-mph
ddir=~/workspace/c.chrodoma/0.prepare/2.RNA.CH22.HNP/
wdir=~/workspace/9.NT-ChIP/1.align/9.RNA_hisat2/f.siMIT-mph
mkd $wdir/logs
for i in $ddir/*mph*R1.fastq.gz; do o=`basename $i .R1.fastq.gz`;
nohup hisat2 --dta-cufflinks --no-discordant --no-mixed -p 4 --no-unal -x /home/share/hisat2_index/new_index/mm10_tran -1 $i -2 ${i/R1./R2.} -S $wdir/$o.sam > $wdir/logs/$o.hisat2.log 2>&1 &
done
wait-m hisat2
for i in $wdir/*.sam; do
(samtools view -Shb $i -o ${i/sam/bam}
samtools sort ${i/sam/bam} -o ${i/sam/sorted.bam}
samtools index ${i/sam/sorted.bam}
samtools flagstat ${i/sam/sorted.bam} > ${i/sam/flagstat}
) &
done
wait-m samtools


