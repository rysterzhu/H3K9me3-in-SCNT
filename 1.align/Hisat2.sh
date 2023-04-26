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

