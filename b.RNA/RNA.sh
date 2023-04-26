nohup featureCounts -T 64 -p -B -C -M -t exon -g gene_id -a /home/share/gtf/mm10.gtf -o $wdir/TE.featureCount $ddir1/te*sorted.bam $ddir2/RWTE*sorted.bam >> $wdir/featureCount.log &

nohup featureCounts -T 64 -p -B -C -M -t exon -g gene_id -a /home/share/gtf/mm10.gtf -o $wdir/ICM2.featureCount $ddir1/icm*sorted.bam $ddir2/RWICM*sorted.bam >> $wdir/featureCount.log &

