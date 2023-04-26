trap "exec 1000>&-;exec 1000<&-;exit 0" 2 
mkfifo temp.fifo 
exec 1000<>temp.fifo 
for((i=0;i<8;i++));do echo >&1000; done
p=8

ddir=~/workspace/8.NT-HiC/b.RNA/0.data/2.Gaorui_RNA/2.cutadapt
wdir=~/workspace/8.NT-HiC/b.RNA/1.align/2.Gaorui_RNA
mkdir -p $wdir/logs $wdir/4.bg_50k
for i in $wdir/1.sortedBam/*.sorted.bam; do o=${i##*/};o=${o%%.*}
read -u 1000 
{
	nohup hisat2 --dta-cufflinks --no-discordant --no-mixed -p $p --no-unal -x /home/share/hisat2_index/new_index/mm10_tran -1 $i -2 ${i/R1./R2.} -S $wdir/$o.sam > $wdir/logs/$o.hisat2.log 2>&1
	samtools view -Sbh -@ $p $wdir/$o.sam -o $wdir/$o.bam
	samtools sort -@ $p $wdir/$o.bam -o $wdir/$o.sorted.bam
	samtools index $wdir/$o.sorted.bam
	rm $wdir/$o.sam $wdir/$o.bam
	nohup bamCoverage -p $p -b $wdir/1.sortedBam/$o.sorted.bam -of bedgraph -o $wdir/4.bg_50k/$o.bw --normalizeUsing CPM --binSize 50000 --smoothLength 200000 > $wdir/logs/$o.bw.log
	echo "nohup bamCoverage -p 5 -b $wdir/1.sortedBam/$o.sorted.bam -of bedgraph -o $wdir/4.bg_50k/$o.bw --normalizeUsing CPM --binSize 50000 --smoothLength 200000 > $wdir/logs/$o.bw.log" >> $wdir/logs/run.log 
	echo >&1000 
	echo $i" done"
}&   
done &

wait
echo "All threads done."
rm -rf temp.fifo