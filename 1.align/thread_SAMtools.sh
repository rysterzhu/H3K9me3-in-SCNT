#trap "exec 1000>&-;exec 1000<&-;exit 0" 2 
mkfifo temp.fifo 
exec 1000<>temp.fifo 
for((i=0;i<4;i++));do echo >&1000; done

thread=4
wdir=~/workspace/9.NT-ChIP/1.align/5.BWA
mkdir -p $wdir/flagstats $wdir/logs
for i in $wdir/*sam; do 
read -u 1000 
{
	k=$(basename $i .sam)
	samtools view -@ $thread -Shb $i -o ${i/sam/bam}
	samtools flagstat -@ $thread ${i/sam/bam} > $wdir/logs/$k.flag.txt 
	#samtools view -@ $thread -F 4 -q 30 -o ${i/sam/bam} ${i/sam/filtered.bam}
	#samtools sort -@ $thread ${i/sam/filtered.bam} -o ${i/sam/sorted.bam}
	#samtools index -@ $thread ${i/sam/sorted.bam}
	#samtools flagstat -@ $thread ${i/sam/bam} > $wdir/flagstats/$k.txt
	#rm $i 
	echo >&1000 
	echo $out" done"
}&   
done

wait
echo "All threads done."
rm -rf temp.fifo 


