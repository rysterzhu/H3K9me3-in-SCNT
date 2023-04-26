#!/bin/sh

#*****************************************************************************************
#    Filename:  BWA.sh
#     Creator:  Ryster Zhu
# Create Time:  2020
# Description:  BWA
#     Version:  1
#*****************************************************************************************
set -e

help()
{
    cat << HELP
    This is a shell script used to BWA and samtools 

HELP
    exit 0
}

thread=8
p=4
index=mm10

while [ -n "$1" ]; do
case "$1" in
    -h) help; shift 1;;
    -i) in_file=$2; shift 2;;
	-w) wdir=$2; shift 2;;
    -d) ddir=$2; shift 2;;
    -p) p=$2; shift 2;;
    -t) thread=$2; shift 2;;
	-x) index=$2; shift 2;;
    *) echo "error: no such option $1. -h for help";exit 1;;
esac
done

#trap "exec 1000>&-;exec 1000<&-;exit 0" 2 
mkfifo temp.$$.fifo 
exec 1000<>temp.$$.fifo 
for((i=0;i<$p;i++));do echo >&1000; done

mkdir -p $wdir/logs $wdir/a.raw/logs $wdir/b.filtered/logs

for name in $ddir/*R1.fastq.gz; do 
read -u 1000 
{
	k=$(basename $name .R1.fastq.gz)
	bwa mem -t $thread $BWA_INDEXES/${index}.fa ${name} ${name/R1/R2} > $wdir/$k.sam 2>$wdir/logs/$k.log
	
	samtools flagstat -@ $thread $wdir/$k.sam > $wdir/logs/$k.flag.txt
	
	samtools sort -@ $thread $wdir/$k.sam -o $wdir/a.raw/$k.temp.bam
	samtools flagstat -@ $thread $wdir/a.raw/$k.temp.bam > $wdir/logs/$k.flag.txt.2
	
	#picard MarkDuplicates INPUT=$wdir/a.raw/$k.temp.bam OUTPUT=$wdir/a.raw/$k.sorted.bam METRICS_FILE=$wdir/a.raw/logs/$k.metrics > $wdir/a.raw/logs/$k.log 2>&1
	sambamba markdup -t $thread $wdir/a.raw/$k.temp.bam $wdir/a.raw/$k.sorted.bam > $wdir/a.raw/logs/$k.log 2>&1
	samtools flagstat -@ $thread $wdir/a.raw/$k.sorted.bam > $wdir/a.raw/logs/$k.txt
	
	samtools view -@ $thread -F 4 $wdir/a.raw/$k.sorted.bam -o $wdir/b.filtered/$k.sorted.bam
	samtools index -@ $thread $wdir/b.filtered/$k.sorted.bam
	samtools flagstat -@ $thread $wdir/b.filtered/$k.sorted.bam > $wdir/b.filtered/logs/$k.txt
	#rm $wdir/$k.sam $wdir/a.raw/$k.temp.bam 
	echo >&1000 
	echo $k" done"
}&   
done

wait
echo "All threads done."
rm -rf temp.$$.fifo 
