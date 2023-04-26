#!/bin/sh

#*****************************************************************************************
#    Filename:  thread_cutrun.sh
#     Creator:  Ryster Zhu
# Create Time:  2020/01/26
# Description:  thread cut run bwa
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

while [ -n "$1" ]; do
case "$1" in
    -h) help; shift 1;;
	-w) wdir=$2; shift 2;;
    -d) ddir=$2; shift 2;;
    -p) p=$2; shift 2;;
    -t) thread=$2; shift 2;;
    -*) echo "error: no such option $1. -h for help";exit 1;;
     *) break;
esac
done

mkdir -p $wdir
mkfifo $wdir/temp.$$.fifo
exec 1000<>$wdir/temp.$$.fifo
for((i=0;i<$p;i++));do echo >&1000; done

mkdir -p $wdir/logs $wdir/a.raw/logs $wdir/b.filtered/logs

for name in $ddir/*R1.fastq.gz; do 
read -u 1000 
{
    k=$(basename $name .R1.fastq.gz)
	index=mm10
    bwa mem -t $thread $BWA_INDEXES/${index}.fa ${name} ${name/R1/R2} > $wdir/$k.sam 2>$wdir/logs/$k.log
    samtools flagstat -@ $thread $wdir/$k.sam > $wdir/logs/$k.flag.txt
    samtools sort -@ $thread $wdir/$k.sam -o $wdir/a.raw/$k.temp.bam
    sambamba markdup -t $thread $wdir/a.raw/$k.temp.bam $wdir/a.raw/$k.sorted.bam
    samtools flagstat -@ $thread $wdir/a.raw/$k.sorted.bam > $wdir/a.raw/logs/$k.flag.txt
    samtools view -@ $thread -F 1036 -q 1 $wdir/a.raw/$k.sorted.bam -o $wdir/b.filtered/$k.sorted.bam
    samtools index -@ $thread $wdir/b.filtered/$k.sorted.bam
    samtools flagstat -@ $thread $wdir/b.filtered/$k.sorted.bam > $wdir/b.filtered/logs/$k.flag.txt
	
	index=sacCer3
    bwa mem -t $thread $BWA_INDEXES/${index}.fa ${name} ${name/R1/R2} > $wdir/s$k.sam 2>$wdir/logs/s$k.log
    samtools flagstat -@ $thread $wdir/s$k.sam > $wdir/logs/s$k.flag.txt
    samtools sort -@ $thread $wdir/s$k.sam -o $wdir/a.raw/s$k.temp.bam
    sambamba markdup -t $thread $wdir/a.raw/s$k.temp.bam $wdir/a.raw/s$k.sorted.bam
    samtools flagstat -@ $thread $wdir/a.raw/s$k.sorted.bam > $wdir/a.raw/logs/s$k.flag.txt
    samtools view -@ $thread -F 1036 -q 1 $wdir/a.raw/s$k.sorted.bam -o $wdir/b.filtered/s$k.sorted.bam
    samtools index -@ $thread $wdir/b.filtered/s$k.sorted.bam
    samtools flagstat -@ $thread $wdir/b.filtered/s$k.sorted.bam > $wdir/b.filtered/logs/s$k.flag.txt
	echo >&1000 
	echo $k" done"
}&   
done

wait
echo "All threads done."
rm -rf $wdir/temp.$$.fifo 
rm $wdir/*sam $wdir/*/*temp.bam




