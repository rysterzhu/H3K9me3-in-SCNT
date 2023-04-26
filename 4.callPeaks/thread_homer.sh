#!/bin/sh

#*****************************************************************************************
#    Filename:  thread_homer.sh
#     Creator:  Ryster Zhu
# Create Time:  2021/03/9
#     Version:  1
#*****************************************************************************************
set -e

thread=32
p=8
wdir=~/workspace/9.NT-ChIP/4.callPeaks/c.homer-reps
ddir1=~/workspace/9.NT-ChIP/1.align/f.K9me3_keep_lowqual/1.merge-reps
ddir2=~/workspace/9.NT-ChIP/1.align/f.K9me3_keep_lowqual/0.links-raw-reps
mkdir -p $wdir/logs

mkfifo $wdir/temp.$$.fifo
exec 1000<>$wdir/temp.$$.fifo
for((i=0;i<$p;i++));do echo >&1000; done



for i in $ddir1/[C6]*_input.redun.bam; do
read -u 1000 
{
k=$(basename $i );k=${k%%.*}
samtools view -@ $thread -F 1036 -h $i -o $wdir/$k.sam
nohup makeTagDirectory $wdir/$k.tag/ $wdir/$k.sam -keepOne -format sam -genome mm10 > $wdir/logs/$k.make.log
rm -rf $wdir/$k.sam
echo >&1000 
echo $k" done"
}&  
done


wait
echo "All threads done."
rm -rf $wdir/temp.$$.fifo 



