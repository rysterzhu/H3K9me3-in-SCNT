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
mkfifo temp.$$.fifo
exec 1000<>$wdir/temp.$$.fifo
for((i=0;i<$p;i++));do echo >&1000; done


for k in `for i in $ddir/*.sorted.bam; do o=$(basename $i );echo ${o%%_rep*};done | sort | uniq`;do 
read -u 1000 
{   samtools merge -@ 8 -f $wdir/${k}.redun.bam $ddir/${k}*sorted.bam;
	samtools flagstat -@ 8 $wdir/${k}.redun.bam > $wdir/logs/${k}.redun.flag.txt; 
	samtools index -@ 8 $wdir/${k}.redun.bam
	echo >&1000 
	echo $k" done"
}&   
done

wait
echo "All threads done."
rm -rf temp.$$.fifo 

 


