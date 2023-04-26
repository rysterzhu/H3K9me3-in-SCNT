ddir=~/workspace/9.NT-ChIP/2.public/1.BWA/c.merge_rep
wdir=~/workspace/9.NT-ChIP/2.public/1.BWA/h.downsample-80M
mkdir -p $wdir/logs
SEED=1003; downSize=80000000 #downsample的目标数据量和种子
for i in $ddir/*.merge.bam; do (k=$(basename $i .merge.bam);
P=`awk -v d=$downSize '/in total/{if($1>d){p=d/$1}else{p=1.0};printf("%g",p)}' $ddir/logs/$k.merge.flag.txt`
nohup picard DownsampleSam I=$i O=$wdir/$k.downsample.bam R=$SEED P=$P A=0.00001 S=ConstantMemory > $wdir/logs/$k.log 
samtools sort -@ 8 $wdir/$k.downsample.bam -o $wdir/$k.sorted.bam
samtools index -@ 8 $wdir/$k.sorted.bam
samtools flagstat -@ 8 $wdir/$k.sorted.bam > $wdir/logs/$k.flag.txt
) &
done
wait
ddir=~/workspace/9.NT-ChIP/2.public/1.BWA/h.downsample-80M
wdir=~/workspace/9.NT-ChIP/2.public/6.callPeaks/3.K9me3/
for k in ${NFC[@]}; do 
odir=$wdir/f.MACS2-broad-default-merge-downsample-80M; mkdir -p $odir/logs
nohup macs2 callpeak --broad -t $ddir/${k}_K9me3.sorted.bam -c $ddir/${k}_input.sorted.bam --outdir $odir -f BAMPE -g mm -n $k --to-large > $odir/logs/$k.log &
done
