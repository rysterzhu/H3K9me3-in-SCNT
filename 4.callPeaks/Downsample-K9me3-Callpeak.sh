################0
#downsample redun K9me3和input到80M test
ddir=~/workspace/9.NT-ChIP/1.align/[01]*/c.merge_rep
wdir=~/workspace/9.NT-ChIP/1.align/1.K9me3_BWA/h.downsample-80M
mkdir -p $wdir/logs
SEED=1005; downSize=80000000 #downsample的目标数据量和种子
for i in $ddir/*.redun.bam; do (k=$(basename $i .redun.bam);
P=`awk -v d=$downSize '/in total/{if($1>d){p=d/$1}else{p=1.0};printf("%g",p)}' $ddir/logs/$k.redun.flag.txt`
nohup picard DownsampleSam I=$i O=$wdir/$k.downsample.bam R=$SEED P=$P A=0.00001 S=ConstantMemory > $wdir/logs/$k.log 
samtools sort -@ 8 $wdir/$k.downsample.bam -o $wdir/$k.sorted.bam
samtools index -@ 8 $wdir/$k.sorted.bam
samtools flagstat -@ 8 $wdir/$k.sorted.bam > $wdir/logs/$k.flag.txt
rm $wdir/$k.downsample.bam
) &
done
wait
ddir=~/workspace/9.NT-ChIP/1.align/1.K9me3_BWA/h.downsample-80M
wdir=~/workspace/9.NT-ChIP/4.callPeaks/4.K9me3/
for k in ${NTC[@]}; do 
odir=$wdir/h.MACS2-broad-default-redun-downsample-80M; mkdir -p $odir/logs
nohup macs2 callpeak --broad -t $ddir/${k}_K9me3.sorted.bam -c $ddir/${k}_input.sorted.bam --outdir $odir -f BAMPE -g mm -n $k --keep-dup all --to-large > $odir/logs/$k.log &
done
wait-m 


#downsample merge K9me3和input到80M test
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
wait-m 80M

ddir=~/workspace/9.NT-ChIP/2.public/5.deeptools/4.bamCoverage/K9me3
wdir=~/workspace/9.NT-ChIP/2.public/6.callPeaks/3.K9me3/f.MACS2-broad-default-merge-downsample-80M
for k in ${NFC[@]}; do 
multiBigwigSummary BED-file -b $ddir/${k}_K9me3.bw --BED $wdir/${k}_peaks.broadPeak -o $wdir/${k}.npz --outRawCounts $wdir/${k}.rpkm.bg -p 8 &
done


########Kdm4b reps 验证
ddir=~/workspace/9.NT-ChIP/1.align/g.K9me3_kdm4b/b.filtered
wdir=~/workspace/9.NT-ChIP/4.callPeaks/1.MACS2/8.kdm4b-reps-test
for i in $ddir/*K9me3*bam; do mkdir -p $wdir/logs;k=${i##*/};
nohup macs2 callpeak --broad -t $i -c $ddir/${k%%_K9me3_rep*}_input_rep1.sorted.bam --outdir $wdir -f BAMPE -g mm -n $k --keep-dup all --to-large > $wdir/logs/$k.log &
done
wait-m 


for i in *_peaks.broadPeak; do k=$(basename $i _peaks.broadPeak); nohup annotatePeaks.pl $i mm10 -annStats annotate/$k.txt > annotate/$k.log & done


