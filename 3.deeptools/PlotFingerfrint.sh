ddir=~/workspace/9.NT-ChIP/1.align/0.input_usable
wdir=~/workspace/9.NT-ChIP/3.deeptools/4.GCBias/0.input_usable
mkdir -p $wdir
for i in $ddir/*.sorted.bam; do k=$wdir/$(basename $i .sorted.bam)
computeGCBias -p 4 -b $i --effectiveGenomeSize 2150570000 -g ~/ann/mm10.2bit -o $k.txt --biasPlot $k.png &
done 

wdir=~/workspace/9.NT-ChIP/3.deeptools/3.plotFingerprint/0.input_usable
ddir=~/workspace/9.NT-ChIP/1.align/0.input_usable
keys=(cc 6h 14h 2cell 4cell 8cell morula icm te blast)
mkdir -p $wdir
for k in ${keys[@]}; do 
plotFingerprint -b $ddir/$k*sorted.bam -o $wdir/$k.png --ignoreDuplicates -T "" -p 4 --samFlagInclude 2 -bs 5000 --numberOfSamples 500000 &
done



