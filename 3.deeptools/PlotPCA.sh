
wdir=~/workspace/9.NT-ChIP/3.deeptools/6.plotPCA/6.K9me3_bwa_merged
ddir1=~/workspace/9.NT-ChIP/4.callPeaks/4.K9me3/5.MACS-nomodel-bwa
ddir2=~/workspace/9.NT-ChIP/2.public/6.callPeaks/3.K9me3/4.MACS-nomodel-bwa-rawdata-downsample

cat $ddir1/*_peaks.bed $ddir2/*_peaks.bed | sort -k1,1 -k2n,2 > test
mergeBed -i test -d 1000 > merge.1k.bed #246525ge 

wdir=~/workspace/9.NT-ChIP/3.deeptools/6.plotPCA/6.K9me3_bwa_merged
ddir1=~/workspace/9.NT-ChIP/3.deeptools/5.compareByInput/4.log2_bwa
ddir2=~/workspace/9.NT-ChIP/2.public/5.deeptools/3.bamCompare/3.log2_bwa
multiBigwigSummary BED-file --BED $wdir/merge.1k.bed \
-b `for i in ${NFC[@]}; do echo -n $ddir2/${i}_K9me3.bw" ";done; for i in ${NTC[@]}; do echo -n $ddir1/${i}_K9me3.bw" ";done`\
-o $wdir/all.npz --outRawCounts $wdir/all.tab -p 32 \
-l `for i in ${NFC[@]}; do echo -n NF-${i}" ";done; for i in ${NTC[@]}; do echo -n NT-${i}" ";done` &

plotPCA -in $wdir/all.npz -o $wdir/1kb-transpose.png --transpose --outFileNameData $wdir/1kb-transpose.tab --plotWidth 10 &
plotPCA -in $wdir/all.npz -o $wdir/1kb-rowCenter.png --rowCenter --outFileNameData $wdir/1kb-rowCenter.tab --plotWidth 10 &



multiBigwigSummary bins -bs 20000 \
-b `for i in ${NFC[@]}; do echo -n $ddir2/${i}_K9me3.bw" ";done; for i in ${NTC[@]}; do echo -n $ddir1/${i}_K9me3.bw" ";done`\
-o $wdir/bins.npz --outRawCounts $wdir/bins.tab -p 32 \
-l `for i in ${NFC[@]}; do echo -n NF-${i}" ";done; for i in ${NTC[@]}; do echo -n NT-${i}" ";done` &
plotPCA -in $wdir/bins.npz -o $wdir/bins-transpose.png --transpose --outFileNameData $wdir/bins-transpose.tab --plotWidth 10 &
plotPCA -in $wdir/bins.npz -o $wdir/bins-rowCenter.png --rowCenter --outFileNameData $wdir/bins-rowCenter.tab --plotWidth 10 &
