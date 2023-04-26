wdir=~/workspace/9.NT-ChIP/6.K9me3_analysis/5.YY1_map/4.log2ratio-befor2cell
pdir=~/workspace/9.NT-ChIP/3.deeptools/c.TF-blastocyst/mergedpeak
ddir=~/workspace/9.NT-ChIP/3.deeptools/5.compareByInput/7.log2-rpkm


key=all
computeMatrix scale-regions -S $ddir/CC*bw $ddir/6h*bw $ddir/14h*bw $ddir/NT-2cell*bw $ddir/O*bw $ddir/S*bw $ddir/PN*bw $ddir/NF-2cell*bw --smartLabels \
-R $pdir/YY1.bed -o $wdir/$key.gz --outFileNameMatrix $wdir/$key.tab \
-m 1000 -a 2000 -b 2000 -bs 50 --sortRegions descend --sortUsing mean --sortUsingSamples 1 -p 32 &
wait
k=5
plotHeatmap -m $wdir/$key.gz -out $wdir/$key.kmeans$k.png --dpi 360 --kmeans $k --outFileSortedRegions $wdir/$key.kmeans$k.tab \
--sortRegions descend --sortUsing mean --colorList "blue,white,red" --colorNumber 100000 --whatToShow "plot, heatmap and colorbar" --zMin -3 --zMax 3 \
--startLabel "" --endLabel "" -T "K9me3 to YY1" &

#
ddir=~/workspace/9.NT-ChIP/3.deeptools/5.compareByInput/a.keep-lowqual
wdir=~/workspace/9.NT-ChIP/6.K9me3_analysis/5.YY1_map/5.log2ratio-kdm4b
pdir=~/workspace/9.NT-ChIP/3.deeptools/c.TF-blastocyst/mergedpeak
mkcd $wdir
key=kdm4b
computeMatrix scale-regions -S $ddir/CC*bw $ddir/6h*bw $ddir/14h*bw $ddir/NT-2cell*bw $ddir/NT-14h_kdm4b*bw $ddir/NT-l2cell_kdm4b*bw $ddir/Oocyte*bw $ddir/Sperm*bw $ddir/PN*bw $ddir/NF-2cell*bw --smartLabels \
-R $pdir/YY1.bed -o $wdir/$key.gz --outFileNameMatrix $wdir/$key.tab \
-m 1000 -a 2000 -b 2000 -bs 50 --sortRegions descend --sortUsing mean --sortUsingSamples 1 -p 32 &
wait
k=5
plotHeatmap -m $wdir/$key.gz -out $wdir/$key.kmeans$k.pdf --dpi 360 --kmeans $k --outFileSortedRegions $wdir/$key.kmeans$k.tab \
--sortRegions descend --sortUsing mean --colorList "blue,white,red" --colorNumber 100000 --whatToShow "plot, heatmap and colorbar" --zMin -3 --zMax 3 \
--startLabel "" --endLabel "" -T "K9me3 to YY1" &

