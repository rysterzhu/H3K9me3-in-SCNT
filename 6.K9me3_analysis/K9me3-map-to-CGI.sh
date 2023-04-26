wdir=~/workspace/9.NT-ChIP/6.K9me3_analysis/1.CpG_analysis_20180608/4.log2ratio-before2cell
pdir=~/ann/mm10.CGI.bed
ddir=~/workspace/9.NT-ChIP/3.deeptools/5.compareByInput/7.log2-rpkm

mkd $wdir

key=before14h
computeMatrix scale-regions -S $ddir/CC*bw $ddir/6h*bw $ddir/14h*bw $ddir/O*bw $ddir/S*bw $ddir/PN*bw --smartLabels \
-R $pdir -o $wdir/$key.gz --outFileNameMatrix $wdir/$key.tab \
-m 1000 -a 2000 -b 2000 -bs 50 --sortRegions descend --sortUsing mean --sortUsingSamples 1 -p 32 &
wait
plotHeatmap -m $wdir/$key.gz -out $wdir/$key.order-6h.png --dpi 360 --outFileSortedRegions $wdir/$key.order-6h.tab \
--sortRegions descend --sortUsing mean --sortUsingSamples 2 --colorList "blue,white,red" --colorNumber 100000 --whatToShow "plot, heatmap and colorbar" --zMin -3 --zMax 3 \
--startLabel "" --endLabel "" -T "K9me3 to CGI" &
plotHeatmap -m $wdir/$key.gz -out $wdir/$key.order-6h.pdf --dpi 360 --outFileSortedRegions $wdir/$key.order-6h.tab \
--sortRegions descend --sortUsing mean --sortUsingSamples 2 --colorList "blue,white,red" --colorNumber 100000 --whatToShow "plot, heatmap and colorbar" --zMin -3 --zMax 3 \
--startLabel "" --endLabel "" -T "K9me3 to CGI" &

###############sih2 
###############################
wdir=~/workspace/9.NT-ChIP/6.K9me3_analysis/1.CpG_analysis_20180608/6.sih2-bamCoverage
ddir=~/workspace/9.NT-ChIP/d.cut_run/2.bamCoverage
mkcd $wdir
key=before14h-rep34
computeMatrix scale-regions -S $ddir/6hpa*[rep34].bw --smartLabels \
-R ~/ann/mm10.CGI.bed -o $wdir/$key.gz --outFileNameMatrix $wdir/$key.tab --maxThreshold 10 \
-m 1000 -a 2000 -b 2000 -bs 50 --sortRegions descend --sortUsing mean --sortUsingSamples 1 -p 32 
plotProfile -m $wdir/$key.gz -out $wdir/$key.profile.pdf --yMin 0 --dpi 360 --startLabel "" --endLabel "" -y "mean RPKM" -T "K9me3 to CGI" --perGroup &
plotHeatmap -m $wdir/$key.gz -out $wdir/$key.order-mean.pdf --dpi 360 --outFileSortedRegions $wdir/$key.order-mean.tab --sortRegions descend --sortUsing mean --colorList "blue,white,red" --colorNumber 100000 --whatToShow "plot, heatmap and colorbar" --zMin -3 --zMax 3 --startLabel "" --endLabel "" -T "K9me3 to CGI" &
