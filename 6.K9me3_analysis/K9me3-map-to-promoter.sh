wdir=~/workspace/9.NT-ChIP/6.K9me3_analysis/4.promoter_map/6.log2ratio-before2cell-map
pdir=~/ann/Uniq.mm10/all.gene/promoter_1k-1k.bed
ddir=~/workspace/9.NT-ChIP/3.deeptools/5.compareByInput/7.log2-rpkm
mkd $wdir

key=all
computeMatrix scale-regions -S $ddir/CC*bw $ddir/6h*bw $ddir/14h*bw $ddir/NT-2cell*bw $ddir/O*bw $ddir/S*bw $ddir/PN*bw $ddir/NF-2cell*bw --smartLabels \
-R $pdir -o $wdir/$key.gz --outFileNameMatrix $wdir/$key.tab \
-m 2000 -a 2000 -b 2000 -bs 50 --sortRegions descend --sortUsing mean --sortUsingSamples 1 -p 32 &
wait
k=5
plotHeatmap -m $wdir/$key.gz -out $wdir/$key.kmeans$k.png --dpi 360 --kmeans $k --outFileSortedRegions $wdir/$key.kmeans$k.tab \
--sortRegions descend --sortUsing mean --colorList "blue,white,red" --colorNumber 100000 --whatToShow "plot, heatmap and colorbar" --zMin -3 --zMax 3 \
--startLabel "" --endLabel "" -T "K9me3 to promoter" &



sort -k9,9 -k8gr,8 ~/ann/promoter.CpG-density.class/mm10_refseq_promoterclass_NM_number_gene_ID_5most.bed > $wdir/mm10.CGd.promoter.bed
key=CGd.promoter
computeMatrix scale-regions -S $ddir/CC*bw $ddir/6h*bw $ddir/14h*bw $ddir/NT-2cell*bw $ddir/O*bw $ddir/S*bw $ddir/PN*bw $ddir/NF-2cell*bw --smartLabels \
-R $pdir -o $wdir/$key.gz --outFileNameMatrix $wdir/$key.tab \
-m 2000 -a 2000 -b 2000 -bs 50 --sortRegions keep -p 32 &
wait
plotHeatmap -m $wdir/$key.gz -out $wdir/$key.png --dpi 360 --outFileSortedRegions $wdir/$key.tab \
--sortRegions keep --colorList "blue,white,red" --colorNumber 100000 --whatToShow "plot, heatmap and colorbar" --zMin -3 --zMax 3 \
--startLabel "" --endLabel "" -T "K9me3 to promoter" &

