ddir=~/workspace/9.NT-ChIP/3.deeptools/1.bamCoverage/7.K9me3
wdir=~/workspace/9.NT-ChIP/3.deeptools/2.map_to_tss/3.K9me3
(n=20190303
computeMatrix reference-point -S $ddir/*bw -R ~/ann/mm10.RefSeq.bed -p 32 -a 3000 -b 3000 --referencePoint TSS -bs 50 --skipZeros -o $wdir/$n.gz --smartLabels 
plotHeatmap -m $wdir/$n.gz -out $wdir/$n.png --dpi 360 --refPointLabel TSS \
--sortRegions descend --sortUsing mean --colorList "white,blue" --whatToShow "plot, heatmap and colorbar" \
--refPointLabel TSS --regionsLabel "Sort by all" -T "") &