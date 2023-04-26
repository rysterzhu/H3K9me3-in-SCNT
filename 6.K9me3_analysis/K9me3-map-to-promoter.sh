key=all
computeMatrix scale-regions -S $ddir/CC*bw $ddir/6h*bw $ddir/14h*bw $ddir/NT-2cell*bw $ddir/O*bw $ddir/S*bw $ddir/PN*bw $ddir/NF-2cell*bw --smartLabels \
-R $pdir -o $wdir/$key.gz --outFileNameMatrix $wdir/$key.tab \
-m 2000 -a 2000 -b 2000 -bs 50 --sortRegions descend --sortUsing mean --sortUsingSamples 1 -p 32 &
wait

sort -k9,9 -k8gr,8 ~/ann/promoter.CpG-density.class/mm10_refseq_promoterclass_NM_number_gene_ID_5most.bed > $wdir/mm10.CGd.promoter.bed
key=CGd.promoter
computeMatrix scale-regions -S $ddir/CC*bw $ddir/6h*bw $ddir/14h*bw $ddir/NT-2cell*bw $ddir/O*bw $ddir/S*bw $ddir/PN*bw $ddir/NF-2cell*bw --smartLabels \
-R $pdir -o $wdir/$key.gz --outFileNameMatrix $wdir/$key.tab \
-m 2000 -a 2000 -b 2000 -bs 50 --sortRegions keep -p 32 &
wait
plotHeatmap -m $wdir/$key.gz -out $wdir/$key.png --dpi 360 --outFileSortedRegions $wdir/$key.tab \
--sortRegions keep --colorList "blue,white,red" --colorNumber 100000 --whatToShow "plot, heatmap and colorbar" --zMin -3 --zMax 3 \
--startLabel "" --endLabel "" -T "K9me3 to promoter" &
