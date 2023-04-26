
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
