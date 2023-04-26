wdir=~/workspace/9.NT-ChIP/6.K9me3_analysis/6.Transposal_map/4.merge-MERVL
ddir=~/workspace/9.NT-ChIP/3.deeptools/5.compareByInput/a.keep-lowqual
cd $wdir
awk '$1~/^chr[1-9X]+$/{a=int($2/100)*100;b=int($3/100)*100;print $1,a,a+100,10;print$1,b,b+100,10}' MERVL.bed | sort -k1,1 -k2n,2 | uniq > temp.bg
bedGraphToBigWig temp.bg ~/ann/mm10.shortChromSizes mervl.boundary.bw

key=K9me3-center-sortlen-10k
computeMatrix reference-point -S $ddir/CC*bw $ddir/6h*bw $ddir/14h*bw $ddir/NT-2cell*bw $ddir/NT-4cell*bw $ddir/NT-14h_kdm4b_K9me3.bw $ddir/NT-l2cell_kdm4b_K9me3.bw $ddir/O*bw $ddir/Sperm*bw $ddir/PN*bw $ddir/NF-2cell*bw $ddir/NF-4cell*bw \
 --smartLabels -R $wdir/MERVL.bed -o $wdir/$key.gz --outFileNameMatrix $wdir/$key.tab \
--referencePoint center -a 10000 -b 10000 -bs 50 --sortRegions keep --sortUsing mean -p 32 
plotHeatmap -m $wdir/$key.gz -out $wdir/$key.png --dpi 360 --outFileSortedRegions $wdir/$key.tab \
--sortRegions keep --sortUsing mean --colorList "blue,white,red" --colorNumber 100000 --whatToShow "plot, heatmap and colorbar" --zMin -3 --zMax 3 \
--startLabel "" --endLabel "" -T "K9me3 to MERVL"

key=boundary-center-sortlen-10k
computeMatrix reference-point -S $wdir/mervl.boundary.bw \
 --smartLabels -R $wdir/MERVL.bed -o $wdir/$key.gz --outFileNameMatrix $wdir/$key.tab \
--referencePoint center -a 10000 -b 10000 -bs 50 --sortRegions keep --sortUsing mean -p 32 
plotHeatmap -m $wdir/$key.gz -out $wdir/$key.png --dpi 360 --outFileSortedRegions $wdir/$key.tab \
--sortRegions keep --sortUsing mean --colorList "blue,white,red" --colorNumber 100000 --whatToShow "plot, heatmap and colorbar" --zMin -3 --zMax 3 \
--startLabel "" --endLabel "" -T "boundary to MERVL"













