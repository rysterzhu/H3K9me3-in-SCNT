################################################################################################################################
##############################################K9me3 BWA callpeaks
ddir=~/workspace/9.NT-ChIP/1.align/1.K9me3_BWA/b.filtered
wdir=~/workspace/9.NT-ChIP/3.deeptools/7.correlation/3.K9me3_BWA/1.callPeaks
for i in $ddir/*bam; do k=$(basename $i .sorted.bam)
nohup macs2 callpeak -g mm -n $k -B -q 0.05 --nomodel --broad --shift=73 --SPMR -t $i >> $wdir/logs/$k.log &
done


for i in *_peaks.broadPeak; do k=$(basename $i _peaks.broadPeak)
awk '{a=int($2/1000+0.5);b=int($3/1000+0.5)
for(i=a;i<b;i++){print $1,i*1000,i*1000+1000}
}' $i | sort -k1,1 -k2n,2 | uniq > $k.fraction &
done

cat *fraction| sort|uniq > merge.fraction

NT2=("cc:CC" "6h:6h" "14h:14h" "2cell:NT-2cell" "4cell:NT-4cell" "8cell:NT-8cell" "morula:NT-Morula" "icm:NT-ICM" "te:NT-TE")

n=merge.fraction
(multiBamSummary BED-file --BED $wdir/merge.fraction -b $ddir/[!Pb]*.sorted.bam  -out $wdir/$n.npz -p 48 --outRawCounts $wdir/$n.tab --ignoreDuplicates --centerReads
plotCorrelation -in $wdir/$n.npz -o $wdir/$n.pearson.heatmap.pdf --outFileCorMatrix $wdir/$n.pearson.heatmap.tab -c pearson -p heatmap --skipZeros --removeOutliers --plotNumbers -T "Correlation of K9me3") &



ddir=~/workspace/9.NT-ChIP/1.align/1.K9me3_BWA/b.filtered
wdir=~/workspace/9.NT-ChIP/3.deeptools/7.correlation/3.K9me3_BWA/4.broadPeaks-eachSample

NT2=("cc:CC" "6h:6h" "14h:14h" "2cell:NT-2cell" "4cell:NT-4cell" "8cell:NT-8cell" "morula:NT-Morula" "icm:NT-ICM" "te:NT-TE")

for i in ${NT2[@]}; do k1=${i%%:*};k2=${i##*:};
(multiBamSummary BED-file --BED $pdir/${k2}.fraction -b $ddir/$k1*.sorted.bam  -out $wdir/$k2.npz -p 48 --outRawCounts $wdir/$k2.tab --ignoreDuplicates --centerReads
plotCorrelation -in $wdir/$k2.npz -o $wdir/$k2.pearson.heatmap.pdf --outFileCorMatrix $wdir/$k2.pearson.heatmap.tab -c pearson -p heatmap --skipZeros --removeOutliers --plotNumbers -T "Correlation of K9me3") &
done



pdir=~/workspace/9.NT-ChIP/c.K9me3/1.MACS2-80M-P3
ddir=~/workspace/9.NT-ChIP/1.align/a.links-K9me3-reps
wdir=~/workspace/9.NT-ChIP/3.deeptools/7.correlation/3.K9me3_BWA/5.broadPeaks-merge-fraction
mkd $wdir
cat $pdir/[C16]*fraction $pdir/NT*fraction | sort -k1,1 -k2n,2 | uniq > $wdir/merge.fraction
k2=merge.fraction
multiBamSummary BED-file --BED $wdir/merge.fraction -b $ddir/[C16]*bam $ddir/NT*bam -out $wdir/$k2.npz -p 48 --outRawCounts $wdir/$k2.tab --ignoreDuplicates --centerReads
plotCorrelation -in $wdir/$k2.npz -o $wdir/$k2.pearson.heatmap.pdf --outFileCorMatrix $wdir/$k2.pearson.heatmap.tab -c pearson -p heatmap --skipZeros --removeOutliers --plotNumbers -T "Correlation of K9me3"
plotCorrelation -in $wdir/$k2.npz -o $wdir/$k2.pearson.scatterplot.pdf --outFileCorMatrix $wdir/$k2.pearson.scatterplot.tab -c pearson -p scatterplot --skipZeros --removeOutliers --plotNumbers -T "Correlation of K9me3"

sed "s/.sorted.bam//g" merge.fraction.pearson.heatmap.tab | tail -n +2 > merge.fraction.pearson.heatmap.txt
sed -i "s/_K9me3//g" merge.fraction.pearson.heatmap.txt
sed -i "s/'//g" merge.fraction.pearson.heatmap.txt


#############################################
#hclust K9me3 fraction， 2021年6月2日
#filter rpkm 1
cd ~/workspace/9.NT-ChIP/3.deeptools/b.correlation-NT-NF/a.hclust-K9me3
cat ~/workspace/9.NT-ChIP/c.K9me3/1.MACS2-60M-P3/filter-rpkm1/*fraction | sort -u -k1,1 -k2n,2 -k3n,3 > merge.fraction

awk 'NR==FNR{a[$1"\t"$2"\t"$3]} NR>FNR&&(($1"\t"$2"\t"$3 in a)||FNR==1){print}' merge.fraction ~/workspace/9.NT-ChIP/3.deeptools/9.multiBamSummary/1.merge-rep/merge.1k.rpkm > merge.1k.rpkm
sed -i '1s/_K9me3//g' merge.1k.rpkm

awk 'NR==FNR{a[$1"\t"$2"\t"$3]} NR>FNR&&(($1"\t"$2"\t"$3 in a)||FNR==1){print}' merge.fraction ~/workspace/9.NT-ChIP/3.deeptools/9.multiBamSummary/1.merge-rep/merge.1k.m1.log2ratio > merge.1k.m1.log2ratio
sed -i '1s/_K9me3//g' merge.1k.m1.log2ratio










