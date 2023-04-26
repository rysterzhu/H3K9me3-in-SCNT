ddir1=~/workspace/9.NT-ChIP/4.callPeaks/4.K9me3/g.MACS2-broad-default-merge-downsample-80M/P3
ddir2=~/workspace/9.NT-ChIP/2.public/6.callPeaks/3.K9me3/f.MACS2-broad-default-merge-downsample-80M/P3
wdir=~/workspace/9.NT-ChIP/5.Peaks_analysis/2.K9me3/7.peaks-trans/2.4-cluster
cat $ddir2/[MzP2]*broadPeak | sort -k1,1 -k2n,2 - | mergeBed -i - | subtractBed -a $ddir1/cc_peaks.broadPeak -b - | \
intersectBed -a - -b $ddir1/6h_peaks.broadPeak | intersectBed -a - -b $ddir1/14h_peaks.broadPeak | \
intersectBed -a - -b $ddir1/2cell_peaks.broadPeak > 2cell-unreprog.bed &

cat $ddir2/[MzP2]*broadPeak $ddir1/2cell_peaks.broadPeak | sort -k1,1 -k2n,2 - | mergeBed -i - | subtractBed -a $ddir1/cc_peaks.broadPeak -b - | \
intersectBed -a - -b $ddir1/6h_peaks.broadPeak | intersectBed -a - -b $ddir1/14h_peaks.broadPeak > 14h-unreprog.bed &

cat $ddir2/[MzP2]*broadPeak $ddir1/2cell_peaks.broadPeak $ddir1/14h_peaks.broadPeak | sort -k1,1 -k2n,2 - | mergeBed -i - | \
subtractBed -a $ddir1/cc_peaks.broadPeak -b - | intersectBed -a - -b $ddir1/6h_peaks.broadPeak > 6h-reprog.bed &

cat $ddir2/[MzP2]*broadPeak $ddir1/2cell_peaks.broadPeak $ddir1/14h_peaks.broadPeak $ddir1/6h_peaks.broadPeak | sort -k1,1 -k2n,2 - | mergeBed -i - | \
subtractBed -a $ddir1/cc_peaks.broadPeak -b - > cc-reprog.bed &
for i in *bed; do
awk '{a+=1;b+=$3-$2} END{print FILENAME,a,b/27e8,b/a}' $i
done

14h-unreprog.bed	25370	0.00263708	280.651
2cell-unreprog.bed	52241	0.00925651	478.409
6h-reprog.bed	33591	0.0029258	235.172
cc-reprog.bed	158011	0.0245769	419.955

for i in *bed; do k=${i/.bed/.motif}
nohup findMotifsGenome.pl $i mm10 $k -size 500 -p 16 &
done
#plot

awk 'BEGIN{print "Cluster\tMotif\tPvalue"} FNR>1{split(FILENAME,temp,"-");split($3,temp2,"e");print temp[1],$1,temp2[2]}' *.motif/knownResults.txt > merge.knownResults.txt
awk 'BEGIN{print "Cluster\tMotif\tPvalue\tTarget"} FNR>1{split(FILENAME,temp,"-");split($3,temp2,"e");print temp[1],$1,temp2[2],$7}' *.motif/knownResults.txt > merge2.knownResults.txt
