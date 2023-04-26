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


for i in *bed; do k=${i/.bed/.motif}
nohup findMotifsGenome.pl $i mm10 $k -size 500 -p 16 &
done

awk 'BEGIN{print "Cluster\tMotif\tPvalue"} FNR>1{split(FILENAME,temp,"-");split($3,temp2,"e");print temp[1],$1,temp2[2]}' *.motif/knownResults.txt > merge.knownResults.txt
awk 'BEGIN{print "Cluster\tMotif\tPvalue\tTarget"} FNR>1{split(FILENAME,temp,"-");split($3,temp2,"e");print temp[1],$1,temp2[2],$7}' *.motif/knownResults.txt > merge2.knownResults.txt
