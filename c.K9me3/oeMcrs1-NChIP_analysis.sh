
s1=NT-ICM_oeMcrs1;s2=NT-TE_oeMcrs1
echo -e "TF\toverlap1\toverlap2\ttotal1\ttotal2" > overlap-number.tab
for i in ~/workspace/9.NT-ChIP/3.deeptools/c.TF-blastocyst/mergedpeak/*bed; do 
k=$(basename $i .bed); 
(o1=`intersectBed -a $i -b $s1.bed -r -e -f 0.5 -wa -u | wc -l`; 
o2=`intersectBed -a $i -b $s2.bed -r -e -f 0.5 -wa -u | wc -l`; 
t1=`cat $s1.bed |wc -l`; 
t2=`cat $s2.bed | wc -l`; 
echo -e $k"\t"$o1"\t"$o2"\t"$t1"\t"$t2 >> overlap-number.tab)& 
done

plot-point-TF.Rmd

for i in NF NT oeMcrs1; do awk -v i=$i 'FNR>1{print i,$8,$6,$7}' ${i}-ICM.${i}-TE/results.tab >> oeMcrs1-ICM.oeMcrs1-TE/merge.result.ta
b ;done




##################################################expression in peaks
subtractBed -a $pdir/NT-ICM_oeMcrs1_peaks.broadPeak -b $pdir/NT-ICM_peaks.broadPeak | awk '$3-$2>200'> $wdir/NT-ICM_oeMcrs1.subtract.peak
subtractBed -b $pdir/NT-ICM_oeMcrs1_peaks.broadPeak -a $pdir/NT-ICM_peaks.broadPeak | awk '$3-$2>200'> $wdir/NT-ICM.subtract.peak
intersectBed -a $pdir/NT-ICM_oeMcrs1_peaks.broadPeak -b $pdir/NT-ICM_peaks.broadPeak > $wdir/NT-ICM.intersect.peak

for i in *peak; do 
multiBigwigSummary BED-file --BED $i -b $ddir/NT*control* $ddir/NT*Mcrs1* --smartLabels  --outRawCounts ${i/peak}.exp.tab -o ${i/peak}.exp.gz -p 32 &
done

for i in *peak; do nohup annotatePeaks.pl $i mm10 -annStats annStats/${i/peak/txt} > annStats/${i/peak/log} & done


awk 'BEGIN{FS="\t";a=0} FNR==1{a=0;sub(".log","",FILENAME)} /^PeakID/{a=1;next} a==1&&$10<10000&&$10>-10000&&$16!=""{print FILENAME,$16}' *.log | sort -k1,1 -k2,2 | uniq > nearest.10k.tab
awk 'BEGIN{FS="\t";a=0} FNR==1{a=0;sub(".log","",FILENAME)} /^PeakID/{a=1;next} a==1&&$16!=""{print FILENAME,$16,$10}' *.log | sort -k1,1 -k2,2 | uniq > nearest.tab
