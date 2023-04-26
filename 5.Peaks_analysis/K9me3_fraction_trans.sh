
cp $pdir/[61POSC]*.fraction $wdir
cp $pdir/*2cell.fraction $wdir 
cd $wdir
awk '{split(FILENAME,temp,".");a[$0]=a[$0]"\t"temp[1]} END{for(i in a){print i""a[i]}}' *.fraction > marked-sample.bed
awk -v samples="CC 6h 14h NT-2cell Oocyte Sperm PN3 PN5 NF-2cell" 'BEGIN{len=split(samples,S," ");printf("chr\tstart\tend");for(i=1;i<=len;i++)printf("\t%s",S[i]);printf("\n")} 
{printf("%s\t%s\t%s",$1,$2,$3);for(i=1;i<=len;i++){match($0,S[i]);if(RSTART==0){printf("\tunmarked")}else{printf("\tmarked")}};printf("\n")}' marked-sample.bed > all-flags.tab


sed '1d' all-flags.tab | cut -f 4- | sort | uniq -c | awk -v samples="CC 6h 14h NT-2cell Oocyte Sperm PN3 PN5 NF-2cell" 'BEGIN{$0=samples" count";$1=$1;print $0} {for(i=2;i<=NF;i++)printf("%s\t",$i);print $1}' > cdata.tab
awk '$4=="marked"&&$8=="unmarked"&&$9=="unmarked"&&$10=="unmarked"&&$11=="unmarked"&&$12=="unmarked"{print $4,$5,$6,$7}' all-flags.tab |  sort | uniq -c | awk -v samples="CC 6h 14h NT-2cell" 'BEGIN{$0=samples" count";$1=$1;print $0} {for(i=2;i<=NF;i++)printf("%s\t",$i);print $1}' > NT-before-2cell.tab

awk 'NR>1{a="";for(i=4;i<=NF;i++){if($i=="unmarked"){a=a"0"}else{a=a"1"}};
if(a=="100000000"){print $1,$2,$3>"CC-rpg.bed"}
else if(a=="110000000"){print $1,$2,$3>"6h-rpg.bed"}
else if(a=="111000000"){print $1,$2,$3>"14h-unrpg.bed"}
else if(a=="111100000"){print $1,$2,$3>"2cell-unrpg.bed"}
}' all-flags.tab


for i in *rpg.bed; do i=${i/.bed}
nohup findMotifsGenome.pl $i.bed mm10 motifs/$i.motif -size 500 -p 16 &
done

cat *rpg.bed |sort -k1,1 -k2n,2> Rpg.bed
cat *unrpg.bed |sort -k1,1 -k2n,2> UnRpg.bed
nohup findMotifsGenome.pl Rpg.bed mm10 motifs/Rpg.motif -size 500 -p 16 &
nohup findMotifsGenome.pl UnRpg.bed mm10 motifs/UnRpg.motif -size 500 -p 16 &

awk 'BEGIN{print "Cluster\tMotif\tPvalue\tTarget"} FNR>1{split(FILENAME,temp,".");split($3,temp2,"e");gsub("%","",$7);print temp[1],$1,temp2[2],$7}' *rpg.motif/knownResults.txt > merge.cluster4.motif
awk 'BEGIN{print "Cluster\tMotif\tPvalue\tTarget"} FNR>1{split(FILENAME,temp,".");split($3,temp2,"e");gsub("%","",$7);print temp[1],$1,temp2[2],$7}' *Rpg.motif/knownResults.txt > merge.cluster2.motif


#############
#map DNase to peaks #Zhangyi, CC 1cell


key=DNase-bw.cluster2
computeMatrix reference-point -R $wdir/*Rpg.bed -S $ddir/*.bw -out $wdir/$key.gz --referencePoint center \
-a 10000 -b 10000 -bs 10 --sortRegions descend
plotProfile -m $wdir/$key.gz -out $wdir/$key.perGroup.pdf --refPointLabel "regions center" --perGroup
plotProfile -m $wdir/$key.gz -out $wdir/$key.pdf --refPointLabel "regions center"

key=DNase-bw.cluster4
computeMatrix reference-point -R $wdir/CC-rpg.bed $wdir/6h-rpg.bed $wdir/14h-unrpg.bed $wdir/2cell-unrpg.bed -S $ddir/*.bw -out $wdir/$key.gz --referencePoint center \
-a 10000 -b 10000 -bs 10 --sortRegions descend
plotProfile -m $wdir/$key.gz -out $wdir/$key.perGroup.pdf --refPointLabel "regions center" --perGroup
plotProfile -m $wdir/$key.gz -out $wdir/$key.pdf --refPointLabel "regions center"

for i in $wdir/*bed; do o=$(basename $i .bed)
multiBamSummary BED-file -b $ddir/cc_dnase.redun.bam $ddir/1cell_dnase.redun.bam --BED $i -out $wdir/DNase-$o.gz --outRawCounts $wdir/DNase-$o.tab -l CC 1cell -p 64 --extendReads 150 &
done
awk 'FNR>1{print $0,FILENAME}' *tab > DNase.boxplot.txt

