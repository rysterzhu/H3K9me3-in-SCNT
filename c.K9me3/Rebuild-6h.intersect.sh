#1. TSS -1k~1k区域，K9me3 peaks占%33(coverage 1kb)以上为K9me3 markered promoter
#2. 6h markered promoter并且CC markered promoter为6h-inherit
#3. 6h markered promoter并且CC unmarkered promoter为6h-rebuild
#4. 计算两者分别与high-CG promoter和low-CG promoter的交集和补集数量

wdir=~/workspace/9.NT-ChIP/c.K9me3/2.6h-rebuild/8.promoter-with-peaks
ddir=~/workspace/9.NT-ChIP/c.K9me3/1.MACS2-60M-P3/filter-rpkm1
cd $wdir
#1
for i in $ddir/*broadPeak; do
coverageBed -a ~/ann/Uniq.mm10/all.gene/promoter_1k-1k.bed -b $i | awk '$NF>=0.3333{print $4}' > $wdir/$(basename $i _peaks.broadPeak).promoter &
done
wait
#2,3
awk 'NR==FNR{a[$1]} NR>FNR{type=($1 in a)?"inherit":"rebuild";print type,$1}' CC.promoter 6h.promoter > flags-6h.promoter
cut -f 9,10 ~/ann/promoter.CpG-density.class/mm10_refseq_promoterclass_NM_number_gene_ID_5most.bed | grep -v "intCG" > flags-CG.promoter

#4
cat flags-6h.promoter flags-CG.promoter | sort -k2,2 | datamash -g 2 collapse 1 | datamash -s -g 2 count 1

#5.rebuild gene和Oocyte gene的GO
awk 'NR==FNR{print "Oocyte",$1} NR>FNR{print $0}' Oocyte.promoter flags-6h.promoter > flags-oocyte.promoter

#6. CC marked promoter, 6h ummarked为6h-erase
awk 'NR==FNR{a[$1]} NR>FNR{if(!($1 in a))print "erase",$1}' 6h.promoter CC.promoter | cat flags-oocyte.promoter - > flags-oocyte2.promoter

clusterPM flags-oocyte2.promoter flags-oocyte2 GO BP &
Rscript ~/R/clusterProfiler/plot-dot-clusterProfiler.R

awk 'ARGIND==1{a[$1]} ARGIND==2{b[$1]} ARGIND==3{print "PN3",$1} 
END{for(i in b){print (i in a)?"Inherited":"Denovo",i};for(i in a){if(!(i in b))print "Erased",i}}' CC.promoter 6h.promoter PN3.promoter > flags-PN3.promoter
clusterPM flags-PN3.promoter flags-PN3 GO BP &
Rscript ~/R/clusterProfiler/plot-dot-clusterProfiler.R

#########################2021年7月17日 CC-6h K9me3 in Oocyte Sperm Zygote
wdir=~/workspace/9.NT-ChIP/c.K9me3/2.6h-rebuild/9.CC-6h_K9me3_in_OSZ
ddir=~/workspace/9.NT-ChIP/c.K9me3/1.MACS2-60M-P3/filter-rpkm1
ddir2=~/workspace/9.NT-ChIP/3.deeptools/9.multiBamSummary/1.merge-rep

echo -e "peak\tsample\tvalue" > merge.log2ratio
for i in Oocyte Sperm PN3 PN5; do 
for j in CC 6h; do 
awk -v i=$i -v j=$j 'NR==FNR{a[$1"\t"$2]} NR>FNR&&($1"\t"$2 in a){print i,j,$4}' $ddir/$i.fraction $ddir2/$j.1k.m1.log2ratio >> merge.log2ratio
done;done


ddir=~/workspace/9.NT-ChIP/1.align/d.sih2-K9me3/1.mm10/c.merge_rep/rep34
mkcd ~/workspace/9.NT-ChIP/c.K9me3/2.6h-rebuild/6.sih2/merge_rep

awk 'NR==FNR{a[$2]=$1} NR>FNR{if($4 in a){print $0,a[$4]}}' ../../8.promoter-with-peaks/flags-oocyte2.promoter ~/ann/Uniq.mm10/all.gene/promoter_1k-1k.bed > merge.promoter.bed
for i in control sih2; do
nohup rpkmSummary -b merge.promoter.bed -t $ddir/6hpa_${i}.redun.bam -o ./$i --id 7 -p 8 -v --deeptools "--ignoreDuplicates --samFlagInclude 2" > $i.log &
done;wait
awk 'FNR==1{gsub(".rpkm.bg","",FILENAME)} FNR>1{print FILENAME,$5,$4}' *rpkm.bg > merge34.rpkm.tab

