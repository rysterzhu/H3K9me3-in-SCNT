cols=("#540d6eff" "#a1286aff" "#ee4266ff" "#ffd23fff" "#9dd076ff" "#25be8bff" "#0ead69ff")

s1=NT-TE;s2=NT-ICM
echo -e "TF\toverlap1\toverlap2\ttotal1\ttotal2" > overlap-number.tab
for i in ~/workspace/9.NT-ChIP/3.deeptools/c.TF-blastocyst/mergedpeak/*bed; do 
k=$(basename $i .bed); 
(o1=$(intersectBed -a $i -b $s1.bed -r -e -f 0.5 -wa -u | wc -l); 
o2=$(intersectBed -a $i -b $s2.bed -r -e -f 0.5 -wa -u | wc -l); 
t1=$(cat $s1.bed |wc -l); 
t2=$(cat $s2.bed | wc -l); 
echo -e $k"\t"$o1"\t"$o2"\t"$t1"\t"$t2 >> overlap-number.tab)& 
done

awk 'NR>1&&FNR<NR{if(substr($1,1,1)""tolower(substr($1,2)) in a){print $1,substr($1,1,1)""tolower(substr($1,2))}} NR==FNR{a[$4]}' ~/ann/Uniq.mm10/mm10.uniq-gene.bed overlap-number.tab > ../TF2gene.tab

cd ~/workspace/9.NT-ChIP/3.deeptools/c.TF-blastocyst/1.fisher-test 
awk -v s="NF-ICM.NF-TE" 'FNR==1{next}
ARGIND==1{g[$2]=$1;t[$1]=$2}
(ARGIND==2)&&($8==s)&&($7 in g){exp1[$7]=$2"\t"$5"\t"$6}
(ARGIND==3)&&($8 in t){res[$8]=$7"\t"$5"\t"$6}
END{print "gene\tTF\tlog2FC.exp\tpvalue.exp\tpadj.exp\tlog2FC.TF\tpvalue.TF\tlog10padj.TF";
for(i in g){print i,g[i],exp1[i],res[g[i]]}}' TF2gene.tab ~/workspace/9.NT-ChIP/b.RNA/b.DESeq2-blastocyst/results.NT-NF.tab \
1.NF-ICM-TE/results.tab > 1.NF-ICM-TE/exp.TF.pvalue.tab

awk -v s="NT-ICM.NF-ICM" 'FNR==1{next}
ARGIND==1{g[$2]=$1;t[$1]=$2}
(ARGIND==2)&&($8==s)&&($7 in g){$2=-$2;exp1[$7]=$2"\t"$5"\t"$6}
(ARGIND==3)&&($8 in t){res[$8]=$7"\t"$5"\t"$6}
END{print "gene\tTF\tlog2FC.exp\tpvalue.exp\tpadj.exp\tlog2FC.TF\tpvalue.TF\tlog10padj.TF";
for(i in g){print i,g[i],exp1[i],res[g[i]]}}' TF2gene.tab ~/workspace/9.NT-ChIP/b.RNA/b.DESeq2-blastocyst/results.NT-NF.tab \
4.NF-NT-ICM/results.tab > 4.NF-NT-ICM/exp.TF.pvalue.tab


awk -v s="NT-TE.NF-TE" 'FNR==1{next}
ARGIND==1{g[$2]=$1;t[$1]=$2}
(ARGIND==2)&&($8==s)&&($7 in g){$2=-$2;exp1[$7]=$2"\t"$5"\t"$6}
(ARGIND==3)&&($8 in t){res[$8]=$7"\t"$5"\t"$6}
END{print "gene\tTF\tlog2FC.exp\tpvalue.exp\tpadj.exp\tlog2FC.TF\tpvalue.TF\tlog10padj.TF";
for(i in g){print i,g[i],exp1[i],res[g[i]]}}' TF2gene.tab ~/workspace/9.NT-ChIP/b.RNA/b.DESeq2-blastocyst/results.NT-NF.tab \
2.NF-NT-TE/results.tab > 2.NF-NT-TE/exp.TF.pvalue.tab

awk -v s="NT-ICM.NT-TE" 'FNR==1{next}
ARGIND==1{g[$2]=$1;t[$1]=$2}
(ARGIND==2)&&($8==s)&&($7 in g){exp1[$7]=$2"\t"$5"\t"$6}
(ARGIND==3)&&($8 in t){res[$8]=$7"\t"$5"\t"$6}
END{print "gene\tTF\tlog2FC.exp\tpvalue.exp\tpadj.exp\tlog2FC.TF\tpvalue.TF\tlog10padj.TF";
for(i in g){print i,g[i],exp1[i],res[g[i]]}}' TF2gene.tab ~/workspace/9.NT-ChIP/b.RNA/b.DESeq2-blastocyst/results.NT-NF.tab \
3.NT-ICM-TE/results.tab > 3.NT-ICM-TE/exp.TF.pvalue.tab

awk 'NR==1{print $0,"constrast"} FNR>1{split(FILENAME,temp,"[\\./]");print $0,temp[2]}' */exp.TF.pvalue.tab > all.exp.TF.qvalue.tab

mkd annStats
awk '{print $1,$2,$3 > $7".bed"}' cluster-ICM-TE.tab
for i in NT NF; do 
nohup annotatePeaks.pl $i.bed mm10 -annStats annStats/$i.stat > annStats/$i.txt 2> annStats/$i.log &
done;wait;
cd annStats
awk 'NR==1{print $0,"sample"} FNR==1{gsub(".stat","",FILENAME)} FNR>14&&$1!~/Annotation/&&$1!~/\?/{print $0,FILENAME}' *stat  > all.stats.tab

for i in NT NF; do 
intersectBed -a ~/ann/Uniq.mm10/all.gene/promoter_2k-2k.bed -b ../$i.bed -wa | awk -v i=$i '{print i,$4}'
done | sort | uniq > cluster.gene

clusterPM cluster.gene cluster.BP GO BP


intersectBed -a ~/ann/Uniq.mm10/all.gene/promoter_2k-2k.bed -b ../cluster-ICM-TE.tab -wo |sort -k12rg,12 > cluster-ICM-TE.promoter.tab
#############################
for i in Cluster*bed; do nohup annotatePeaks.pl $i mm10 -annStats annStat/${i/bed/stat} > annStat/${i/bed/txt} 2> annStat/${i/bed/log} & done
awk 'NR==1{print $0,"sample"} FNR==1{gsub(".stat","",FILENAME)} FNR>14&&$1!~/Annotation/&&$1!~/\?/{print $0,FILENAME}' *stat  > all.stats.tab

tail -n +2 Cluster4.bed | cat Cluster6.bed - | tail -n +2 |cat Cluster3.bed - | tail -n +2 | sort -k1,1 -k2n,2  > Down.bed
tail -n +2 Cluster5.bed > Up.bed
for i in Down.bed Up.bed; do nohup annotatePeaks.pl $i mm10 -annStats annStat2/${i/bed/stat} > annStat2/${i/bed/txt} 2> annStat2/${i/bed/log} & done



