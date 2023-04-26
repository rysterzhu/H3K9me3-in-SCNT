wdir=~/workspace/9.NT-ChIP/c.K9me3/3.domain_cluster/a.CC-inherit-in-filter-NF/1.1k
pdir=~/workspace/9.NT-ChIP/c.K9me3/1.MACS2-80M-P3
ddir=~/workspace/9.NT-ChIP/3.deeptools/9.multiBamSummary/1.merge-rep
cat $pdir/PN*fraction $pdir/NF-2cell.fraction | awk 'NR==FNR{a[$1"\t"$2]} NR>FNR&&!($1"\t"$2 in a){print $0}' - $pdir/CC.fraction > target.bed
awk 'NR==FNR{a[$1"\t"$2]} NR>FNR&&(($1"\t"$2 in a)||FNR==1){print $0} ' target.bed $ddir/merge.1k.m1.log2ratio > target.log2ratio.tab

#kdm4b heatmap
awk 'NR==FNR{a[$1"\t"$2]} NR>FNR&&(($1"\t"$2 in a)||FNR==1){print $0} ' target.bed ~/workspace/9.NT-ChIP/3.deeptools/9.multiBamSummary/7.kdm4b-K9me3/merge.1k.m1.log2ratio > target.kdm4b.log2ratio.tab
sed -i 's/_K9me3//g' target.kdm4b.log2ratio.tab

results.NT-oe.tab
#将就
ddir3=~/workspace/9.NT-ChIP/2.public/c.Xiewei/5.RNA/4.hisat2/3.bg-1k
ddir2=~/workspace/9.NT-ChIP/1.align/9.RNA_hisat2/3.bg-1k

tail -n +2 cluster-kdm4b.bed|intersectBed -a ~/ann/Uniq.mm10/all.gene/promoter_2k-2k.bed -b - -wao | awk '$NF>0' | cut -f 4,10,11 | sort -k1,1 -k3n,3 -k2nr,2 | awk '{a[$1]=$2} END{for(i in a)print i,a[i]}' > cluster.gene.tab

wdir=~/workspace/9.NT-ChIP/c.K9me3/3.domain_cluster/a.CC-inherit-in-filter-NF/2.500bp
pdir=~/workspace/9.NT-ChIP/c.K9me3/1.MACS2-80M-P3/2.500bp-frations
ddir=~/workspace/9.NT-ChIP/3.deeptools/9.multiBamSummary/1.merge-rep/2.500
cat $pdir/Oocyte.fraction $pdir/Sperm.fraction $pdir/PN*fraction $pdir/NF-2cell.fraction | awk 'NR==FNR{a[$1"\t"$2]} NR>FNR&&!($1"\t"$2 in a){print $0}' - $pdir/CC.fraction > target.bed
awk 'NR==FNR{a[$1"\t"$2]} NR>FNR&&(($1"\t"$2 in a)||FNR==1){print $0} ' target.bed $ddir/merge.log2ratio > target.log2ratio.tab


awk 'NR>1{print $1,$2,$3 > "motifs/Cluster"$4".bed"}' cluster5.bed; cd motifs
for i in Cluster*bed; do
nohup findMotifsGenome.pl $i mm10 ${i/bed/motif} -size 500 -mask -p 8 &  done
#motif 仍然不行









