wdir=~/workspace/9.NT-ChIP/c.K9me3/3.domain_cluster/a.CC-inherit-in-filter-NF/1.1k
pdir=~/workspace/9.NT-ChIP/c.K9me3/1.MACS2-80M-P3
ddir=~/workspace/9.NT-ChIP/3.deeptools/9.multiBamSummary/1.merge-rep
cat $pdir/PN*fraction $pdir/NF-2cell.fraction | awk 'NR==FNR{a[$1"\t"$2]} NR>FNR&&!($1"\t"$2 in a){print $0}' - $pdir/CC.fraction > target.bed
awk 'NR==FNR{a[$1"\t"$2]} NR>FNR&&(($1"\t"$2 in a)||FNR==1){print $0} ' target.bed $ddir/merge.1k.m1.log2ratio > target.log2ratio.tab

#kdm4b heatmap
awk 'NR==FNR{a[$1"\t"$2]} NR>FNR&&(($1"\t"$2 in a)||FNR==1){print $0} ' target.bed ~/workspace/9.NT-ChIP/3.deeptools/9.multiBamSummary/7.kdm4b-K9me3/merge.1k.m1.log2ratio > target.kdm4b.log2ratio.tab
sed -i 's/_K9me3//g' target.kdm4b.log2ratio.tab









