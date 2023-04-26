wdir=~/workspace/9.NT-ChIP/3.deeptools/c.TF-blastocyst/3.each-TF-K9me3/c.log2ratio
pdir=~/workspace/9.NT-ChIP/3.deeptools/c.TF-blastocyst/mergedpeak
ddir=~/workspace/9.NT-ChIP/e.Keep_lowqual_K9me3/1.merge-reps
ddir2=~/workspace/9.NT-ChIP/1.align/h.K9me3_oeMcrs1/c.merge_rep
ddir3=~/workspace/9.NT-ChIP/1.align/i.K9me3_sertoli/c.merge_rep
mkd $wdir/logs
#计算TF的peaks上的K9me3信号
mult 20
for i in $pdir/*bed; do read -u 9; { 
key=$(basename $i .bed); 
nohup multiBamSummary BED-file --BED $i -b $ddir/N*[IT][CE]*redun.bam $ddir2/*oeMcrs1*redun.bam $ddir3/*redun.bam -o $wdir/logs/$key.gz --smartLabels -p 4 --centerReads --outRawCounts $wdir/$key.tab > $wdir/logs/$key.log 
echo>&9; }&done;
nohup multiBamSummary BED-file --BED ~/ann/makewindows/mm10.chrom.bed -b $ddir/N*[IT][CE]*redun.bam $ddir2/*oeMcrs1*redun.bam $ddir3/*redun.bam -o $wdir/logs/mm10.chrom.gz --smartLabels -p 4 --centerReads --outRawCounts $wdir/mm10.chrom.tab > $wdir/logs/mm10.chrom.log & 
wait
sed -i "1s/[#']//g" mm10.chrom.tab;sed -i "1s/.redun//g" mm10.chrom.tab;
for i in *tab; do cut -f -15,17- $i > $i.temp & done 
#rpkm
for i in $pdir/*bed; do
k=$(basename $i .bed); 
(awk 'NR==1{print;next} NR==FNR{for(i=4;i<=NF;i++)a[i]+=$i} NR>FNR&&FNR>1{for(i=4;i<=NF;i++)$i=($i*1e9)/(($3-$2)*a[i]);print $0}' $wdir/1.tabs/mm10.chrom.tab $wdir/1.tabs/$k.tab > $wdir/$k.rpkm 
awk -v m=1 'NR==1{out=$1"\t"$2"\t"$3;for(i=4;i<NF;i+=2){out=out"\t"$(i+1)};print out} NR!=1{out=$1"\t"$2"\t"$3;for(i=4;i<NF;i+=2){r=(log(($(i+1)+m)/($i+m))/log(2));out=out"\t"r};print out}' $wdir/$k.rpkm > $wdir/$k.m1.log2ratio) &
done;wait

mkd 1.tabs 2.rpkms 3.ratios;
mv *tab 1.tabs;mv *rpkm 2.rpkms;mv *ratio 3.ratios
mkd 4.filters 5.pdfs
#筛选与NTNF ICMTE相关的peak上的K9me3信号，画boxplot（在这4种细胞中有K9me3 peaks的那些TF peaks）
for i in *log2ratio; do sed -i "1s/_K9me3//g" $i;done 
Rscript ~/R/9.NTChIP/5.Deeptools/plot-boxplot-TF.R 3.ratios/MCRS1.m1.log2ratio MCRS1.m1.pdf

intersectBed -a $i -b ~/workspace/9.NT-ChIP/c.K9me3/b.fraction-blastocyst/1.cluster7-from-pca/merge.fraction -r -e -f 0.5 -wa -u -header > "${i/rpkm/filter}"; done


#计算kruskal test的pvalue和FC
#for i in *filter; do  Rscript ~/R/9.NTChIP/5.Deeptools/dunnTest.R $i ; done
wdir=~/workspace/9.NT-ChIP/3.deeptools/c.TF-blastocyst/3.each-TF-K9me3
pdir=~/workspace/9.NT-ChIP/3.deeptools/c.TF-blastocyst/mergedpeak
ddir=~/workspace/9.NT-ChIP/1.align/6.links-merge-redun
trap "exec 1000>&-;exec 1000<&-;exit 0" 2;mkfifo /tmp/$$.fifo;exec 1000<>/tmp/$$.fifo;for((i=0;i<16;i++));do echo >&1000; done
for i in *filter; do read -u 1000{; Rscript ~/R/9.NTChIP/5.Deeptools/dunnTest.R $i ; echo >&1000}& done

for i in *filter; do Rscript ~/R/9.NTChIP/5.Deeptools/plot-boxplot-TF.R $i ${i/filter/pdf} & done