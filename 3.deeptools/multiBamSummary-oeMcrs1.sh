#######################################################################################
#############################NT oeMcrs1 NChIP 
ddir=~/workspace/9.NT-ChIP/1.align/h.K9me3_oeMcrs1/c.merge_rep
wdir=~/workspace/9.NT-ChIP/3.deeptools/9.multiBamSummary/8.oeMcrs1-K9me3
mkd $wdir/logs;cd $wdir
conditions=("NT-ICM_oeMcrs1" "NT-TE_oeMcrs1")
mult 10
for i in ${conditions[@]}; do for k in 10k 1k chrom; do read -u 9; { key=$i.$k
nohup multiBamSummary BED-file --BED ~/ann/makewindows/mm10.$k.bed -b $ddir/${i}_K9me3.redun.bam $ddir/${i}_input.redun.bam -o $wdir/logs/$key.gz --smartLabels -p 8 --centerReads --outRawCounts $wdir/$key.tab --ignoreDuplicates --samFlagInclude 2 --samFlagExclude 1024 > $wdir/logs/$key.log
echo >&9;}&
done;done;wait-m;for i in $wdir/*chrom.tab; do sed -i "1s/['#]//g" $i; sed -i "1s/.sorted//g" $i; sed -i "1s/.redun//g" $i; done

k=promoter_10k-10k
bed=~/ann/Uniq.mm10/all.gene/$k.bed
for i in ${conditions[@]}; do key=$i.$k
nohup multiBamSummary BED-file --BED $bed -b $ddir/${i}_K9me3.redun.bam $ddir/${i}_input.redun.bam -o $wdir/logs/$key.gz --smartLabels -p 8 --centerReads --ignoreDuplicates --samFlagInclude 2 --samFlagExclude 1024 --outRawCounts $wdir/$key.tab > $wdir/logs/$key.log &
done;wait
mult 20
for i in ${conditions[@]}; do for k in 10k 1k chrom promoter_10k-10k; do read -u 9; { 
awk 'NR==1{print;next} NR==FNR{for(i=4;i<=NF;i++)a[i]+=$i} NR>FNR&&FNR>1{for(i=4;i<=NF;i++)$i=($i*1e9)/(($3-$2)*a[i]);print $0}' $wdir/$i.chrom.tab $wdir/$i.$k.tab > $wdir/$i.$k.rpkm 
awk -v m=1 'NR==1{NF-=1;print $0;next} NR!=1{for(i=4;i<NF;i++) $i=(log(($i+m)/($NF+m))/log(2));NF-=1;print $0}' $wdir/$i.$k.rpkm > $wdir/$i.$k.m1.log2ratio
echo>&9; }&done;done;wait
#merge

ln -s ~/workspace/9.NT-ChIP/3.deeptools/9.multiBamSummary/1.merge-rep/*ICM* ./
ln -s ~/workspace/9.NT-ChIP/3.deeptools/9.multiBamSummary/1.merge-rep/*TE* ./
conditions=("NF-ICM" "NF-TE" "NT-ICM" "NT-TE" "NT-ICM_oeMcrs1" "NT-TE_oeMcrs1")
for k in 1k.m1.log2ratio 1k.rpkm;do
cut -f 1,2,3 NF-ICM.$k | sorth -k1V,1 -k2n,2 -k3n,3 > bed.temp &
for i in ${conditions[@]}; do sorth -k1V,1 -k2n,2 -k3n,3 < $i.$k | cut -f 4 > $i.temp & done;wait
paste `for i in bed ${conditions[@]}; do echo -n $i.temp" ";done` > merge.$k
rm *temp
done 
#merge promoter
for k in promoter_10k-10k.m1.log2ratio promoter_10k-10k.rpkm;do
awk 'BEGIN{print "chr\tstart\tend\tgene"} NR==FNR&&FNR>1{a[$1"\t"$2]} NR>FNR&&($1"\t"$2 in a){print $1,$2,$3,$4}' NF-ICM.$k ~/ann/Uniq.mm10/all.gene/promoter_10k-10k.bed | sorth -k1V,1 -k2n,2 -k3n,3 -k4,4 > bed.temp &
for i in ${conditions[@]}; do sort -k1V,1 -k2n,2 -k3n,3 -k4,4 $i.$k | cut -f 4 > $i.temp & done;wait 
paste `for i in bed ${conditions[@]}; do echo -n $i.temp" ";done` > merge.$k
wc -l *temp
rm *temp
done


wdir=~/workspace/9.NT-ChIP/3.deeptools/9.multiBamSummary/9.oeMcrs1-reps
ddir=~/workspace/9.NT-ChIP/1.align/h.K9me3_oeMcrs1/b.filtered
ddir2=~/workspace/9.NT-ChIP/1.align/h.K9me3_oeMcrs1/c.merge_rep
mkd $wdir/logs

mult 10
for i in NT-ICM_oeMcrs1 NT-TE_oeMcrs1; do 
for k in chrom 1k 10k; do key=$i.$k
nohup multiBamSummary BED-file --BED ~/ann/makewindows/mm10.$k.bed -b $ddir/$i*K9me3*.sorted.bam $ddir2/${i}_input.redun.bam -o $wdir/logs/$key.gz --smartLabels -p 8 --centerReads --outRawCounts $wdir/$key.tab --samFlagInclude 2 --samFlagExclude 1036 > $wdir/logs/$key.log &
done;done;wait
for i in $wdir/*chrom.tab; do sed -i "1s/['#]//g" $i; sed -i "1s/.sorted//g" $i; sed -i "1s/.redun//g" $i; done

for k in 1k 10k; do
for i in NT-ICM_oeMcrs1 NT-TE_oeMcrs1; do key=$i.$k
awk 'NR==1{print;next} NR==FNR{for(i=4;i<=NF;i++)a[i]+=$i} 
NR>FNR&&FNR>1{for(i=4;i<=NF;i++)$i=($i*1e9)/(($3-$2)*a[i]);print $0}' $wdir/$i.chrom.tab $wdir/$i.$k.tab > $wdir/$i.$k.rpkm &
done;done;wait

for k in 1k 10k; do
for i in NT-ICM_oeMcrs1 NT-TE_oeMcrs1; do
awk -v m=1 'NR==1{NF-=1;print $0;next} NR!=1{for(i=4;i<NF;i++) $i=(log(($i+m)/($NF+m))/log(2));NF-=1;print $0}' $wdir/$i.$k.rpkm > $wdir/$i.$k.m1.log2ratio &
done
done
 
cat ~/workspace/9.NT-ChIP/c.K9me3/1.MACS2-60M-P3/*[IE]*broadPeak | sort -k1,1 -k2n,2 | mergeBed -i -| awk '{a=int($2/1000+0.5);b=int($3/1000+0.5)
for(i=a;i<b;i++){print $1,i*1000,i*1000+1000}
}' - | sort -k1,1 -k2n,2 | uniq > merge.fraction
