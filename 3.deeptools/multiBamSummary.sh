ddir=~/workspace/9.NT-ChIP/1.align/f.K9me3_keep_lowqual/1.merge-reps
wdir=~/workspace/9.NT-ChIP/3.deeptools/9.multiBamSummary/1.merge-rep
mkd $wdir/logs;cd $wdir
mult 10
for i in ${NTF[@]}; do for k in 10k 1k chrom; do read -u 9; { key=$i.$k
nohup multiBamSummary BED-file --BED ~/ann/makewindows/mm10.$k.bed -b $ddir/${i}_K9me3.redun.bam $ddir/${i}_input.redun.bam -o $wdir/logs/$key.gz --smartLabels -p 8 --centerReads --outRawCounts $wdir/$key.tab --ignoreDuplicates --samFlagInclude 2 --samFlagExclude 1024 > $wdir/logs/$key.log
echo >&9;}&
done;done&wait-m;for i in $wdir/*chrom.tab; do sed -i "1s/['#]//g" $i; sed -i "1s/.sorted//g" $i; sed -i "1s/.redun//g" $i; done

k=promoter_10k-10k
bed=~/ann/Uniq.mm10/all.gene/$k.bed
for i in ${NTF[@]}; do key=$i.$k
nohup multiBamSummary BED-file --BED $bed -b $ddir/${i}_K9me3.redun.bam $ddir/${i}_input.redun.bam -o $wdir/logs/$key.gz --smartLabels -p 8 --centerReads --ignoreDuplicates --samFlagInclude 2 --samFlagExclude 1024 --outRawCounts $wdir/$key.tab > $wdir/logs/$key.log &
done


mult 20
for i in ${NTF[@]}; do for k in 10k 1k chrom; do read -u 9; { 
awk 'NR==1{print;next} NR==FNR{for(i=4;i<=NF;i++)a[i]+=$i} NR>FNR&&FNR>1{for(i=4;i<=NF;i++)$i=($i*1e9)/(($3-$2)*a[i]);print $0}' $wdir/$i.chrom.tab $wdir/$i.$k.tab > $wdir/$i.$k.rpkm 
awk -v m=1 'NR==1{NF-=1;print $0;next} NR!=1{for(i=4;i<NF;i++) $i=(log(($i+m)/($NF+m))/log(2));NF-=1;print $0}' $wdir/$i.$k.rpkm > $wdir/$i.$k.m1.log2ratio
echo>&9; }&done;done;wait

#merge
for k in 1k.m1.log2ratio 1k.rpkm;do
cut -f 1,2,3 CC.$k | sorth -k1V,1 -k2n,2 -k3n,3 > bed.temp &
for i in ${NTF[@]}; do sorth -k1V,1 -k2n,2 -k3n,3 < $i.$k | cut -f 4 > $i.temp & done;wait
paste `for i in bed ${NTF[@]}; do echo -n $i.temp" ";done` > merge.$k
rm *temp
done 
#merge promoter
for k in promoter_10k-10k.m1.log2ratio promoter_10k-10k.rpkm;do
awk 'BEGIN{print "chr\tstart\tend\tgene"} NR==FNR&&FNR>1{a[$1"\t"$2]} NR>FNR&&($1"\t"$2 in a){print $1,$2,$3,$4}' CC.$k ~/ann/Uniq.mm10/all.gene/promoter_10k-10k.bed | sorth -k1V,1 -k2n,2 -k3n,3 -k4,4 > bed.temp &
for i in ${NTF[@]}; do sort -k1V,1 -k2n,2 -k3n,3 -k4,4 $i.$k | cut -f 4 > $i.temp & done;wait
paste `for i in bed ${NTF[@]}; do echo -n $i.temp" ";done` > merge.$k
wc -l *temp
rm *temp
done

