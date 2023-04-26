wdir=~/workspace/9.NT-ChIP/d.cut_run/2.bamCoverage/
ddir=~/workspace/9.NT-ChIP/d.cut_run/0.align/b.filtered
mkd $wdir/logs;mult 8;
for i in $ddir/N*sorted.bam; do read -u 9;{
k=$(basename $i .sorted.bam)
nohup bamCoverage -b $i -o $wdir/$k.bw --normalizeUsing RPKM -p 8 -bs 50 --smoothLength 150 --ignoreDuplicates --samFlagInclude 2 --samFlagExclude 1036 --maxFragmentLength 1000 > $wdir/logs/$k.log 
echo>&9;}&done&&wait

wdir=~/workspace/9.NT-ChIP/d.cut_run/3.bamCoverage-merge
ddir=~/workspace/9.NT-ChIP/d.cut_run/0.align/c.merge_rep
mkd $wdir/logs;mult 8;
for i in $ddir/N*redun.bam; do read -u 9;{
k=$(basename $i .redun.bam)
nohup bamCoverage -b $i -o $wdir/$k.bw --normalizeUsing RPKM -p 8 -bs 50 --smoothLength 150 --ignoreDuplicates --samFlagInclude 2 --samFlagExclude 1036 --maxFragmentLength 1000 > $wdir/logs/$k.log 
echo>&9;}&done&&wait

ddir=~/workspace/9.NT-ChIP/1.align/d.sih2-K9me3/1.mm10/b.filtered
wdir=~/workspace/9.NT-ChIP/d.cut_run/2.bamCoverage/
mkd $wdir/logs;
for i in $ddir/*rep[34].sorted.bam; do 
k=$(basename $i .sorted.bam)
nohup bamCoverage -b $i -o $wdir/$k.bw --normalizeUsing RPKM -p 8 -bs 50 --smoothLength 150 --ignoreDuplicates --samFlagInclude 2 --samFlagExclude 1036 --maxFragmentLength 1000 > $wdir/logs/$k.log &
done

wdir=~/workspace/9.NT-ChIP/d.cut_run/3.bamCoverage-merge
ddir=~/workspace/9.NT-ChIP/1.align/d.sih2-K9me3/1.mm10/c.merge_rep
mkd $wdir/logs;
for i in $ddir/*redun.bam; do
k=$(basename $i .redun.bam)
nohup bamCoverage -b $i -o $wdir/$k.bw --normalizeUsing RPKM -p 8 -bs 50 --smoothLength 150 --ignoreDuplicates --samFlagInclude 2 --samFlagExclude 1036 --maxFragmentLength 1000 > $wdir/logs/$k.log &
done



wdir=~/workspace/9.NT-ChIP/d.cut_run/3.bamCoverage-merge/calibration
ddir=~/workspace/9.NT-ChIP/1.align/d.sih2-K9me3/1.mm10/b.filtered
mkd $wdir/logs;
for i in $ddir/*rep2.sorted.bam; do 
k=$(basename $i .sorted.bam)
sfactor=`awk -v scale=1e6 'NR==1{print scale/$1}' ~/workspace/9.NT-ChIP/1.align/d.sih2-K9me3/2.spiking/b.filtered/logs/$k.txt`
nohup bamCoverage -b $i -o $wdir/$k.bw --scaleFactor $sfactor --normalizeUsing None -p 8 -bs 50 --smoothLength 150 --ignoreDuplicates --samFlagInclude 2 --samFlagExclude 1036 --maxFragmentLength 1000 > $wdir/logs/$k.log  &
done




