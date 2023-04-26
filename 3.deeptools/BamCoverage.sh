ddir=~/workspace/9.NT-ChIP/1.align/f.K9me3_keep_lowqual/1.merge-reps
ddir=~/workspace/9.NT-ChIP/1.align/g.K9me3_kdm4b/c.merge_rep
#ddir=~/workspace/9.NT-ChIP/1.align/h.K9me3_oeMcrs1/c.merge_rep
ddir=~/workspace/9.NT-ChIP/1.align/i.K9me3_sertoli/c.merge_rep
wdir=~/workspace/9.NT-ChIP/3.deeptools/1.bamCoverage/c.K9me3_keep_lowqual
mkd $wdir/logs
for i in $ddir/Sertoli*redun.bam; do k=$(basename $i .redun.bam)
nohup bamCoverage -b $i -o $wdir/$k.bw --normalizeUsing RPKM -p 8 --ignoreDuplicates -bs 50 --smoothLength 150 --samFlagInclude 2 > $wdir/logs/$k.log &
done 


ddir=~/workspace/9.NT-ChIP/1.align/f.K9me3_keep_lowqual/0.links-raw-reps
ddir=~/workspace/9.NT-ChIP/1.align/g.K9me3_kdm4b/a.raw
#ddir=~/workspace/9.NT-ChIP/1.align/h.K9me3_oeMcrs1/a.raw
ddir=~/workspace/9.NT-ChIP/1.align/i.K9me3_sertoli/a.raw
wdir=~/workspace/9.NT-ChIP/3.deeptools/1.bamCoverage/d.K9me3_keep_lowqual-reps
mkd $wdir/logs
for i in $ddir/Sertoli*sorted.bam; do k=$(basename $i .sorted.bam)
nohup bamCoverage -b $i -o $wdir/$k.bw --normalizeUsing RPKM -p 8 --ignoreDuplicates -bs 50 --smoothLength 150 --samFlagInclude 2 --samFlagExclude 12 > $wdir/logs/$k.log &
done 


