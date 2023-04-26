

ddir=~/workspace/9.NT-ChIP/1.align/f.K9me3_keep_lowqual/1.merge-reps
ddir=~/workspace/9.NT-ChIP/1.align/g.K9me3_kdm4b/c.merge_rep
ddir=~/workspace/9.NT-ChIP/1.align/h.K9me3_oeMcrs1/c.merge_rep
ddir=~/workspace/9.NT-ChIP/1.align/i.K9me3_sertoli/c.merge_rep
wdir=~/workspace/9.NT-ChIP/3.deeptools/5.compareByInput/a.keep-lowqual
mkd $wdir/logs
for i in $ddir/*K9me3.redun.bam; do
o=$(basename $i .redun.bam);k=${o%_*};
nohup bamCompare -b1 $i -b2 $ddir/${k}_input.redun.bam -o $wdir/${o}.bw --centerReads --operation log2 --scaleFactorsMethod None --pseudocount 1 -bs 50 --normalizeUsing RPKM --smoothLength 150 -p 32 --samFlagInclude 2 >$wdir/logs/$o.log &
done 

#ddir=~/workspace/9.NT-ChIP/1.align/f.K9me3_keep_lowqual/0.links-raw-reps
#ddir=~/workspace/9.NT-ChIP/1.align/g.K9me3_kdm4b/b.filtered
#wdir=~/workspace/9.NT-ChIP/3.deeptools/5.compareByInput/b.keep-lowqual-reps
#mkd $wdir/logs
#for i in $ddir/*K9me3_rep*sorted.bam; do o=$(basename $i .sorted.bam);k=${o%_*};
#nohup bamCompare -b1 $i -b2 $ddir/${k}_input_rep1.sorted.bam -o $wdir/${o}.bw --centerReads --operation log2 --scaleFactorsMethod None --pseudocount 1 -bs 50 --normalizeUsing RPKM --smoothLength 150 -p 2 --samFlagInclude 2 --samFlagExclude 12 >$wdir/logs/$o.log &
#done 






