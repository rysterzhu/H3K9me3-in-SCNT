for i in $ddir/*redun.bam; do
k=$(basename $i .redun.bam)
nohup bamCoverage -b $i -o $wdir/$k.bw --normalizeUsing RPKM -p 8 -bs 50 --smoothLength 150 --ignoreDuplicates --samFlagInclude 2 --samFlagExclude 1036 --maxFragmentLength 1000 > $wdir/logs/$k.log &
done



