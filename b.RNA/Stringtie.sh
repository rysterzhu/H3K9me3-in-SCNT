
for i in $ddir/*sorted.bam; do read -u 9;{ 
k=$(basename $i .sorted.bam);k=${k/_RNA}
nohup stringtie $i -p 8 -G /home/share/gtf/mm10.gtf -l $k -A $wdir/1.expressions/$k.exp -C $wdir/2.gtfs/$k.gtf -o $wdir/3.denovo/$k.gtf > $wdir/0.logs/$k.log 
echo>&9;}&done&&wait-d "stringtie"
