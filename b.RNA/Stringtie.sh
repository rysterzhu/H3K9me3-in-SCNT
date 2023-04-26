
ddir=~/workspace/9.NT-ChIP/1.align/9.RNA_hisat2/[df]*
wdir=~/workspace/9.NT-ChIP/b.RNA/i.stringtie-siMIT
mkd $wdir/1.expressions $wdir/2.gtfs $wdir/3.denovo $wdir/0.logs
mult 8
for i in $ddir/*sorted.bam; do read -u 9;{ 
k=$(basename $i .sorted.bam);k=${k/_RNA}
nohup stringtie $i -p 8 -G /home/share/gtf/mm10.gtf -l $k -A $wdir/1.expressions/$k.exp -C $wdir/2.gtfs/$k.gtf -o $wdir/3.denovo/$k.gtf > $wdir/0.logs/$k.log 
echo>&9;}&done&&wait-d "stringtie"

awk 'FNR==1{split(FILENAME,temp,"[./_]");key=temp[3]"_"temp[4]"_"temp[5];next} $1!~/_/{print $1,$9,key}' j.uploads/*exp | datamash -s crosstab 1,3 mean 2 |awk 'NR==1{print "gene"$0} NR>1'> all-uploads.reps.TPM.txt

#####Mcrs1 
awk 'FNR==1{split(FILENAME,temp,"[./_]");key=temp[3]"_"temp[4]"_"temp[5];next} $1!~/_/{print $1,$9,key}' 1.expressions/NT-*_old6**_rep[123456].exp 1.expressions/NT-*_oeMcrs1*_rep[1-6].exp 1.expressions/NF-*_mphControl_rep*exp | datamash -s crosstab 1,3 mean 2 |awk 'NR==1{print "gene"$0} NR>1'> blastocyst-oeMcrs1.reps.TPM.txt


