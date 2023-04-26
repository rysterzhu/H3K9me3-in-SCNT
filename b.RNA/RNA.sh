ddir1=~/workspace/9.NT-ChIP/1.align/9.RNA_hisat2
ddir2=~/workspace/9.NT-ChIP/2.public/i.RNA-seq_wang/1.hisat2
wdir=~/workspace/9.NT-ChIP/b.RNA/9.DESeq2-ICM
nohup featureCounts -T 64 -p -B -C -M -t exon -g gene_id -a /home/share/gtf/mm10.gtf -o $wdir/TE.featureCount $ddir1/te*sorted.bam $ddir2/RWTE*sorted.bam >> $wdir/featureCount.log &

nohup featureCounts -T 64 -p -B -C -M -t exon -g gene_id -a /home/share/gtf/mm10.gtf -o $wdir/ICM2.featureCount $ddir1/icm*sorted.bam $ddir2/RWICM*sorted.bam >> $wdir/featureCount.log &



wdir=~/workspace/9.NT-ChIP/1.align/b.links-RNA-reps
ddir=~/workspace/9.NT-ChIP/1.align/9.RNA_hisat2/a.NT
NT3=("cc:CC" "6h:6h" "14h:14h" "e2cell:NT-e2cell" "l2cell:NT-l2cell" "4cell:NT-4cell" "8cell:NT-8cell" "morula:NT-Morula" "icm:NT-ICM" "te:NT-TE")
for i in ${NT3[@]}; do a=${i%%:*};b=${i##*:}
for j in $ddir/$a*sorted.bam; do k=$(basename $j .sorted.bam);k=${k/$a/$b};k=${k/_RNA}
ln -s $j $k.sorted.bam;done;done

ddir=~/workspace/9.NT-ChIP/2.public/c.Xiewei/5.RNA/4.hisat2
for i in $ddir/NF*sorted.bam; do k=$(basename $i .sorted.bam);k=${k/NF_/NF-}
ln -s $i $k.sorted.bam;done

wdir=~/workspace/9.NT-ChIP/1.align/b.links-RNA-reps/siMIT
ddir=~/workspace/9.NT-ChIP/1.align/9.RNA_hisat2/d.siMIT
for i in $ddir/*sorted.bam; do k=$(basename $i .sorted.bam);k=NF-${k/_RNA}
ln -s $i $k.sorted.bam;done
wdir=~/workspace/9.NT-ChIP/1.align/b.links-RNA-reps/siMIT
ddir=~/workspace/9.NT-ChIP/1.align/9.RNA_hisat2/e.siMIT-NT
for i in $ddir/*sorted.bam; do k=$(basename $i .sorted.bam);k=${k/_RNA}
ln -s $i $k.sorted.bam;done


ddir=~/workspace/9.NT-ChIP/2.public/i.RNA-seq_wang/1.hisat2
for i in $ddir/RW*sorted.bam; do k=$(basename $i .sorted.bam);k=${k/-/_rep};k=${k/RW/RW-}
ln -s $i $k.sorted.bam;done

mkd a.remove
mv NT-4cell_rep5* a.remove/
mv RW-* a.remove/
mv *rep1[1-9]* a.remove/

for i in TE_*; do mv $i ${i/TE_/TE-};done
for i in ICM_*; do mv $i ${i/ICM_/ICM-};done
for i in Morula_*; do mv $i ${i/Morula_/Morula-};done




awk 'NR==1{j=1;for(i=2;i<=NF;i++){if(index($i,"ontrol")>0){a[j]=i;j++}}} {out=$1;for(i=1;i<j;i++){out=out"\t"$(a[i])};print out}'

awk 'NR==1{j=1;for(i=2;i<=NF;i++){if(index($i,"ontrol")>0){a[j]=i;j++}}} {out=$1;for(i=1;i<j;i++){out=out"\t"$(a[i])};print out}' all.reps.TPM.txt |hd


awk 'NR==1{j=1;for(i=2;i<=NF;i++){if(match($i,"ICM|TE")>0){a[j]=i;j++}}} {out=$1;for(i=1;i<j;i++){out=out"\t"$(a[i])};print out}' addsiMIT.reps.TPM.txt | grep "Max|Mcrs1|Kat5|gene"


