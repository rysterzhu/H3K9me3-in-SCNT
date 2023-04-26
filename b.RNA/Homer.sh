ddir=~/workspace/9.NT-ChIP/1.align/9.RNA_hisat2/a.NT
wdir=~/workspace/9.NT-ChIP/b.RNA/7.homer-keepOne
mkd $wdir/logs $wdir/1.Tags
for i in $ddir/*sorted.bam; do k=$(basename $i .sorted.bam);k=NT-${k/RNA_}
(samtools view -@ 4 $i -o $wdir/$k.sam
nohup makeTagDirectory $wdir/1.Tags/$k.tag/ $wdir/$k.sam -format sam -keepOne -genome mm10 > $wdir/logs/$k.make.log) &
done
wait
ddir2=~/workspace/9.NT-ChIP/2.public/c.Xiewei/5.RNA/4.hisat2
for i in $ddir2/*sorted.bam; do k=$(basename $i .sorted.bam);k=${k/NF_/NF-}
(samtools view -@ 4 $i -o $wdir/$k.sam
nohup makeTagDirectory $wdir/1.Tags/$k.tag/ $wdir/$k.sam -format sam -keepOne -genome mm10 > $wdir/logs/$k.make.log)&
done
wait
cd $wdir/1.Tags
for i in L1 L2 L3; do 
nohup analyzeRepeats.pl repeats mm10 -d *tag -condense$i -tpm > ../all.Lall.condense$i.tpm.txt 2>>../logs/all.analyzeRepeats.$i.tpm.logs &
nohup analyzeRepeats.pl repeats mm10 -d *tag -condense$i -rpkm > ../all.Lall.condense$i.rpkm.txt 2>>../logs/all.analyzeRepeats.$i.rpkm.logs &
nohup analyzeRepeats.pl repeats mm10 -d *tag -condense$i > ../all.Lall.condense$i.txt 2>>../logs/all.analyzeRepeats.$i.default.logs &
nohup analyzeRepeats.pl repeats mm10 -d *tag -condense$i -rlog > ../all.Lall.condense$i.rlog.txt 2>>../logs/all.analyzeRepeats.$i.rlog.logs &
done
wait-m analyzeRepeats






###########################kdm4d  
ddir=~/workspace/9.NT-ChIP/2.public/i.RNA-seq_wang/1.hisat2
wdir=~/workspace/9.NT-ChIP/b.RNA/l.homer-wang
mkd $wdir/logs $wdir/1.Tags
for i in $ddir/R4BB4*sorted.bam $ddir/R4DB4*sorted.bam $ddir/R4B5BB4*sorted.bam $ddir/RW4*sorted.bam ; do k=$(basename $i .sorted.bam)
(samtools view -@ 8 $i -o $wdir/$k.sam
nohup makeTagDirectory $wdir/1.Tags/$k.tag/ $wdir/$k.sam -format sam -keepOne -genome mm10 > $wdir/logs/$k.make.log) &
done
wait
cd $wdir/1.Tags
i=L1
nohup analyzeRepeats.pl repeats mm10 -d *tag -condense$i -rpkm > ../all.Lall.condense$i.rpkm.txt 2>>../logs/all.analyzeRepeats.$i.rpkm.logs &
nohup analyzeRepeats.pl repeats mm10 -d *tag -L1 MERVL-int -rpkm > ../MERVL.rpkm.txt 2>../logs/MERVL.analyzeRepeats.rpkm.log &
nohup analyzeRepeats.pl repeats mm10 -d *tag -L1 MT2_Mm -rpkm > ../MT2.rpkm.txt 2>../logs/MT2.analyzeRepeats.rpkm.log &

ddir=~/workspace/9.NT-ChIP/2.public/f.Zhangyi/7.hisat2
wdir=~/workspace/9.NT-ChIP/b.RNA/m.homer-zhangyi
mkd $wdir/logs $wdir/1.Tags
for i in $ddir/*bam ; do k=$(basename $i .bam)
(samtools view -@ 8 $i -o $wdir/$k.sam
nohup makeTagDirectory $wdir/1.Tags/$k.tag/ $wdir/$k.sam -format sam -keepOne -genome mm10 > $wdir/logs/$k.make.log) &
done
wait
cd $wdir/1.Tags
i=L1
nohup analyzeRepeats.pl repeats mm10 -d *tag -condense$i -rpkm > ../all.Lall.condense$i.rpkm.txt 2>>../logs/all.analyzeRepeats.$i.rpkm.logs &
nohup analyzeRepeats.pl repeats mm10 -d *tag -L1 MERVL-int -rpkm > ../MERVL.rpkm.txt 2>../logs/MERVL.analyzeRepeats.rpkm.log &
nohup analyzeRepeats.pl repeats mm10 -d *tag -L1 MT2_Mm -rpkm > ../MT2.rpkm.txt 2>../logs/MT2.analyzeRepeats.rpkm.log &
