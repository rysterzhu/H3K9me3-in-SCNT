wdir=~/workspace/9.NT-ChIP/4.callPeaks/3.homer
ddir=~/workspace/9.NT-ChIP/1.align/f.K9me3_keep_lowqual/1.merge-reps
mkd $wdir/logs
for i in $ddir/*.redun.bam; do k=NT-$(basename $i .redun.bam)
(samtools view -@ 8 -h $i -o $wdir/$k.sam
nohup makeTagDirectory $wdir/$k.tag/ $wdir/$k.sam -format sam -genome mm10 > $wdir/logs/$k.make.log
rm $wdir/$k.sam)&
done
wait

nohup analyzeRepeats.pl repeats mm10 -d *tag -condenseL2 > ../all.Lall.condenseL2.txt 2>../logs/all.condenseL2.logs &
nohup analyzeRepeats.pl repeats mm10 -d *tag -condenseL1 > ../all.Lall.condenseL1.txt 2>../logs/all.condenseL1.logs &
nohup analyzeRepeats.pl repeats mm10 -d *tag -condenseL3 > ../all.Lall.condenseL3.txt 2>../logs/all.condenseL3.logs &
nohup analyzeRepeats.pl repeats mm10 -d *tag -condenseL2 -rlog > ../all.Lall.condenseL2-rlog.txt 2>../logs/all.condenseL2-rlog.logs &
nohup analyzeRepeats.pl repeats mm10 -d *tag -condenseL1 -rlog > ../all.Lall.condenseL1-rlog.txt 2>../logs/all.condenseL1-rlog.logs &
nohup analyzeRepeats.pl repeats mm10 -d *tag -condenseL3 -rlog > ../all.Lall.condenseL3-rlog.txt 2>../logs/all.condenseL3-rlog.logs &
wait-m 

for i in L1 L2 L3; do 
nohup analyzeRepeats.pl repeats mm10 -d *tag -condense$i -rpkm > ../all.Lall.condense$i.rpkm.txt 2>>../logs/all.analyzeRepeats.$i.rpkm.logs &
done



nohup analyzeRepeats.pl repeats mm10 -d *tag -L1 MERVL-int -rpkm > ../MERVL.rpkm.txt 2>../logs/MERVL.analyzeRepeats.rpkm.log &
nohup analyzeRepeats.pl repeats mm10 -d *tag -L3 LTR -rpkm > ../LTR.rpkm.txt 2>../logs/LTR.analyzeRepeats.rpkm.log &
nohup analyzeRepeats.pl repeats mm10 -d *tag -L3 ERVL -rpkm > ../ERVL.rpkm.txt 2>../logs/ERVL.analyzeRepeats.rpkm.log &

################2021年10月14日 kdm4b
wdir=~/workspace/9.NT-ChIP/4.callPeaks/d.homer-reps_kdm4b
ddir=~/workspace/9.NT-ChIP/1.align/g.K9me3_kdm4b/b.filtered
mkd $wdir/logs $wdir/1.Tags
for i in $ddir/*.sorted.bam; do
k=$(basename $i );k=${k%%.*};k=${k/_input_rep1/_input};
(samtools view -@ 8 -F 1036 -h $i -o $wdir/$k.sam
nohup makeTagDirectory $wdir/1.Tags/$k.tag/ $wdir/$k.sam -keepOne -format sam -genome mm10 > $wdir/logs/$k.make.log) &
done
cd $wdir/1.Tags
for i in L1 L2 L3; do 
nohup analyzeRepeats.pl repeats mm10 -d *tag -condense$i -rpkm > ../all.Lall.condense$i.rpkm.txt 2>>../logs/all.analyzeRepeats.$i.rpkm.logs &
nohup analyzeRepeats.pl repeats mm10 -d *tag -condense$i > ../all.Lall.condense$i.txt 2>../logs/all.condense$i.logs &
nohup analyzeRepeats.pl repeats mm10 -d *tag -condense$i -rlog > ../all.Lall.condense$i.rlog.txt 2>../logs/all.condense$i.rlog.logs &
done
nohup analyzeRepeats.pl repeats mm10 -d *tag -L1 MERVL-int -rpkm > ../MERVL.rpkm.txt 2>../logs/MERVL.analyzeRepeats.rpkm.log &
nohup analyzeRepeats.pl repeats mm10 -d *tag -L3 LTR -rpkm > ../LTR.rpkm.txt 2>../logs/LTR.analyzeRepeats.rpkm.log &
nohup analyzeRepeats.pl repeats mm10 -d *tag -L3 ERVL -rpkm > ../ERVL.rpkm.txt 2>../logs/ERVL.analyzeRepeats.rpkm.log &


######################2022年4月16日 oeMcrs1
wdir=~/workspace/9.NT-ChIP/4.callPeaks/i.homer-reps_oeMcrs1
ddir=~/workspace/9.NT-ChIP/1.align/h.K9me3_oeMcrs1/b.filtered
ddir2=~/workspace/9.NT-ChIP/1.align/h.K9me3_oeMcrs1/c.merge_rep
mkd $wdir/logs $wdir/1.Tags
for i in $ddir/*K9me3*.sorted.bam $ddir2/*input.redun.bam; do
k=$(basename $i );k=${k%%.*};
(samtools view -@ 8 -F 1036 -h $i -o $wdir/$k.sam
nohup makeTagDirectory $wdir/1.Tags/$k.tag/ $wdir/$k.sam -keepOne -format sam -genome mm10 > $wdir/logs/$k.make.log) &
done
cd $wdir/1.Tags
ln -s ~/workspace/9.NT-ChIP/4.callPeaks/c.homer-reps/1.Tags/NT-[IT]*tag ./
for i in L1 L2 L3; do 
nohup analyzeRepeats.pl repeats mm10 -d *tag -condense$i -rpkm > ../all.Lall.condense$i.rpkm.txt 2>>../logs/all.analyzeRepeats.$i.rpkm.logs &
nohup analyzeRepeats.pl repeats mm10 -d *tag -condense$i > ../all.Lall.condense$i.txt 2>../logs/all.condense$i.logs &
nohup analyzeRepeats.pl repeats mm10 -d *tag -condense$i -rlog > ../all.Lall.condense$i.rlog.txt 2>../logs/all.condense$i.rlog.logs &
done
nohup analyzeRepeats.pl repeats mm10 -d *tag -L1 MERVL-int -rpkm > ../MERVL.rpkm.txt 2>../logs/MERVL.analyzeRepeats.rpkm.log &
nohup analyzeRepeats.pl repeats mm10 -d *tag -L3 LTR -rpkm > ../LTR.rpkm.txt 2>../logs/LTR.analyzeRepeats.rpkm.log &
nohup analyzeRepeats.pl repeats mm10 -d *tag -L3 ERVL -rpkm > ../ERVL.rpkm.txt 2>../logs/ERVL.analyzeRepeats.rpkm.log &

nohup analyzeRepeats.pl repeats mm10 -d *tag -L2 LINE -rpkm > ../LINE.rpkm.txt 2>../logs/LINE.analyzeRepeats.rpkm.log &
##merge rep
wdir=~/workspace/9.NT-ChIP/4.callPeaks/j.homer-mergeRep_oeMcrs1
ddir=~/workspace/9.NT-ChIP/1.align/h.K9me3_oeMcrs1/c.merge_rep
mkd $wdir/logs $wdir/1.Tags
for i in $ddir/*.redun.bam; do
k=$(basename $i );k=${k%%.*};
(samtools view -@ 8 -F 1036 -h $i -o $wdir/$k.sam
nohup makeTagDirectory $wdir/1.Tags/$k.tag/ $wdir/$k.sam -keepOne -format sam -genome mm10 > $wdir/logs/$k.make.log) &
done






