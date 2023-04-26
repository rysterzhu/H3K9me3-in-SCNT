#################################################################################################################
###########################################################2018-03-01 ChIP H3K4me3 H3K27me3 H3K9me3
ddir1=/home/share/DATA/clean_reads/GaoSR_mouse_nuclear_transfer_ChIP_20180301/20180301/clean_data_20180115
ddir2=/home/share/DATA/clean_reads/GaoSR_mouse_nuclear_transfer_ChIP_20180301/20180301/4cell_raw_data_20171023
ddir3=/home/share/DATA/clean_reads/GaoSR_mouse_nuclear_transfer_ChIP_20180301/2cell_20180302/Clean
cd ~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links

ln -s $ddir1/LC-0115-10_L4_A010.R1.clean.fastq.gz 8cell_K9me3.R1.fastq.gz
ln -s $ddir1/LC-0115-10_L4_A010.R2.clean.fastq.gz 8cell_K9me3.R2.fastq.gz 
ln -s $ddir1/LC-0115-11_L4_A011.R1.clean.fastq.gz 8cell_K4me3.R1.fastq.gz
ln -s $ddir1/LC-0115-11_L4_A011.R2.clean.fastq.gz 8cell_K4me3.R2.fastq.gz
ln -s $ddir1/LC-0115-12_L4_A012.R1.clean.fastq.gz 8cell_K27me3.R1.fastq.gz
ln -s $ddir1/LC-0115-12_L4_A012.R2.clean.fastq.gz 8cell_K27me3.R2.fastq.gz
ln -s $ddir1/LC-0115-1_L3_A001.R1.clean.fastq.gz 6h_input.R1.fastq.gz
ln -s $ddir1/LC-0115-1_L3_A001.R2.clean.fastq.gz 6h_input.R2.fastq.gz
ln -s $ddir1/LC-0115-2_L3_A002.R1.clean.fastq.gz 6h_K9me3.R1.fastq.gz
ln -s $ddir1/LC-0115-2_L3_A002.R2.clean.fastq.gz 6h_K9me3.R2.fastq.gz
ln -s $ddir1/LC-0115-4_L3_A004.R1.clean.fastq.gz 14h_input.R1.fastq.gz
ln -s $ddir1/LC-0115-4_L3_A004.R2.clean.fastq.gz 14h_input.R2.fastq.gz
ln -s $ddir1/LC-0115-6_L1_A006.R1.clean.fastq.gz 14h_K4me3.R1.fastq.gz
ln -s $ddir1/LC-0115-6_L1_A006.R2.clean.fastq.gz 14h_K4me3.R2.fastq.gz
ln -s $ddir1/LC-0115-7_L4_A007.R1.clean.fastq.gz 4cell_K4me3.R1.fastq.gz
ln -s $ddir1/LC-0115-7_L4_A007.R2.clean.fastq.gz 4cell_K4me3.R2.fastq.gz
ln -s $ddir1/LC-0115-8_L4_A008.R1.clean.fastq.gz 4cell_K27me3.R1.fastq.gz
ln -s $ddir1/LC-0115-8_L4_A008.R2.clean.fastq.gz 4cell_K27me3.R2.fastq.gz
ln -s $ddir1/LC-0115-9_L1_A009.R1.clean.fastq.gz 8cell_input.R1.fastq.gz
ln -s $ddir1/LC-0115-9_L1_A009.R2.clean.fastq.gz 8cell_input.R2.fastq.gz

ln -s $ddir2/XRM-1023-1_H72JKCCXY_L1_1.fq.gz 4cell_K9me3.R1.fastq.gz
ln -s $ddir2/XRM-1023-1_H72JKCCXY_L1_2.fq.gz 4cell_K9me3.R2.fastq.gz
ln -s $ddir2/XRM-1023-2_H72JKCCXY_L1_1.fq.gz 4cell_input.R1.fastq.gz
ln -s $ddir2/XRM-1023-2_H72JKCCXY_L1_2.fq.gz 4cell_input.R2.fastq.gz

ln -s $ddir3/LC-0925-1_R1.fq.gz 2cell_K9me3.R1.fastq.gz
ln -s $ddir3/LC-0925-1_R2.fq.gz 2cell_K9me3.R2.fastq.gz
ln -s $ddir3/LC-0925-2_R1.fq.gz 2cell_input.R1.fastq.gz
ln -s $ddir3/LC-0925-2_R2.fq.gz 2cell_input.R2.fastq.gz
ln -s $ddir3/LC-0925-3_R1.fq.gz 2cell_K4me3.R1.fastq.gz
ln -s $ddir3/LC-0925-3_R2.fq.gz 2cell_K4me3.R2.fastq.gz
ln -s $ddir3/LC-0925-4_R1.fq.gz 2cell_K27me3.R1.fastq.gz
ln -s $ddir3/LC-0925-4_R2.fq.gz 2cell_K27me3.R2.fastq.gz


nohup fastqc -o 1.qc/1.20180301/ 0.links/*gz -t 12 &

ddir=~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links
wdir=~/workspace/9.NT-ChIP/0.cutadapt-qc/2.fastp
for i in $ddir/*R1.fastq.gz; do o=${i##*/}
nohup fastp -i $i -I ${i/R1/R2} -o $wdir/$o -O $wdir/${o/R1/R2} -p -w 8 -l 75 \
-h $wdir/${o/R1.fastq.gz/fastp.html} -j $wdir/${o/R1.fastq.gz/fastp.json} > $wdir/${o/R1.fastq.gz/log} &
done

mkdir 1.qc
nohup fastqc -o 1.qc/ *gz -t 12 &

##############fastp去掉的adapter不干净 (如果样本的adapter较多，fastp不合适)
wdir=~/workspace/9.NT-ChIP/0.cutadapt-qc/3.cutadapt
ddir=~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links
for i in $ddir/2cell*R1.fastq.gz; do o=$wdir/${i##*/}
nohup cutadapt --trim-n -q 25,25 -m 75 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o $o -p ${o/R1/R2} $i ${i/R1/R2} > ${o/R1.fastq.gz/log} &
done




#################################################################################################################
###########################################################2018-03-20 ChIP H3K4me3 H3K27me3 H3K9me3
ddir=/home/share/DATA/clean_reads/GaoSR_mouse_nuclear_transfer_ChIP_20180301/20180320
ln -s $ddir/XRM-0202-11_L2_A022.R1.clean.fastq.gz morula_input.R1.fastq.gz
ln -s $ddir/XRM-0202-11_L2_A022.R2.clean.fastq.gz morula_input.R2.fastq.gz
ln -s $ddir/XRM-0202-12_L2_A023.R1.clean.fastq.gz morula_K4me3.R1.fastq.gz
ln -s $ddir/XRM-0202-12_L2_A023.R2.clean.fastq.gz morula_K4me3.R2.fastq.gz
ln -s $ddir/XRM-0202-14_L2_A026.R1.clean.fastq.gz morula_K27me3.R1.fastq.gz
ln -s $ddir/XRM-0202-14_L2_A026.R2.clean.fastq.gz morula_K27me3.R2.fastq.gz
ln -s $ddir/XRM-0202-5_L1_A010.R1.clean.fastq.gz 8cell_input_rep2.R1.fastq.gz
ln -s $ddir/XRM-0202-5_L1_A010.R2.clean.fastq.gz 8cell_input_rep2.R2.fastq.gz
ln -s $ddir/XRM-0202-6_L1_A011.R1.clean.fastq.gz 8cell_K4me3_rep2.R1.fastq.gz
ln -s $ddir/XRM-0202-6_L1_A011.R2.clean.fastq.gz 8cell_K4me3_rep2.R2.fastq.gz
ln -s $ddir/XRM-0202-7_L1_A012.R1.clean.fastq.gz 8cell_K9me3_rep2.R1.fastq.gz
ln -s $ddir/XRM-0202-7_L1_A012.R2.clean.fastq.gz 8cell_K9me3_rep2.R2.fastq.gz
ln -s $ddir/XRM-0202-8_L1_A021.R1.clean.fastq.gz 8cell_K9me3_rep3.R1.fastq.gz
ln -s $ddir/XRM-0202-8_L1_A021.R2.clean.fastq.gz 8cell_K9me3_rep3.R2.fastq.gz
ln -s $ddir/XRM-0202-9_L1_A003.R1.clean.fastq.gz 6h_K4me3.R1.fastq.gz
ln -s $ddir/XRM-0202-9_L1_A003.R2.clean.fastq.gz 6h_K4me3.R2.fastq.gz

nohup fastqc -o 1.qc/b.20180320/ 0.links/b.20180320/*gz -t 36 &


wdir=~/workspace/9.NT-ChIP/0.cutadapt-qc/3.cutadapt/b.20180320
ddir=~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links/b.20180320
for i in $ddir/*R1.fastq.gz; do o=$wdir/${i##*/}
nohup cutadapt --trim-n -q 25,25 -m 75 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o $o -p ${o/R1/R2} $i ${i/R1/R2} > ${o/R1.fastq.gz/log} &
done

wait 
nohup fastqc -o qc/ *gz -t 12 &


#######################################################################################################################
##############################2018-5-28 6h 14h 2cell K9me3
ddir=/home/share/DATA/clean_reads/GaoSR_mouse_nuclear_transfer_ChIP_20180301/20180528
wdir=~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links/c.20180528

for i in $ddir/*/*_1.fq.gz; do key=${i%/*};key=${key##*/}
ln -s $i $wdir/$key.R1.fastq.gz
done
for i in $ddir/*/*_2.fq.gz; do key=${i%/*};key=${key##*/}
ln -s $i $wdir/$key.R2.fastq.gz
done


nohup fastqc -o 1.qc/c.20180528 0.links/c.20180528/*gz -t 12 &

wdir=~/workspace/9.NT-ChIP/0.cutadapt-qc/3.cutadapt/c.20180528
ddir=~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links/c.20180528
mkdir -p $wdir/qc
for i in $ddir/*R1.fastq.gz; do o=$wdir/${i##*/}
nohup cutadapt --trim-n -q 25,25 -m 75 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o $o -p ${o/R1/R2} $i ${i/R1/R2} > ${o/R1.fastq.gz/log} &
done
nohup fastqc -o $wdir/qc $wdir/*gz -t 12 &

#######################################################################################################################
##############################2018-6-25 4cell 8cell morula + 加测的数据
ddir=/home/share/DATA/clean_reads/GaoSR_mouse_nuclear_transfer_ChIP_20180301/20180625
wdir=~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links/d.20180626
mkdir -p $wdir ${wdir/0.links/1.qc} ${wdir/0.links/3.cutadapt} ${wdir/0.links/4.qc_cutadapt}
#加测的数据
ln -s $ddir/s0427-xrm-1_TKD180402516_1.fq.gz $wdir/2cell_input.R1.fastq.gz
ln -s $ddir/s0427-xrm-1_TKD180402516_2.fq.gz $wdir/2cell_input.R2.fastq.gz
ln -s $ddir/s0507-XRM-1_TKD180500666_1.fq.gz $wdir/2cell_K9me3.R1.fastq.gz
ln -s $ddir/s0507-XRM-1_TKD180500666_2.fq.gz $wdir/2cell_K9me3.R2.fastq.gz
ln -s $ddir/s0507-XRM-2_TKD180500675_1.fq.gz $wdir/6h_input.R1.fastq.gz
ln -s $ddir/s0507-XRM-2_TKD180500675_2.fq.gz $wdir/6h_input.R2.fastq.gz
ln -s $ddir/s0507-XRM-3_TKD180500676_1.fq.gz $wdir/6h_K9me3.R1.fastq.gz
ln -s $ddir/s0507-XRM-3_TKD180500676_2.fq.gz $wdir/6h_K9me3.R2.fastq.gz
ln -s $ddir/s0510-XRM-1_TKD180500933_1.fq.gz $wdir/14h_input.R1.fastq.gz
ln -s $ddir/s0510-XRM-1_TKD180500933_2.fq.gz $wdir/14h_input.R2.fastq.gz
ln -s $ddir/s0510-XRM-2_TKD180500934_1.fq.gz $wdir/14h_K9me3.R1.fastq.gz
ln -s $ddir/s0510-XRM-2_TKD180500934_2.fq.gz $wdir/14h_K9me3.R2.fastq.gz
#后上机的数据
ln -s $ddir/s0530-XRM-1_TKD180600339_1.fq.gz $wdir/4cell_input.R1.fastq.gz
ln -s $ddir/s0530-XRM-1_TKD180600339_2.fq.gz $wdir/4cell_input.R2.fastq.gz
ln -s $ddir/s0530-XRM-2_TKD180600340_1.fq.gz $wdir/4cell_K9me3.R1.fastq.gz
ln -s $ddir/s0530-XRM-2_TKD180600340_2.fq.gz $wdir/4cell_K9me3.R2.fastq.gz
ln -s $ddir/s0530-XRM-3_TKD180600341_1.fq.gz $wdir/8cell_input.R1.fastq.gz
ln -s $ddir/s0530-XRM-3_TKD180600341_2.fq.gz $wdir/8cell_input.R2.fastq.gz
ln -s $ddir/s0530-XRM-4_TKD180600342_1.fq.gz $wdir/8cell_K9me3.R1.fastq.gz
ln -s $ddir/s0530-XRM-4_TKD180600342_2.fq.gz $wdir/8cell_K9me3.R2.fastq.gz
ln -s $ddir/s0530-XRM-5_TKD180600343_1.fq.gz $wdir/morula_input.R1.fastq.gz
ln -s $ddir/s0530-XRM-5_TKD180600343_2.fq.gz $wdir/morula_input.R2.fastq.gz
ln -s $ddir/s0530-XRM-6_TKD180600344_1.fq.gz $wdir/morula_K9me3_rep1.R1.fastq.gz
ln -s $ddir/s0530-XRM-6_TKD180600344_2.fq.gz $wdir/morula_K9me3_rep1.R2.fastq.gz
ln -s $ddir/s0530-XRM-7_TKD180600345_1.fq.gz $wdir/morula_K9me3_rep2.R1.fastq.gz
ln -s $ddir/s0530-XRM-7_TKD180600345_2.fq.gz $wdir/morula_K9me3_rep2.R2.fastq.gz

nohup fastqc -o ${wdir/0.links/1.qc} $wdir/*gz -t 20 &
py35
for i in $wdir/*R1.fastq.gz; do o=${i/0.links/3.cutadapt}
nohup cutadapt -j 6 --trim-n -q 25,25 -m 75 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o $o -p ${o/R1/R2} $i ${i/R1/R2} > ${o/R1.fastq.gz/log} &
done
wait

nohup fastqc -o ${wdir/0.links/4.qc_cutadapt} ${wdir/0.links/3.cutadapt}/*.gz -t 20 &

###########################################################################################
##2018年7月8日 
ddir=/home/share/DATA/clean_reads/GaoSR_mouse_nuclear_transfer_ChIP_20180301/20180707/1.rawdata
wdir=~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links/e.20180708
mkdir -p $wdir ${wdir/0.links/1.qc} ${wdir/0.links/3.cutadapt} ${wdir/0.links/4.qc_cutadapt}

ln -s $ddir/s0612-XRM-1_TKD180601211_1.fq.gz $wdir/2cell_input.R1.fastq.gz
ln -s $ddir/s0612-XRM-2_TKD180601212_1.fq.gz $wdir/2cell_K4me3.R1.fastq.gz
ln -s $ddir/s0612-XRM-3_TKD180601213_1.fq.gz $wdir/4cell_input.R1.fastq.gz
ln -s $ddir/s0612-XRM-4_TKD180601214_1.fq.gz $wdir/4cell_K4me3.R1.fastq.gz
ln -s $ddir/s0612-XRM-5_TKD180601215_1.fq.gz $wdir/8cell_input.R1.fastq.gz
ln -s $ddir/s0612-XRM-6_TKD180601216_1.fq.gz $wdir/8cell_K4me3.R1.fastq.gz
ln -s $ddir/s0612-XRM-7_TKD180601217_1.fq.gz $wdir/icm_input.R1.fastq.gz
ln -s $ddir/s0612-XRM-8_TKD180601218_1.fq.gz $wdir/icm_K9me3.R1.fastq.gz
ln -s $ddir/s0612-XRM-9_TKD180601219_1.fq.gz $wdir/te_input.R1.fastq.gz
ln -s $ddir/s0612-XRM-10_TKD180601220_1.fq.gz $wdir/te_K9me3.R1.fastq.gz
ln -s $ddir/s0612-XRM-11_TKD180601221_1.fq.gz $wdir/blast_input.R1.fastq.gz
ln -s $ddir/s0612-XRM-12_TKD180601222_1.fq.gz $wdir/blast_K9me3.R1.fastq.gz

ln -s $ddir/s0612-XRM-1_TKD180601211_2.fq.gz $wdir/2cell_input.R2.fastq.gz
ln -s $ddir/s0612-XRM-2_TKD180601212_2.fq.gz $wdir/2cell_K4me3.R2.fastq.gz
ln -s $ddir/s0612-XRM-3_TKD180601213_2.fq.gz $wdir/4cell_input.R2.fastq.gz
ln -s $ddir/s0612-XRM-4_TKD180601214_2.fq.gz $wdir/4cell_K4me3.R2.fastq.gz
ln -s $ddir/s0612-XRM-5_TKD180601215_2.fq.gz $wdir/8cell_input.R2.fastq.gz
ln -s $ddir/s0612-XRM-6_TKD180601216_2.fq.gz $wdir/8cell_K4me3.R2.fastq.gz
ln -s $ddir/s0612-XRM-7_TKD180601217_2.fq.gz $wdir/icm_input.R2.fastq.gz
ln -s $ddir/s0612-XRM-8_TKD180601218_2.fq.gz $wdir/icm_K9me3.R2.fastq.gz
ln -s $ddir/s0612-XRM-9_TKD180601219_2.fq.gz $wdir/te_input.R2.fastq.gz
ln -s $ddir/s0612-XRM-10_TKD180601220_2.fq.gz $wdir/te_K9me3.R2.fastq.gz
ln -s $ddir/s0612-XRM-11_TKD180601221_2.fq.gz $wdir/blast_input.R2.fastq.gz
ln -s $ddir/s0612-XRM-12_TKD180601222_2.fq.gz $wdir/blast_K9me3.R2.fastq.gz

nohup fastqc -o ${wdir/0.links/1.qc} $wdir/*gz -t 8 &
py35
for i in $wdir/*R1.fastq.gz; do o=${i/0.links/3.cutadapt}
nohup cutadapt -j 4 --trim-n -q 25,25 -m 75 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o $o -p ${o/R1/R2} $i ${i/R1/R2} > ${o/R1.fastq.gz/log} &
done
wait

nohup fastqc -o ${wdir/0.links/4.qc_cutadapt} ${wdir/0.links/3.cutadapt}/*.gz -t 16 &

######################################################################
#2018年7月23日
ddir=/home/share/DATA/clean_reads/GaoSR_mouse_nuclear_transfer_ChIP_20180301/20180719/1.rawdata
wdir=~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links/f.20180723
mkdir -p $wdir ${wdir/0.links/1.qc} ${wdir/0.links/3.cutadapt} ${wdir/0.links/4.qc_cutadapt}

ln -s $ddir/s0629-xrm-1_TKD180602711_1.fq.gz $wdir/6h_input.R1.fastq.gz
ln -s $ddir/s0629-xrm-2_TKD180602712_1.fq.gz $wdir/6h_K4me3.R1.fastq.gz
ln -s $ddir/s0629-xrm-3_TKD180602713_1.fq.gz $wdir/morula_input_rep1.R1.fastq.gz
ln -s $ddir/s0629-xrm-4_TKD180602714_1.fq.gz $wdir/morula_K4me3_rep1.R1.fastq.gz
ln -s $ddir/s0629-xrm-5_TKD180602715_1.fq.gz $wdir/morula_input_rep2.R1.fastq.gz
ln -s $ddir/s0629-xrm-6_TKD180602716_1.fq.gz $wdir/morula_K4me3_rep2.R1.fastq.gz
ln -s $ddir/s0629-xrm-1_TKD180602711_2.fq.gz $wdir/6h_input.R2.fastq.gz
ln -s $ddir/s0629-xrm-2_TKD180602712_2.fq.gz $wdir/6h_K4me3.R2.fastq.gz
ln -s $ddir/s0629-xrm-3_TKD180602713_2.fq.gz $wdir/morula_input_rep1.R2.fastq.gz
ln -s $ddir/s0629-xrm-4_TKD180602714_2.fq.gz $wdir/morula_K4me3_rep1.R2.fastq.gz
ln -s $ddir/s0629-xrm-5_TKD180602715_2.fq.gz $wdir/morula_input_rep2.R2.fastq.gz
ln -s $ddir/s0629-xrm-6_TKD180602716_2.fq.gz $wdir/morula_K4me3_rep2.R2.fastq.gz

nohup fastqc -o ${wdir/0.links/1.qc} $wdir/*gz -t 16 &
py35
for i in $wdir/*R1.fastq.gz; do o=${i/0.links/3.cutadapt}
nohup cutadapt -j 4 --trim-n -q 25,25 -m 75 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o $o -p ${o/R1/R2} $i ${i/R1/R2} > ${o/R1.fastq.gz/log} &
done
wait

nohup fastqc -o ${wdir/0.links/4.qc_cutadapt} ${wdir/0.links/3.cutadapt}/*.gz -t 16 &


######################################################################
#2018年8月04日
ddir1=/home/share/DATA/clean_reads/GaoSR_mouse_nuclear_transfer_ChIP_20180301/20180730/rawdata
ddir2=/home/share/DATA/clean_reads/GaoSR_mouse_nuclear_transfer_ChIP_20180301/20180801/rawdata
wdir=~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links/g.20180804
mkdir -p $wdir ${wdir/0.links/1.qc} ${wdir/0.links/3.cutadapt} ${wdir/0.links/4.qc_cutadapt}

ln -s $ddir1/s0713-XRM-1_TKD180701428_1.fq.gz $wdir/6h_input.R1.fastq.gz
ln -s $ddir1/s0713-XRM-2_TKD180701429_1.fq.gz $wdir/6h_K9me3.R1.fastq.gz
ln -s $ddir1/s0713-XRM-3_TKD180701430_1.fq.gz $wdir/morula_input.R1.fastq.gz
ln -s $ddir1/s0713-XRM-4_TKD180701431_1.fq.gz $wdir/morula_K9me3_rep1.R1.fastq.gz
ln -s $ddir1/s0713-XRM-5_TKD180701432_1.fq.gz $wdir/morula_K9me3_rep2.R1.fastq.gz
ln -s $ddir1/s0713-XRM-6_TKD180701433_1.fq.gz $wdir/4cell_input.R1.fastq.gz
ln -s $ddir1/s0713-XRM-7_TKD180701434_1.fq.gz $wdir/4cell_K9me3.R1.fastq.gz
ln -s $ddir2/0717-XRM-1_TKD180701692_1.fq.gz $wdir/2cell_input.R1.fastq.gz
ln -s $ddir2/0717-XRM-2_TKD180701693_1.fq.gz $wdir/2cell_K9me3.R1.fastq.gz
ln -s $ddir1/s0713-XRM-1_TKD180701428_2.fq.gz $wdir/6h_input.R2.fastq.gz
ln -s $ddir1/s0713-XRM-2_TKD180701429_2.fq.gz $wdir/6h_K9me3.R2.fastq.gz
ln -s $ddir1/s0713-XRM-3_TKD180701430_2.fq.gz $wdir/morula_input.R2.fastq.gz
ln -s $ddir1/s0713-XRM-4_TKD180701431_2.fq.gz $wdir/morula_K9me3_rep1.R2.fastq.gz
ln -s $ddir1/s0713-XRM-5_TKD180701432_2.fq.gz $wdir/morula_K9me3_rep2.R2.fastq.gz
ln -s $ddir1/s0713-XRM-6_TKD180701433_2.fq.gz $wdir/4cell_input.R2.fastq.gz
ln -s $ddir1/s0713-XRM-7_TKD180701434_2.fq.gz $wdir/4cell_K9me3.R2.fastq.gz
ln -s $ddir2/0717-XRM-1_TKD180701692_2.fq.gz $wdir/2cell_input.R2.fastq.gz
ln -s $ddir2/0717-XRM-2_TKD180701693_2.fq.gz $wdir/2cell_K9me3.R2.fastq.gz



nohup fastqc -o ${wdir/0.links/1.qc} $wdir/*gz -t 16 &
py35
for i in $wdir/*R1.fastq.gz; do o=${i/0.links/3.cutadapt}
nohup cutadapt -j 4 --trim-n -q 25,25 -m 75 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o $o -p ${o/R1/R2} $i ${i/R1/R2} > ${o/R1.fastq.gz/log} &
done
wait

nohup fastqc -o ${wdir/0.links/4.qc_cutadapt} ${wdir/0.links/3.cutadapt}/*.gz -t 16 &

http://10.11.41.110/new_home/qszhu/9.NT-ChIP/0.cutadapt-qc/1.qc/g.20180804/multiqc_report.html
http://10.11.41.110/new_home/qszhu/9.NT-ChIP/0.cutadapt-qc/4.qc_cutadapt/g.20180804/multiqc_report.html


################################################################################
##2018年9月3日到
cat *txt > temp.md5
md5sum -c temp.md5
for i in *1.fq.gz;do echo "ln -s \$ddir2/$i \$wdir/.R1.fastq.gz";done

ddir1=/home/share/DATA/clean_reads/GaoSR_mouse_nuclear_transfer_ChIP_20180301/20180809/data/rawdata
ddir2=/home/share/DATA/clean_reads/GaoSR_mouse_nuclear_transfer_ChIP_20180301/20180813/data/rawdata
wdir=~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links/h.20180909
mkdir -p $wdir ${wdir/0.links/1.qc} ${wdir/0.links/3.cutadapt} ${wdir/0.links/4.qc_cutadapt}

ln -s $ddir1/s0809-xrm-1_TKD180801015_1.fq.gz $wdir/6h_input_rep1.R1.fastq.gz
ln -s $ddir1/s0809-xrm-2_TKD180801016_1.fq.gz $wdir/14h_K4me3_rep1.R1.fastq.gz
ln -s $ddir1/s0809-xrm-3_TKD180801017_1.fq.gz $wdir/cc_input_rep1.R1.fastq.gz
ln -s $ddir1/s0809-xrm-4_TKD180801018_1.fq.gz $wdir/cc_input_rep2.R1.fastq.gz
ln -s $ddir2/s0813-xrm-1_TKD180801312_1.fq.gz $wdir/6h_K4me3_rep1.R1.fastq.gz
ln -s $ddir2/s0813-xrm-2_TKD180801313_1.fq.gz $wdir/14h_input_rep1.R1.fastq.gz
ln -s $ddir2/s0813-xrm-3_TKD180801314_1.fq.gz $wdir/cc_K27me3_rep1.R1.fastq.gz
ln -s $ddir2/s0813-xrm-4_TKD180801315_1.fq.gz $wdir/cc_K27me3_rep2.R1.fastq.gz
ln -s $ddir2/s0813-xrm-5_TKD180801316_1.fq.gz $wdir/cc_input_rep3.R1.fastq.gz
ln -s $ddir2/s0813-xrm-6_TKD180801317_1.fq.gz $wdir/cc_K9me3_rep1.R1.fastq.gz
ln -s $ddir2/s0813-xrm-7_TKD180801318_1.fq.gz $wdir/cc_input_rep4.R1.fastq.gz
ln -s $ddir2/s0813-xrm-8_TKD180801319_1.fq.gz $wdir/cc_K9me3_rep2.R1.fastq.gz
ln -s $ddir2/s0813-xrm-9_TKD180801320_1.fq.gz $wdir/cc_input_rep5.R1.fastq.gz
ln -s $ddir2/s0813-xrm-10_TKD180801321_1.fq.gz $wdir/cc_K4me3_rep1.R1.fastq.gz
ln -s $ddir2/s0813-xrm-11_TKD180801322_1.fq.gz $wdir/cc_input_rep6.R1.fastq.gz
ln -s $ddir2/s0813-xrm-12_TKD180801323_1.fq.gz $wdir/cc_K4me3_rep2.R1.fastq.gz
ln -s $ddir2/s0813-xrm-13_TKD180801301_1.fq.gz $wdir/te_input_rep1.R1.fastq.gz
ln -s $ddir2/s0813-xrm-14_TKD180801302_1.fq.gz $wdir/te_K9me3_rep1.R1.fastq.gz
ln -s $ddir2/s0813-xrm-15_TKD180801303_1.fq.gz $wdir/icm_input_rep1.R1.fastq.gz
ln -s $ddir2/s0813-xrm-16_TKD180801304_1.fq.gz $wdir/icm_K9me3_rep1.R1.fastq.gz
ln -s $ddir2/s0813-xrm-17_TKD180801305_1.fq.gz $wdir/blast_input_rep1.R1.fastq.gz
ln -s $ddir2/s0813-xrm-18_TKD180801306_1.fq.gz $wdir/blast_K9me3_rep1.R1.fastq.gz
ln -s $ddir2/s0813-xrm-19_TKD180801307_1.fq.gz $wdir/morula_input_rep1.R1.fastq.gz
ln -s $ddir2/s0813-xrm-20_TKD180801308_1.fq.gz $wdir/morula_K27me3_rep1.R1.fastq.gz
ln -s $ddir2/s0813-xrm-21_TKD180801309_1.fq.gz $wdir/morula_K27me3_rep2.R1.fastq.gz

ln -s $ddir1/s0809-xrm-1_TKD180801015_2.fq.gz $wdir/6h_input_rep1.R2.fastq.gz
ln -s $ddir1/s0809-xrm-2_TKD180801016_2.fq.gz $wdir/14h_K4me3_rep1.R2.fastq.gz
ln -s $ddir1/s0809-xrm-3_TKD180801017_2.fq.gz $wdir/cc_input_rep1.R2.fastq.gz
ln -s $ddir1/s0809-xrm-4_TKD180801018_2.fq.gz $wdir/cc_input_rep2.R2.fastq.gz
ln -s $ddir2/s0813-xrm-1_TKD180801312_2.fq.gz $wdir/6h_K4me3_rep1.R2.fastq.gz
ln -s $ddir2/s0813-xrm-2_TKD180801313_2.fq.gz $wdir/14h_input_rep1.R2.fastq.gz
ln -s $ddir2/s0813-xrm-3_TKD180801314_2.fq.gz $wdir/cc_K27me3_rep1.R2.fastq.gz
ln -s $ddir2/s0813-xrm-4_TKD180801315_2.fq.gz $wdir/cc_K27me3_rep2.R2.fastq.gz
ln -s $ddir2/s0813-xrm-5_TKD180801316_2.fq.gz $wdir/cc_input_rep3.R2.fastq.gz
ln -s $ddir2/s0813-xrm-6_TKD180801317_2.fq.gz $wdir/cc_K9me3_rep1.R2.fastq.gz
ln -s $ddir2/s0813-xrm-7_TKD180801318_2.fq.gz $wdir/cc_input_rep4.R2.fastq.gz
ln -s $ddir2/s0813-xrm-8_TKD180801319_2.fq.gz $wdir/cc_K9me3_rep2.R2.fastq.gz
ln -s $ddir2/s0813-xrm-9_TKD180801320_2.fq.gz $wdir/cc_input_rep5.R2.fastq.gz
ln -s $ddir2/s0813-xrm-10_TKD180801321_2.fq.gz $wdir/cc_K4me3_rep1.R2.fastq.gz
ln -s $ddir2/s0813-xrm-11_TKD180801322_2.fq.gz $wdir/cc_input_rep6.R2.fastq.gz
ln -s $ddir2/s0813-xrm-12_TKD180801323_2.fq.gz $wdir/cc_K4me3_rep2.R2.fastq.gz
ln -s $ddir2/s0813-xrm-13_TKD180801301_2.fq.gz $wdir/te_input_rep1.R2.fastq.gz
ln -s $ddir2/s0813-xrm-14_TKD180801302_2.fq.gz $wdir/te_K9me3_rep1.R2.fastq.gz
ln -s $ddir2/s0813-xrm-15_TKD180801303_2.fq.gz $wdir/icm_input_rep1.R2.fastq.gz
ln -s $ddir2/s0813-xrm-16_TKD180801304_2.fq.gz $wdir/icm_K9me3_rep1.R2.fastq.gz
ln -s $ddir2/s0813-xrm-17_TKD180801305_2.fq.gz $wdir/blast_input_rep1.R2.fastq.gz
ln -s $ddir2/s0813-xrm-18_TKD180801306_2.fq.gz $wdir/blast_K9me3_rep1.R2.fastq.gz
ln -s $ddir2/s0813-xrm-19_TKD180801307_2.fq.gz $wdir/morula_input_rep1.R2.fastq.gz
ln -s $ddir2/s0813-xrm-20_TKD180801308_2.fq.gz $wdir/morula_K27me3_rep1.R2.fastq.gz
ln -s $ddir2/s0813-xrm-21_TKD180801309_2.fq.gz $wdir/morula_K27me3_rep2.R2.fastq.gz

# ln -s $ddir2/s0813-yly-1_TKD180801310_1.fq.gz $wdir/.R1.fastq.gz
# ln -s $ddir2/s0813-yly-2_TKD180801311_1.fq.gz $wdir/.R1.fastq.gz

nohup fastqc -o ${wdir/0.links/1.qc} $wdir/*gz -t 16 &
py35
for i in $wdir/*R1.fastq.gz; do o=${i/0.links/3.cutadapt}
nohup cutadapt -j 4 --trim-n -q 25,25 -m 75 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o $o -p ${o/R1/R2} $i ${i/R1/R2} > ${o/R1.fastq.gz/log} &
done
wait

nohup fastqc -o ${wdir/0.links/4.qc_cutadapt} ${wdir/0.links/3.cutadapt}/*.gz -t 32 &

http://10.11.41.110/new_home/qszhu/9.NT-ChIP/0.cutadapt-qc/1.qc/h.20180909/multiqc_report.html
http://10.11.41.110/new_home/qszhu/9.NT-ChIP/0.cutadapt-qc/4.qc_cutadapt/h.20180909/multiqc_report.html



###########################
#名称检验
ll a*/*R1*gz | awk 'BEGIN{FS="[ /.]"} {$1=$1;print $8,$9,$10,$23}' > links.tab
ll b*/*R1*gz | awk 'BEGIN{FS="[ /.]"} {$1=$1;print $9,$10,$11,$23}' >> links.tab
ll c*/*R1*gz | awk 'BEGIN{FS="[ /.]"} {$1=$1;print $9,$10,$11,$24}' >> links.tab
ll d*/*R1*gz | awk 'BEGIN{FS="[ /.]"} {$1=$1;print $9,$10,$11,$23}' >> links.tab
ll e*/*R1*gz | awk 'BEGIN{FS="[ /.]"} {$1=$1;print $10,$11,$12,$26}' >> links.tab
ll f*/*R1*gz | awk 'BEGIN{FS="[ /.]"} {$1=$1;print  $9,$10,$11,$25}' >> links.tab
ll g*/*R1*gz | awk 'BEGIN{FS="[ /.]"} {$1=$1;print $10,$11,$12,$25}' >> links.tab
ll h*/*R1*gz | awk 'BEGIN{FS="[ /.]"} {$1=$1;print $9,$10,$11,$25}' >> links.tab
ll i*/*R1*gz | awk 'BEGIN{FS="[ /.]"} {$1=$1;print $9,$10,$11,$23}' >> links.tab
ll j*/*R1*gz | awk 'BEGIN{FS="[ /.]"} {$1=$1;print $8,$9,$10,$22}' >> links.tab
ll k*/*R1*gz | awk 'BEGIN{FS="[ /.]"} {$1=$1;print $8,$9,$10,$22}' >> links.tab
##################################################

################################################################################
##2018年11月09日到   1005 1006 1008 1017测序
cat *txt > temp.md5
md5sum -c temp.md5
for i in *1.fq.gz;do echo "ln -s \$ddir2/$i \$wdir/.R1.fastq.gz";done

ddir=/home/share/DATA/clean_reads/GaoSR_mouse_nuclear_transfer_ChIP_20180301/20181109get
wdir=~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links/i.20181109
mkdir -p $wdir ${wdir/0.links/1.qc} ${wdir/0.links/3.cutadapt} ${wdir/0.links/4.qc_cutadapt}

ln -s $ddir/s1005-xrm-1_TKD181000143_HLKHTCCXY_L8_1.fq.gz $wdir/8cell_K27me3.R1.fastq.gz
ln -s $ddir/s1006-xrm-1_TKD181000156_1.fq.gz $wdir/2cell_K27me3.R1.fastq.gz
ln -s $ddir/s1006-xrm-2_TKD181000157_1.fq.gz $wdir/4cell_K27me3.R1.fastq.gz
ln -s $ddir/s1008-xrm-1_TKD181000650_HVLGMCCXY_L4_1.fq.gz $wdir/6h_input.R1.fastq.gz
ln -s $ddir/s1008-xrm-2_TKD181000651_HVLGMCCXY_L4_1.fq.gz $wdir/6h_K9me3.R1.fastq.gz
ln -s $ddir/s1017-xrm-1_TKD181001497_HT3MYCCXY_L4_1.fq.gz $wdir/8cell_input.R1.fastq.gz
ln -s $ddir/s1017-xrm-2_TKD181001498_HT3MYCCXY_L4_1.fq.gz $wdir/8cell_K9me3.R1.fastq.gz
ln -s $ddir/s1017-xrm-4_TKD181001499_HT3MYCCXY_L4_1.fq.gz $wdir/4cell_K4me3.R1.fastq.gz

ln -s $ddir/s1005-xrm-1_TKD181000143_HLKHTCCXY_L8_2.fq.gz $wdir/8cell_K27me3.R2.fastq.gz
ln -s $ddir/s1006-xrm-1_TKD181000156_2.fq.gz $wdir/2cell_K27me3.R2.fastq.gz
ln -s $ddir/s1006-xrm-2_TKD181000157_2.fq.gz $wdir/4cell_K27me3.R2.fastq.gz
ln -s $ddir/s1008-xrm-1_TKD181000650_HVLGMCCXY_L4_2.fq.gz $wdir/6h_input.R2.fastq.gz
ln -s $ddir/s1008-xrm-2_TKD181000651_HVLGMCCXY_L4_2.fq.gz $wdir/6h_K9me3.R2.fastq.gz
ln -s $ddir/s1017-xrm-1_TKD181001497_HT3MYCCXY_L4_2.fq.gz $wdir/8cell_input.R2.fastq.gz
ln -s $ddir/s1017-xrm-2_TKD181001498_HT3MYCCXY_L4_2.fq.gz $wdir/8cell_K9me3.R2.fastq.gz
ln -s $ddir/s1017-xrm-4_TKD181001499_HT3MYCCXY_L4_2.fq.gz $wdir/4cell_K4me3.R2.fastq.gz
# ln -s $ddir2/s0813-yly-1_TKD180801310_1.fq.gz $wdir/.R1.fastq.gz
# ln -s $ddir2/s0813-yly-2_TKD180801311_1.fq.gz $wdir/.R1.fastq.gz

nohup fastqc -o ${wdir/0.links/1.qc} $wdir/*gz -t 16 &
py35
for i in $wdir/*R1.fastq.gz; do o=${i/0.links/3.cutadapt}
nohup cutadapt -j 4 --trim-n -q 25,25 -m 75 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o $o -p ${o/R1/R2} $i ${i/R1/R2} > ${o/R1.fastq.gz/log} &
done
wait

nohup fastqc -o ${wdir/0.links/4.qc_cutadapt} ${wdir/0.links/3.cutadapt}/*.gz -t 32 &

http://10.11.41.110/new_home/qszhu/9.NT-ChIP/0.cutadapt-qc/1.qc/i.20181109/multiqc_report.html
http://10.11.41.110/new_home/qszhu/9.NT-ChIP/0.cutadapt-qc/4.qc_cutadapt/i.20181109/multiqc_report.html


#################################################################
##2018年12月14日   1126测序
cat *txt > temp.md5
md5sum -c temp.md5
for i in *1.fq.gz;do echo "ln -s \$ddir2/$i \$wdir/.R1.fastq.gz";done

ddir=/home/share/DATA/clean_reads/GaoSR_mouse_nuclear_transfer_ChIP_20180301/20181214
wdir=~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links/j.20181214
mkdir -p $wdir ${wdir/0.links/1.qc} ${wdir/0.links/3.cutadapt} ${wdir/0.links/4.qc_cutadapt}
ln -s $ddir/s1028-XRM-2_TKD181102559_htgwgccxy_L3_1.fq.gz $wdir/4cell_K4me3_rep3.R1.fastq.gz
ln -s $ddir/s1028-XRM-4_TKD181102560_htgwgccxy_L3_1.fq.gz $wdir/8cell_K9me3_rep3.R1.fastq.gz
ln -s $ddir/s1028-XRM-5_TKD181102588_htgwgccxy_L2_1.fq.gz $wdir/8cell_K4me3_rep2.R1.fastq.gz
ln -s $ddir/s1122-XRM-10_TKD181102589_htgwgccxy_L2_1.fq.gz $wdir/icm_K4me3_rep1.R1.fastq.gz
ln -s $ddir/s1122-XRM-12_TKD181102584_htgwgccxy_L2_1.fq.gz $wdir/14h_K9me3_rep2.R1.fastq.gz
ln -s $ddir/s1122-XRM-14_TKD181102585_htgwgccxy_L2_1.fq.gz $wdir/8cell_K9me3_rep4.R1.fastq.gz
ln -s $ddir/s1122-XRM-2_TKD181102586_htgwgccxy_L2_1.fq.gz $wdir/2cell_K4me3_rep2.R1.fastq.gz
ln -s $ddir/s1122-XRM-4_TKD181102587_htgwgccxy_L2_1.fq.gz $wdir/4cell_K4me3_rep4.R1.fastq.gz
ln -s $ddir/s1122-XRM-6_TKD181102590_htgwgccxy_L2_1.fq.gz $wdir/blast_K4me3_rep1.R1.fastq.gz
ln -s $ddir/s1122-XRM-8_TKD181102561_htgwgccxy_L3_1.fq.gz $wdir/te_K4me3_rep1.R1.fastq.gz
#ln -s $ddir/s1123-XRM-4_TKD181102562_htgwgccxy_L3_1.fq.gz $wdir/6h_K4me3_rep3.R1.fastq.gz #这个命名错误


ln -s $ddir/s1028-XRM-2_TKD181102559_htgwgccxy_L3_2.fq.gz $wdir/4cell_K4me3_rep3.R2.fastq.gz
ln -s $ddir/s1028-XRM-4_TKD181102560_htgwgccxy_L3_2.fq.gz $wdir/8cell_K9me3_rep3.R2.fastq.gz
ln -s $ddir/s1028-XRM-5_TKD181102588_htgwgccxy_L2_2.fq.gz $wdir/8cell_K4me3_rep2.R2.fastq.gz
ln -s $ddir/s1122-XRM-10_TKD181102589_htgwgccxy_L2_2.fq.gz $wdir/icm_K4me3_rep1.R2.fastq.gz
ln -s $ddir/s1122-XRM-12_TKD181102584_htgwgccxy_L2_2.fq.gz $wdir/14h_K9me3_rep2.R2.fastq.gz
ln -s $ddir/s1122-XRM-14_TKD181102585_htgwgccxy_L2_2.fq.gz $wdir/8cell_K9me3_rep4.R2.fastq.gz
ln -s $ddir/s1122-XRM-2_TKD181102586_htgwgccxy_L2_2.fq.gz $wdir/2cell_K4me3_rep2.R2.fastq.gz
ln -s $ddir/s1122-XRM-4_TKD181102587_htgwgccxy_L2_2.fq.gz $wdir/4cell_K4me3_rep4.R2.fastq.gz
ln -s $ddir/s1122-XRM-6_TKD181102590_htgwgccxy_L2_2.fq.gz $wdir/blast_K4me3_rep1.R2.fastq.gz
ln -s $ddir/s1122-XRM-8_TKD181102561_htgwgccxy_L3_2.fq.gz $wdir/te_K4me3_rep1.R2.fastq.gz
ln -s $ddir/s1123-XRM-4_TKD181102562_htgwgccxy_L3_2.fq.gz $wdir/6h_K4me3_rep3.R2.fastq.gz

nohup fastqc -o ${wdir/0.links/1.qc} $wdir/*gz -t 16 &
py35
for i in $wdir/*R1.fastq.gz; do o=${i/0.links/3.cutadapt}
nohup cutadapt -j 4 --trim-n -q 25,25 -m 75 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o $o -p ${o/R1/R2} $i ${i/R1/R2} > ${o/R1.fastq.gz/log} &
done
wait
nohup fastqc -o ${wdir/0.links/4.qc_cutadapt} ${wdir/0.links/3.cutadapt}/*.gz -t 32 &

#################################################################
##2019年02月26日   20181229 20190104 20190115 20190119 20190121
cat *txt | md5sum -c -
for i in *1.fq.gz;do i=${i/_1.fq.gz};echo "ln -s \$ddir/${i}_\${i}.fq.gz \$wdir/.R\${i}.fastq.gz";done

ddir=/home/share/DATA/clean_reads/GaoSR_mouse_nuclear_transfer_ChIP_20180301/20190226get
wdir=~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links/k.20190228
mkdir -p $wdir ${wdir/0.links/1.qc} ${wdir/0.links/3.cutadapt} ${wdir/0.links/4.qc_cutadapt}
for i in 1 2; do 
ln -s $ddir/0117-xrm-14_FKDL190724273-1a_HY2FWCCXY_L3_${i}.fq.gz $wdir/blast_K4me3_rep2.R${i}.fastq.gz
ln -s $ddir/s0115-xrm-6_FKDL190664690-1a_HY2FWCCXY_L1_${i}.fq.gz $wdir/14h_input_rep2.R${i}.fastq.gz
ln -s $ddir/s0117-xrm-1_FKDL190665516-1a_HYF72CCXY_L6_${i}.fq.gz $wdir/14h_K9me3_rep3.R${i}.fastq.gz
ln -s $ddir/s0117-xrm-2_FKDL190665520-1a_HYF72CCXY_L6_${i}.fq.gz $wdir/6h_input_rep2.R${i}.fastq.gz
ln -s $ddir/s0117-xrm-3_FKDL190665521-1a_HYF72CCXY_L6_${i}.fq.gz $wdir/6h_K27me3_rep1.R${i}.fastq.gz
ln -s $ddir/s0117-xrm-4_FKDL190665522-1a_HYF72CCXY_L6_${i}.fq.gz $wdir/te_K4me3_rep2.R${i}.fastq.gz
ln -s $ddir/s0117-xrm-5_FKDL190665523-1a_HYF72CCXY_L6_${i}.fq.gz $wdir/icm_K4me3_rep2.R${i}.fastq.gz
ln -s $ddir/s1129-XRM-1_TKD181203597_HWH33CCXY_L4_${i}.fq.gz $wdir/6h_K4me3_rep4.R${i}.fastq.gz
ln -s $ddir/s1129-XRM-2_TKD181203598_HWH33CCXY_L4_${i}.fq.gz $wdir/2cell_K27me3_rep2.R${i}.fastq.gz
ln -s $ddir/s1129-XRM-3_TKD181203599_HWH33CCXY_L4_${i}.fq.gz $wdir/8cell_K27me3_rep2.R${i}.fastq.gz
ln -s $ddir/s1230-xrm-1_TKD190100358_HWGGKCCXY_L8_${i}.fq.gz $wdir/14h_input_rep3.R${i}.fastq.gz
ln -s $ddir/s1230-xrm-2_TKD190100359_HWGGKCCXY_L8_${i}.fq.gz $wdir/14h_K4me3_rep3.R${i}.fastq.gz
done

nohup fastqc -o ${wdir/0.links/1.qc} $wdir/*gz -t 16 &
py35
for i in $wdir/*R1.fastq.gz; do o=${i/0.links/3.cutadapt}
nohup cutadapt -j 4 --trim-n -q 25,25 -m 75 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o $o -p ${o/R1/R2} $i ${i/R1/R2} > ${o/R1.fastq.gz/log} &
done
wait
nohup fastqc -o ${wdir/0.links/4.qc_cutadapt} ${wdir/0.links/3.cutadapt}/*.gz -t 32 &




###################################
#2019年3月11日 K9me3数据完整，后面重新比对使用此命名
wdir=~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links
keys=(cc 6h 14h 2cell 4cell 8cell morula icm te blast)
for k in ${keys[@]}; do let n=1
for i in $wdir/[a-k]*/${k}_K9me3*.R1.fastq.gz; do 
cp -r $i $wdir/1.K9me3/${k}_K9me3_rep${n}.R1.fastq.gz
cp -r ${i/R1/R2} $wdir/1.K9me3/${k}_K9me3_rep${n}.R2.fastq.gz
let n++
done
done


wdir=~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links
keys=(cc 6h 14h 2cell 4cell 8cell morula icm te blast)
for k in ${keys[@]}; do let n=1
for i in $wdir/[a-k]*/${k}_input*.R1.fastq.gz; do 
cp -r $i $wdir/0.input/${k}_input_rep${n}.R1.fastq.gz
cp -r ${i/R1/R2} $wdir/0.input/${k}_input_rep${n}.R2.fastq.gz
let n++
done
done



###################2019年3月29日 K27 K4 
cat *txt | md5sum -c -
for i in *1.fq.gz;do i=${i/_R1.fq.gz};echo "ln -s \$ddir/${i}_\${i}.fq.gz \$wdir/.R\${i}.fastq.gz";done

ddir=/home/share/DATA/clean_reads/GaoSR_mouse_nuclear_transfer_ChIP_20180301/20190325
wdir=~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links/l.20190325
mkdir -p $wdir ${wdir/0.links/1.qc} ${wdir/0.links/3.cutadapt} ${wdir/0.links/4.qc_cutadapt}
for i in 1 2; do 
ln -s $ddir/0126-xrm-4_R${i}.fq.gz $wdir/6h_K27me3_rep2.R${i}.fastq.gz
ln -s $ddir/0126-xrm-5_R${i}.fq.gz $wdir/te_K27me3_rep1.R${i}.fastq.gz
ln -s $ddir/0126-xrm-6_R${i}.fq.gz $wdir/icm_K27me3_rep1.R${i}.fastq.gz
ln -s $ddir/0126-xrm-7_R${i}.fq.gz $wdir/blast_K27me3_rep1.R${i}.fastq.gz
ln -s $ddir/0222-xrm-1_R${i}.fq.gz $wdir/te_K27me3_rep2.R${i}.fastq.gz
ln -s $ddir/0222-xrm-2_R${i}.fq.gz $wdir/icm_K27me3_rep2.R${i}.fastq.gz
ln -s $ddir/0222-xrm-3_R${i}.fq.gz $wdir/blast_K27me3_rep2.R${i}.fastq.gz
ln -s $ddir/0222-xrm-4_R${i}.fq.gz $wdir/14h_K27me3_rep1.R${i}.fastq.gz
ln -s $ddir/0303-xrm-1_R${i}.fq.gz $wdir/14h_K27me3_rep2.R${i}.fastq.gz
ln -s $ddir/0303-xrm-2_R${i}.fq.gz $wdir/6h_K27me3_rep3.R${i}.fastq.gz
ln -s $ddir/0312-xrm-1_R${i}.fq.gz $wdir/14h_K4me3_rep4.R${i}.fastq.gz
ln -s $ddir/0312-xrm-2_R${i}.fq.gz $wdir/blast_K4me3_rep3.R${i}.fastq.gz
done

nohup fastqc -o ${wdir/0.links/1.qc} $wdir/*gz -t 64 &
py35
#不用去adapt，跑BWA
######################################重新命名所有K27me3
wdir=~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links
keys=(cc 6h 14h 2cell 4cell 8cell morula icm te blast)
for k in ${keys[@]}; do let n=1
for i in $wdir/[a-l]*/${k}_K27me3*.R1.fastq.gz; do 
cp -r $i $wdir/3.K27me3/${k}_K27me3_rep${n}.R1.fastq.gz
cp -r ${i/R1/R2} $wdir/3.K27me3/${k}_K27me3_rep${n}.R2.fastq.gz
let n++
done
done
#更正一个命名错误的数据
cp -r $wdir/j.20181214/6h_K4me3_rep3.R1.fastq.gz $wdir/3.K27me3/4cell_K27me3_rep3.R1.fastq.gz
cp -r $wdir/j.20181214/6h_K4me3_rep3.R2.fastq.gz $wdir/3.K27me3/4cell_K27me3_rep3.R2.fastq.gz
#统计数据名称关系
ll *R1*gz | awk 'BEGIN{FS="[ /.]"} {$1=$1;print $8,$9,$10,$22}' >> links.tab
ll l*/*R1*gz | awk 'BEGIN{FS="[ /.]"} {$1=$1;print $8,$9,$10,$22}' >> links.tab

############################################重新命名所有K4me3
wdir=~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links
keys=(cc 6h 14h 2cell 4cell 8cell morula icm te blast)
for k in ${keys[@]}; do let n=1
for i in $wdir/[a-l]*/${k}_K4me3*.R1.fastq.gz; do 
cp -r $i $wdir/2.K4me3/${k}_K4me3_rep${n}.R1.fastq.gz
cp -r ${i/R1/R2} $wdir/2.K4me3/${k}_K4me3_rep${n}.R2.fastq.gz
let n++
done
done
#更正一个命名错误的数据
rm $wdir/2.K4me3/6h_K4me3_rep4*
mv $wdir/2.K4me3/6h_K4me3_rep5.R1.fastq.gz $wdir/2.K4me3/6h_K4me3_rep4.R1.fastq.gz
mv $wdir/2.K4me3/6h_K4me3_rep5.R2.fastq.gz $wdir/2.K4me3/6h_K4me3_rep4.R2.fastq.gz

ll *R1*gz | awk 'BEGIN{FS="[ /._]"} {$1=$1;print $8,$9,$10,$28,$29,$30,$31}'



###################################20190530 H3.3 HA & cc-K9me3
cat *txt | md5sum -c -
for i in *1.fq.gz;do i=${i/_R1.fq.gz};echo "ln -s \$ddir/${i}_\${i}.fq.gz \$wdir/.R\${i}.fastq.gz";done

ddir=/home/share/DATA/clean_reads/GaoSR_mouse_nuclear_transfer_ChIP_20180301/20190530
wdir=~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links/m.20190530
mkdir -p $wdir ${wdir/0.links/1.qc} ${wdir/0.links/3.cutadapt} ${wdir/0.links/4.qc_cutadapt}
for i in 1 2; do 
ln -s $ddir/0427-xrm-1_*${i}.fq.gz $wdir/6h_input_rep9.R${i}.fastq.gz  #HA
ln -s $ddir/0427-xrm-2_*${i}.fq.gz $wdir/6h_HA_rep1.R${i}.fastq.gz
ln -s $ddir/0427-xrm-3_*${i}.fq.gz $wdir/cc_input_rep8.R${i}.fastq.gz  #HA
ln -s $ddir/0427-xrm-4_*${i}.fq.gz $wdir/cc_input_rep9.R${i}.fastq.gz  #HA
ln -s $ddir/0427-xrm-5_*${i}.fq.gz $wdir/cc_HA_rep1.R${i}.fastq.gz
ln -s $ddir/0427-xrm-6_*${i}.fq.gz $wdir/cc_HA_rep2.R${i}.fastq.gz
ln -s $ddir/0427-xrm-7_*${i}.fq.gz $wdir/cc_input_rep7.R${i}.fastq.gz
ln -s $ddir/0427-xrm-8_*${i}.fq.gz $wdir/cc_K9me3_rep3.R${i}.fastq.gz
ln -s $ddir/0427-xrm-9_*${i}.fq.gz $wdir/cc_K9me3_rep4.R${i}.fastq.gz
done
nohup fastqc -o ${wdir/0.links/1.qc} $wdir/*gz -t 64 &



###################################20190530 H3.3 HA & cc-K9me3
cat *txt | md5sum -c -
for i in *1.fq.gz;do i=${i/_1.fq.gz};echo "ln -s \$ddir/${i}_\${i}.fq.gz \$wdir/.R\${i}.fastq.gz";done

ddir=/home/share/DATA/clean_reads/GaoSR_mouse_nuclear_transfer_ChIP_20180301/20190706
wdir=~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links/n.20190709
mkdir -p $wdir ${wdir/0.links/1.qc} ${wdir/0.links/3.cutadapt} ${wdir/0.links/4.qc_cutadapt}
for i in 1 2; do 
ln -s $ddir/0519-xrm-1_FKDL190739577-1a_H2NTGCCX2_L4_${i}.fq.gz $wdir/6h_input_repa.R${i}.fastq.gz #HA
ln -s $ddir/0519-xrm-2_FKDL190739578-1a_H2NTGCCX2_L4_${i}.fq.gz $wdir/6h_HA_rep2.R${i}.fastq.gz
ln -s $ddir/0606-xrm-1_FKDL190742150-1a_H2WFLCCX2_L3_${i}.fq.gz $wdir/PN5_input_rep1.R${i}.fastq.gz
ln -s $ddir/0606-xrm-2_FKDL190742151-1a_H32NHCCX2_L5_${i}.fq.gz $wdir/PN5_K9me3_rep1.R${i}.fastq.gz
ln -s $ddir/0606-xrm-3_FKDL190742152-1a_H32NHCCX2_L5_${i}.fq.gz $wdir/PN5_input_rep2.R${i}.fastq.gz
ln -s $ddir/0606-xrm-4_FKDL190742153-1a_H32NHCCX2_L5_${i}.fq.gz $wdir/PN5_K9me3_rep2.R${i}.fastq.gz
ln -s $ddir/0624-xrm-1_FKDL190745450-1a_H32NHCCX2_L5_${i}.fq.gz $wdir/14h_input_rep7.R${i}.fastq.gz
ln -s $ddir/0624-xrm-2_FKDL190745451-1a_H32NHCCX2_L5_${i}.fq.gz $wdir/14h_HA_rep1.R${i}.fastq.gz
ln -s $ddir/0624-xrm-3_FKDL190745452-1a_H32NHCCX2_L5_${i}.fq.gz $wdir/14h_HA_rep2.R${i}.fastq.gz
done
nohup fastqc -o ${wdir/0.links/1.qc} $wdir/*gz -t 64 &

for i in input K9me3 HA; do 
cp -r ~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links/[mn]*/*$i*gz ~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links/*$i/
done

for i in input K9me3 HA; do 
mkd $wdir/$i
cp -r $wdir/*$i*gz $wdir/$i/
done

####bwa
wdir=~/workspace/9.NT-ChIP/1.align/7.HA_BWA
ddir=~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links/n.20190709/HA
bash ~/workspace/9.NT-ChIP/1.align/thread_bwa.sh -w $wdir -d $ddir -t 8 -p 4 &

wdir=~/workspace/9.NT-ChIP/1.align/1.K9me3_BWA
ddir=~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links/n.20190709/K9me3
bash ~/workspace/9.NT-ChIP/1.align/thread_bwa.sh -w $wdir -d $ddir -t 8 -p 4 & #
wdir=~/workspace/9.NT-ChIP/1.align/0.input_BWA
ddir=~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links/n.20190709/input
bash ~/workspace/9.NT-ChIP/1.align/thread_bwa.sh -w $wdir -d $ddir -t 8 -p 4 & #

#####################################################################################################
###2019年10月21日 RNA smart-seq
cat *txt | md5sum -c -
for i in *1.fq.gz;do i=${i/_1.fq.gz};echo "ln -s \$ddir/${i}_\${i}.fq.gz \$wdir/.R\${i}.fastq.gz";done

ddir=~/data/GaoSR_mouse_nuclear_transfer_ChIP_20180301/20191020-smart-seq
wdir=~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links/o.20191021
mkdir -p $wdir ${wdir/0.links/1.qc} ${wdir/0.links/3.cutadapt} ${wdir/0.links/4.qc_cutadapt}
for i in 1 2; do 
ln -s $ddir/0722-yly-12_FKDL190751173-1a_H32CMCCX2_L8_${i}.fq.gz $wdir/6h_RNA_rep2.R${i}.fastq.gz
ln -s $ddir/0722-yly-16_FKDL190751177-1a_H32CMCCX2_L8_${i}.fq.gz $wdir/6h_RNA_rep3.R${i}.fastq.gz
ln -s $ddir/0722-yly-19_FKDL190751180-1a_H32CMCCX2_L8_${i}.fq.gz $wdir/l2cell_RNA_rep3.R${i}.fastq.gz
ln -s $ddir/0722-yly-1_FKDL190751162-1a_H32CMCCX2_L8_${i}.fq.gz $wdir/6h_RNA_rep1.R${i}.fastq.gz
ln -s $ddir/0722-yly-4_FKDL190751165-1a_H32CMCCX2_L8_${i}.fq.gz $wdir/l2cell_RNA_rep2.R${i}.fastq.gz
ln -s $ddir/0722-yly-7_FKDL190751168-1a_H32CMCCX2_L8_${i}.fq.gz $wdir/l2cell_RNA_rep1.R${i}.fastq.gz
ln -s $ddir/0923-xrm-1_FKDL190765346-1a_${i}.fq.gz $wdir/14h_RNA_rep1.R${i}.fastq.gz
ln -s $ddir/0923-xrm-2_FKDL190765347-1a_${i}.fq.gz $wdir/14h_RNA_rep2.R${i}.fastq.gz
ln -s $ddir/0923-xrm-3_FKDL190765348-1a_${i}.fq.gz $wdir/14h_RNA_rep3.R${i}.fastq.gz
ln -s $ddir/0923-xrm-4_FKDL190765349-1a_${i}.fq.gz $wdir/e2cell_RNA_rep1.R${i}.fastq.gz
ln -s $ddir/0923-xrm-5_FKDL190765350-1a_${i}.fq.gz $wdir/e2cell_RNA_rep2.R${i}.fastq.gz
ln -s $ddir/0923-xrm-6_FKDL190765351-1a_${i}.fq.gz $wdir/e2cell_RNA_rep3.R${i}.fastq.gz
ln -s $ddir/0923-xrm-7_FKDL190765352-1a_${i}.fq.gz $wdir/4cell_RNA_rep1.R${i}.fastq.gz
ln -s $ddir/0923-xrm-8_FKDL190765353-1a_${i}.fq.gz $wdir/4cell_RNA_rep2.R${i}.fastq.gz
ln -s $ddir/0923-xrm-9_FKDL190765354-1a_${i}.fq.gz $wdir/4cell_RNA_rep3.R${i}.fastq.gz
ln -s $ddir/0923-xrm-10_FKDL190765355-1a_${i}.fq.gz $wdir/8cell_RNA_rep1.R${i}.fastq.gz
ln -s $ddir/0923-xrm-11_FKDL190765356-1a_${i}.fq.gz $wdir/8cell_RNA_rep2.R${i}.fastq.gz
ln -s $ddir/0923-xrm-12_FKDL190765357-1a_${i}.fq.gz $wdir/8cell_RNA_rep3.R${i}.fastq.gz
ln -s $ddir/0923-xrm-13_FKDL190765358-1a_${i}.fq.gz $wdir/morula_RNA_rep1.R${i}.fastq.gz
ln -s $ddir/0923-xrm-14_FKDL190765359-1a_${i}.fq.gz $wdir/morula_RNA_rep2.R${i}.fastq.gz
ln -s $ddir/0923-xrm-15_FKDL190765360-1a_${i}.fq.gz $wdir/morula_RNA_rep3.R${i}.fastq.gz
ln -s $ddir/0923-xrm-16_FKDL190765361-1a_${i}.fq.gz $wdir/icm_RNA_rep1.R${i}.fastq.gz
ln -s $ddir/0923-xrm-17_FKDL190765362-1a_${i}.fq.gz $wdir/icm_RNA_rep2.R${i}.fastq.gz
ln -s $ddir/0923-xrm-18_FKDL190765363-1a_${i}.fq.gz $wdir/icm_RNA_rep3.R${i}.fastq.gz
ln -s $ddir/0923-xrm-19_FKDL190765364-1a_${i}.fq.gz $wdir/te_RNA_rep1.R${i}.fastq.gz
ln -s $ddir/0923-xrm-20_FKDL190765365-1a_${i}.fq.gz $wdir/te_RNA_rep2.R${i}.fastq.gz
ln -s $ddir/0923-xrm-21_FKDL190765366-1a_${i}.fq.gz $wdir/te_RNA_rep3.R${i}.fastq.gz
done
nohup fastqc -o ${wdir/0.links/1.qc} $wdir/*gz -t 64 &

cp -r ~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links/o.20191021/*RNA*gz ~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links/5.RNA/

py35
for i in $wdir/*R1.fastq.gz; do o=${i/0.links/3.cutadapt}
nohup cutadapt -j 4 --trim-n -q 25,25 -m 30 -a ATCAACGCAGAGTAC -a AGATCGGAAGAGC -A ATCAACGCAGAGTAC -A AGATCGGAAGAGC -o $o -p ${o/R1/R2} $i ${i/R1/R2} > ${o/R1.fastq.gz/log} &
done
wait
nohup fastqc -o ${wdir/0.links/4.qc_cutadapt} ${wdir/0.links/3.cutadapt}/*.gz -t 32 &
#需要去掉smart kit adapt，不然比对率很低，但最终得到的数据是差不多的

http://guru/new_home/qszhu/9.NT-ChIP/0.cutadapt-qc/1.qc/o.20191021/multiqc_report.html
http://guru/new_home/qszhu/9.NT-ChIP/0.cutadapt-qc/3.cutadapt/o.20191021/multiqc_report.html


########################################################
#2020年6月14日
annoutil ls oss://annoroad-cloud-product/user/project/XS05KF2020050296/PM-XS05KF2020050296-01

annoutil cp oss://annoroad-cloud-product/user/project/XS05KF2020050296/PM-XS05KF2020050296-01/ANNO_XS05KF2020050296_PM-XS05KF2020050296-01_2020-06-11 ~/data/9.NT-ChIP/20200614 -r -f --jobs 3 --parallel 16 &
wdir=~/data/9.NT-ChIP/20200614
md5sum -c md5.txt


ddir=~/data/9.NT-ChIP/20200614/Rawdata
wdir=~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links/q.20200614
mkdir -p $wdir ${wdir/0.links/1.qc} ${wdir/0.links/3.cutadapt} ${wdir/0.links/4.qc_cutadapt}
#RNA
for i in 1 2; do 
ln -s $ddir/0519-xrm-1_R${i}.fq.gz $wdir/kdm4b-e2cell_RNA_rep1.R${i}.fastq.gz
ln -s $ddir/0519-xrm-2_R${i}.fq.gz $wdir/kdm4b-e2cell_RNA_rep2.R${i}.fastq.gz
ln -s $ddir/0519-xrm-3_R${i}.fq.gz $wdir/kdm4b-e2cell_RNA_rep3.R${i}.fastq.gz
ln -s $ddir/0519-xrm-4_R${i}.fq.gz $wdir/6h_RNA_rep4.R${i}.fastq.gz
ln -s $ddir/0519-xrm-5_R${i}.fq.gz $wdir/6h_RNA_rep5.R${i}.fastq.gz
ln -s $ddir/0519-xrm-6_R${i}.fq.gz $wdir/6h_RNA_rep6.R${i}.fastq.gz
ln -s $ddir/0519-xrm-7_R${i}.fq.gz $wdir/14h_RNA_rep4.R${i}.fastq.gz
ln -s $ddir/0519-xrm-8_R${i}.fq.gz $wdir/14h_RNA_rep5.R${i}.fastq.gz
ln -s $ddir/0519-xrm-9_R${i}.fq.gz $wdir/14h_RNA_rep6.R${i}.fastq.gz
ln -s $ddir/0519-xrm-10_R${i}.fq.gz $wdir/e2cell_RNA_rep4.R${i}.fastq.gz
ln -s $ddir/0519-xrm-11_R${i}.fq.gz $wdir/e2cell_RNA_rep5.R${i}.fastq.gz
ln -s $ddir/0519-xrm-12_R${i}.fq.gz $wdir/e2cell_RNA_rep6.R${i}.fastq.gz
ln -s $ddir/0519-xrm-13_R${i}.fq.gz $wdir/l2cell_RNA_rep4.R${i}.fastq.gz
ln -s $ddir/0519-xrm-14_R${i}.fq.gz $wdir/l2cell_RNA_rep5.R${i}.fastq.gz
ln -s $ddir/0519-xrm-15_R${i}.fq.gz $wdir/l2cell_RNA_rep6.R${i}.fastq.gz
ln -s $ddir/0519-xrm-16_R${i}.fq.gz $wdir/4cell_RNA_rep4.R${i}.fastq.gz
ln -s $ddir/0519-xrm-17_R${i}.fq.gz $wdir/4cell_RNA_rep5.R${i}.fastq.gz
ln -s $ddir/0519-xrm-18_R${i}.fq.gz $wdir/4cell_RNA_rep6.R${i}.fastq.gz
ln -s $ddir/0519-xrm-19_R${i}.fq.gz $wdir/8cell_RNA_rep4.R${i}.fastq.gz
ln -s $ddir/0519-xrm-20_R${i}.fq.gz $wdir/8cell_RNA_rep5.R${i}.fastq.gz
ln -s $ddir/0519-xrm-21_R${i}.fq.gz $wdir/8cell_RNA_rep6.R${i}.fastq.gz
ln -s $ddir/0519-xrm-22_R${i}.fq.gz $wdir/morula_RNA_rep4.R${i}.fastq.gz
ln -s $ddir/0519-xrm-23_R${i}.fq.gz $wdir/morula_RNA_rep5.R${i}.fastq.gz
ln -s $ddir/0519-xrm-24_R${i}.fq.gz $wdir/morula_RNA_rep6.R${i}.fastq.gz
ln -s $ddir/0519-xrm-25_R${i}.fq.gz $wdir/icm_RNA_rep4.R${i}.fastq.gz
ln -s $ddir/0519-xrm-26_R${i}.fq.gz $wdir/icm_RNA_rep5.R${i}.fastq.gz
ln -s $ddir/0519-xrm-27_R${i}.fq.gz $wdir/icm_RNA_rep6.R${i}.fastq.gz
ln -s $ddir/0519-xrm-28_R${i}.fq.gz $wdir/te_RNA_rep4.R${i}.fastq.gz
ln -s $ddir/0519-xrm-29_R${i}.fq.gz $wdir/te_RNA_rep5.R${i}.fastq.gz
ln -s $ddir/0519-xrm-30_R${i}.fq.gz $wdir/te_RNA_rep6.R${i}.fastq.gz
done

for i in 1 2; do 
ln -s $ddir/0504-xrm-7_R${i}.fq.gz $wdir/2cell_K9me3_rep5.R${i}.fastq.gz
ln -s $ddir/0504-xrm-8_R${i}.fq.gz $wdir/2cell_K9me3_rep6.R${i}.fastq.gz
ln -s $ddir/0504-xrm-5_R${i}.fq.gz $wdir/morula_K9me3_rep5.R${i}.fastq.gz
ln -s $ddir/0504-xrm-6_R${i}.fq.gz $wdir/morula_K9me3_rep6.R${i}.fastq.gz
ln -s $ddir/0504-xrm-1_R${i}.fq.gz $wdir/e2cell_K9me3_rep1.R${i}.fastq.gz
ln -s $ddir/0504-xrm-2_R${i}.fq.gz $wdir/e2cell_K9me3_rep2.R${i}.fastq.gz
ln -s $ddir/0504-xrm-3_R${i}.fq.gz $wdir/kdm4b-e2cell_K9me3_rep1.R${i}.fastq.gz
ln -s $ddir/0504-xrm-4_R${i}.fq.gz $wdir/kdm4b-e2cell_K9me3_rep2.R${i}.fastq.gz
done
nohup fastqc -o ${wdir/0.links/1.qc} $wdir/*gz -t 64 &

py35
for i in $wdir/*R1.fastq.gz; do o=${i/0.links/3.cutadapt}
nohup cutadapt -j 2 --trim-n -q 25,25 -m 30 -a ATCAACGCAGAGTAC -a AGATCGGAAGAGC -A ATCAACGCAGAGTAC -A AGATCGGAAGAGC -o $o -p ${o/R1/R2} $i ${i/R1/R2} > ${o/R1.fastq.gz/log} &
done
nohup fastqc -o ${wdir/0.links/4.qc_cutadapt} ${wdir/0.links/3.cutadapt}/*.gz -t 32 &

http://10.11.41.110/qszhu/9.NT-ChIP/0.cutadapt-qc/1.qc/q.20200614/multiqc_report.html



######################2020年8月17日
wget https://c.solargenomics.com/download/ossutil/linux/annoutil
annoutil config -e https://oss-cn-beijing.aliyuncs.com -i LTAI4GGEYjJ1FcNyF3rJi9LB -k f4eFcuTFrtRlTqd4U1u2yUKpqxnftZ
#有CM的数据一起
annoutil ls oss://annoroad-cloud-product/user/project/XS05KF2020070134/PM-XS05KF2020070134-07
annoutil cp oss://annoroad-cloud-product/user/project/XS05KF2020070134/PM-XS05KF2020070134-07/ANNO_XS05KF2020070134_PM-XS05KF2020070134-07_2020-08-16 ~/data/9.NT-ChIP/20200817 -r -f --jobs 3 --parallel 16 &

wdir=~/data/9.NT-ChIP/20200817
md5sum -c md5.txt

ddir=~/data/9.NT-ChIP/20200817/Rawdata
wdir=~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links/r.20200817
mkdir -p $wdir ${wdir/0.links/1.qc} ${wdir/0.links/3.cutadapt} ${wdir/0.links/4.qc_cutadapt}

for i in 1 2; do 
ln -s $ddir/*/0723XRM1_R${i}.fq.gz $wdir/control_RNA_rep1.R${i}.fastq.gz
ln -s $ddir/*/0723XRM2_R${i}.fq.gz $wdir/control_RNA_rep2.R${i}.fastq.gz
ln -s $ddir/*/0723XRM3_R${i}.fq.gz $wdir/control_RNA_rep3.R${i}.fastq.gz
ln -s $ddir/*/0723XRM4_R${i}.fq.gz $wdir/IL_RNA_rep1.R${i}.fastq.gz
ln -s $ddir/*/0723XRM5_R${i}.fq.gz $wdir/IL_RNA_rep2.R${i}.fastq.gz
ln -s $ddir/*/0723XRM6_R${i}.fq.gz $wdir/IL_RNA_rep3.R${i}.fastq.gz
ln -s $ddir/*/0723XRM7_R${i}.fq.gz $wdir/control_K9me3_rep1.R${i}.fastq.gz
ln -s $ddir/*/0723XRM8_R${i}.fq.gz $wdir/control_IgG_rep1.R${i}.fastq.gz
ln -s $ddir/*/0723XRM9_R${i}.fq.gz $wdir/IL_K9me3_rep1.R${i}.fastq.gz
ln -s $ddir/*/0723XRM10_R${i}.fq.gz $wdir/IL_IgG_rep1.R${i}.fastq.gz
done
nohup fastqc -o ${wdir/0.links/1.qc} $wdir/*gz -t 64 &

py35
for i in $wdir/*R1.fastq.gz; do o=${i/0.links/3.cutadapt}
nohup cutadapt -j 2 --trim-n -q 25,25 -m 30 -a ATCAACGCAGAGTAC -a AGATCGGAAGAGC -A ATCAACGCAGAGTAC -A AGATCGGAAGAGC -o $o -p ${o/R1/R2} $i ${i/R1/R2} > ${o/R1.fastq.gz/log} &
done
nohup fastqc -o ${wdir/0.links/4.qc_cutadapt} ${wdir/0.links/3.cutadapt}/*.gz -t 32 &

http://10.11.41.110/qszhu/9.NT-ChIP/0.cutadapt-qc/1.qc/r.20200817/multiqc_report.html



##################2020年9月20日 补IL RNA-seq
md5sum -c md5.txt
ddir=~/data/9.NT-ChIP/20200919
wdir=~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links/s.20200920

mkdir -p $wdir ${wdir/0.links/1.qc} ${wdir/0.links/3.cutadapt} ${wdir/0.links/4.qc_cutadapt}

for i in 1 2; do 
ln -s $ddir/B0903-xrm-1*R${i}.fastq.gz $wdir/control_RNA_rep4.R${i}.fastq.gz
ln -s $ddir/B0903-xrm-2*R${i}.fastq.gz $wdir/control_RNA_rep5.R${i}.fastq.gz
ln -s $ddir/B0903-xrm-3*R${i}.fastq.gz $wdir/control_RNA_rep6.R${i}.fastq.gz
ln -s $ddir/B0903-xrm-4*R${i}.fastq.gz $wdir/IL_RNA_rep4.R${i}.fastq.gz
ln -s $ddir/B0903-xrm-5*R${i}.fastq.gz $wdir/IL_RNA_rep5.R${i}.fastq.gz
ln -s $ddir/B0903-xrm-6*R${i}.fastq.gz $wdir/IL_RNA_rep6.R${i}.fastq.gz
done
nohup fastqc -o ${wdir/0.links/1.qc} $wdir/*gz -t 64 &


py35
for i in $wdir/*R1.fastq.gz; do o=${i/0.links/3.cutadapt}
nohup cutadapt -j 8 --trim-n -q 25,25 -m 30 -a ATCAACGCAGAGTAC -a AGATCGGAAGAGC -A ATCAACGCAGAGTAC -A AGATCGGAAGAGC -o $o -p ${o/R1/R2} $i ${i/R1/R2} > ${o/R1.fastq.gz/log} &
done
wait;de;de
nohup fastqc -o ${wdir/0.links/4.qc_cutadapt} ${wdir/0.links/3.cutadapt}/*.gz -t 32 &




#########2020年9月29日 6hpa sih2（siRNA-Suv39h2） 没有IgG
md5sum -c md5.txt
ddir=~/data/9.NT-ChIP/20200921
wdir=~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links/t.20200921

mkdir -p $wdir ${wdir/0.links/1.qc} ${wdir/0.links/3.cutadapt} ${wdir/0.links/4.qc_cutadapt}

for i in 1 2; do 
ln -s $ddir/0917-xrm-1*_${i}.fq.gz $wdir/6hpa_control_rep1.R${i}.fastq.gz
ln -s $ddir/0917-xrm-2*_${i}.fq.gz $wdir/6hpa_sih2_rep1.R${i}.fastq.gz
done
nohup fastqc -o ${wdir/0.links/1.qc} $wdir/*gz -t 64 &

py35
for i in $wdir/*R1.fastq.gz; do o=${i/0.links/3.cutadapt}
nohup cutadapt -j 8 --trim-n -q 25,25 -m 30 -a ATCAACGCAGAGTAC -a AGATCGGAAGAGC -A ATCAACGCAGAGTAC -A AGATCGGAAGAGC -o $o -p ${o/R1/R2} $i ${i/R1/R2} > ${o/R1.fastq.gz/log} &
done
wait;de;de
nohup fastqc -o ${wdir/0.links/4.qc_cutadapt} ${wdir/0.links/3.cutadapt}/*.gz -t 32 &


#2020年10月16日 加测数据
md5sum -c md5.txt
ddir=~/data/9.NT-ChIP/20201014
ddir2=~/data/9.NT-ChIP/20200921
wdir=~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links/t.20200921

mkdir -p $wdir ${wdir/0.links/1.qc} ${wdir/0.links/3.cutadapt} ${wdir/0.links/4.qc_cutadapt}
for i in 1 2; do 
zcat $ddir/0917-xrm-1*_${i}.fq.gz $ddir2/0917-xrm-1*_${i}.fq.gz | gzip -c - > $wdir/6hpa_control_rep1.R${i}.fastq.gz &
zcat $ddir/0917-xrm-2*_${i}.fq.gz $ddir2/0917-xrm-2*_${i}.fq.gz | gzip -c - > $wdir/6hpa_sih2_rep1.R${i}.fastq.gz &
done
nohup fastqc -o ${wdir/0.links/1.qc} $wdir/*gz -t 64 &


######################2020年11月24日 si Kat5 smart-seq
linuxnd login -u X101SC19050463-Z01-J092 -p p65g34f7
linuxnd cp -d oss://CP2019051000025/H101SC19050463/RSCS0500/X101SC19050463-Z01/X101SC19050463-Z01-J092/1.rawdata/ ~/data/9.NT-ChIP/20201124

ddir=~/data/9.NT-ChIP/20201124
wdir=~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links/u.20201124
mkdir -p $wdir ${wdir/0.links/1.qc} ${wdir/0.links/3.cutadapt} ${wdir/0.links/4.qc_cutadapt}


for i in 1 2; do  
ln -s $ddir/1118-XRM-1*_${i}.fq.gz $wdir/morula_control_RNA_rep1.R${i}.fastq.gz; 
ln -s $ddir/1118-XRM-2*_${i}.fq.gz $wdir/morula_sikat5_RNA_rep1.R${i}.fastq.gz; 
ln -s $ddir/1118-XRM-3*_${i}.fq.gz $wdir/ICM_control_RNA_rep1.R${i}.fastq.gz; 
ln -s $ddir/1118-XRM-4*_${i}.fq.gz $wdir/TE_control_RNA_rep1.R${i}.fastq.gz; 
ln -s $ddir/1118-XRM-5*_${i}.fq.gz $wdir/ICM_sikat5_RNA_rep1.R${i}.fastq.gz; 
ln -s $ddir/1118-XRM-6*_${i}.fq.gz $wdir/TE_sikat5_RNA_rep1.R${i}.fastq.gz; done


nohup fastqc -o ${wdir/0.links/1.qc} $wdir/*gz -t 64 &

py35
for i in $wdir/*R1.fastq.gz; do o=${i/0.links/3.cutadapt}
nohup cutadapt -j 8 --trim-n -q 25,25 -m 30 -a ATCAACGCAGAGTAC -a AGATCGGAAGAGC -A ATCAACGCAGAGTAC -A AGATCGGAAGAGC -o $o -p ${o/R1/R2} $i ${i/R1/R2} > ${o/R1.fastq.gz/log} &
done
wait;de;de
nohup fastqc -o ${wdir/0.links/4.qc_cutadapt} ${wdir/0.links/3.cutadapt}/*.gz -t 32 &




####################2020年12月7日 si Kat5 K9me3 cut&run
linuxnd login -u X101SC19050463-Z01-J094 -p fbrf02tn
linuxnd cp -d oss://CP2019051000025/H101SC19050463/RSCS0500/X101SC19050463-Z01/X101SC19050463-Z01-J094/1.rawdata/ ~/data/9.NT-ChIP/20201207
mv */* ./ 
cat *txt > md5.txt 
md5sum -c md5.txt
ddir=~/data/9.NT-ChIP/20201207
wdir=~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links/v.20201207
mkdir -p $wdir ${wdir/0.links/1.qc} ${wdir/0.links/3.cutadapt} ${wdir/0.links/4.qc_cutadapt}

for i in 1 2; do  
ln -s $ddir/XRM-1122-1*_${i}.fq.gz $wdir/morula_control_K9me3_rep1.R${i}.fastq.gz; 
ln -s $ddir/XRM-1122-2*_${i}.fq.gz $wdir/morula_control_K9me3_rep2.R${i}.fastq.gz; 
ln -s $ddir/XRM-1122-3*_${i}.fq.gz $wdir/morula_sikat5_K9me3_rep1.R${i}.fastq.gz; 
ln -s $ddir/XRM-1122-4*_${i}.fq.gz $wdir/morula_sikat5_K9me3_rep2.R${i}.fastq.gz;

ln -s $ddir/XRM-1122-5*_${i}.fq.gz $wdir/ICM_control_K9me3_rep1.R${i}.fastq.gz; 
ln -s $ddir/XRM-1122-6*_${i}.fq.gz $wdir/TE_control_K9me3_rep1.R${i}.fastq.gz; 
ln -s $ddir/XRM-1122-7*_${i}.fq.gz $wdir/ICM_sikat5_K9me3_rep1.R${i}.fastq.gz; 
done





#############2020年12月27日 siKat5 smart-seq 
#原始数据已迁移到 siMIT.NF.smart
ddir=~/data/9.NT-ChIP/20201224
ddir2=~/data/9.NT-ChIP/20201227
wdir=~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links/w.20201227
mkdir -p $wdir ${wdir/0.links/1.qc} ${wdir/0.links/3.cutadapt} ${wdir/0.links/4.qc_cutadapt}

for i in 1 2; do  
ln -s $ddir/XRM-1129-1*R${i}.fastq.gz $wdir/morula_control_RNA_rep3.R${i}.fastq.gz; 
ln -s $ddir/XRM-1129-2*R${i}.fastq.gz $wdir/morula_sikat5_RNA_rep2.R${i}.fastq.gz; 
ln -s $ddir/XRM-1129-3*R${i}.fastq.gz $wdir/ICM_control_RNA_rep2.R${i}.fastq.gz; 
ln -s $ddir/XRM-1129-4*R${i}.fastq.gz $wdir/TE_control_RNA_rep2.R${i}.fastq.gz; 
ln -s $ddir/XRM-1129-5*R${i}.fastq.gz $wdir/TE_sikat5_RNA_rep2.R${i}.fastq.gz; 
ln -s $ddir/XRM-1129-6*R${i}.fastq.gz $wdir/ICM_control_RNA_rep3.R${i}.fastq.gz; 
ln -s $ddir/XRM-1129-7*R${i}.fastq.gz $wdir/ICM_sikat5_RNA_rep3.R${i}.fastq.gz; 
ln -s $ddir/XRM-1129-8*R${i}.fastq.gz $wdir/TE_control_RNA_rep3.R${i}.fastq.gz; 
ln -s $ddir/XRM-1129-9*R${i}.fastq.gz $wdir/TE_sikat5_RNA_rep3.R${i}.fastq.gz; 
ln -s $ddir2/XRM-1210-1_*R${i}.fastq.gz $wdir/morula_control_RNA_rep2.R${i}.fastq.gz; 
ln -s $ddir2/XRM-1210-5*R${i}.fastq.gz $wdir/TE_control_RNA_rep4.R${i}.fastq.gz; 
ln -s $ddir2/XRM-1210-6*R${i}.fastq.gz $wdir/ICM_siKat5_RNA_rep4.R${i}.fastq.gz; 
ln -s $ddir2/XRM-1210-2*R${i}.fastq.gz $wdir/morula_siMcrs1_RNA_rep1.R${i}.fastq.gz; 
ln -s $ddir2/XRM-1210-15*R${i}.fastq.gz $wdir/TE_siMcrs1_RNA_rep2.R${i}.fastq.gz; 
done


nohup fastqc -o ${wdir/0.links/1.qc} $wdir/*gz -t 64 &
py35
for i in $wdir/*R1.fastq.gz; do o=${i/0.links/3.cutadapt}
nohup cutadapt -j 8 --trim-n -q 25,25 -m 30 -a ATCAACGCAGAGTAC -a AGATCGGAAGAGC -A ATCAACGCAGAGTAC -A AGATCGGAAGAGC -o $o -p ${o/R1/R2} $i ${i/R1/R2} > ${o/R1.fastq.gz/log} &
done
wait;de;de
nohup fastqc -o ${wdir/0.links/4.qc_cutadapt} ${wdir/0.links/3.cutadapt}/*.gz -t 32 &
wait
cd ${wdir/0.links/1.qc};multiqc --no-data-dir -f * -dd 10; mv multiqc_report.html ${wdir}/raw.fastqc.html 
cd ${wdir/0.links/4.qc_cutadapt};multiqc --no-data-dir -f * -dd 10; mv multiqc_report.html ${wdir}/cutadapt.fastqc.html 



#漏掉了一些smart-seq记录
##########################################
#ICM TE都是E4
#所有Morula ICM TE的siKat5 siMax siMcrs1和control的cut&run K9me3数据合并处理
#20210110
#20210125
#20201207
mv ... > siMIT.K9.cutrun
ddir=~/data/9.NT-ChIP/siMIT.K9.cutrun
wdir=~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links/7.siMIT-K9
mkdir -p $wdir ${wdir/0.links/1.qc} ${wdir/0.links/3.cutadapt} ${wdir/0.links/4.qc_cutadapt}


for i in 1 2; do  
ln -s $ddir/XRM-1122-1*_${i}.fq.gz $wdir/Morula_control_K9me3_rep1.R${i}.fastq.gz; 
ln -s $ddir/XRM-1122-2*_${i}.fq.gz $wdir/Morula_control_K9me3_rep2.R${i}.fastq.gz; 
ln -s $ddir/XRM-1122-3*_${i}.fq.gz $wdir/Morula_siKat5_K9me3_rep1.R${i}.fastq.gz; 
ln -s $ddir/XRM-1122-4*_${i}.fq.gz $wdir/Morula_siKat5_K9me3_rep2.R${i}.fastq.gz;
ln -s $ddir/XRM-1122-5*_${i}.fq.gz $wdir/ICM_control_K9me3_rep1.R${i}.fastq.gz; 
ln -s $ddir/XRM-1122-6*_${i}.fq.gz $wdir/TE_control_K9me3_rep1.R${i}.fastq.gz; 
ln -s $ddir/XRM-1122-7*_${i}.fq.gz $wdir/ICM_siKat5_K9me3_rep1.R${i}.fastq.gz; 

ln -s $ddir/XRM-1211-10_L3_A011.R${i}.fastq.gz $wdir/ICM_siMcrs1_K9me3_rep1.R${i}.fastq.gz;
ln -s $ddir/XRM-1211-11_L3_A012.R${i}.fastq.gz $wdir/TE_siMcrs1_K9me3_rep1.R${i}.fastq.gz;
ln -s $ddir/XRM-1211-12_L3_A008.R${i}.fastq.gz $wdir/TE_siMcrs1_K9me3_rep2.R${i}.fastq.gz;
ln -s $ddir/XRM-1211-1_L3_A001.R${i}.fastq.gz $wdir/ICM_control_K9me3_rep2.R${i}.fastq.gz;
ln -s $ddir/XRM-1211-2_L3_A002.R${i}.fastq.gz $wdir/TE_control_K9me3_rep2.R${i}.fastq.gz;
ln -s $ddir/XRM-1211-3_L3_A003.R${i}.fastq.gz $wdir/TE_control_K9me3_rep3.R${i}.fastq.gz;
ln -s $ddir/XRM-1211-4_L3_A004.R${i}.fastq.gz $wdir/ICM_siKat5_K9me3_rep2.R${i}.fastq.gz;
ln -s $ddir/XRM-1211-5_L3_A005.R${i}.fastq.gz $wdir/TE_siKat5_K9me3_rep2.R${i}.fastq.gz;
ln -s $ddir/XRM-1211-6_L3_A006.R${i}.fastq.gz $wdir/TE_siKat5_K9me3_rep3.R${i}.fastq.gz;
ln -s $ddir/XRM-1211-7_L3_A007.R${i}.fastq.gz $wdir/ICM_siMax_K9me3_rep1.R${i}.fastq.gz;
ln -s $ddir/XRM-1211-8_L3_A009.R${i}.fastq.gz $wdir/TE_siMax_K9me3_rep1.R${i}.fastq.gz;
ln -s $ddir/XRM-1211-9_L3_A010.R${i}.fastq.gz $wdir/TE_siMax_K9me3_rep2.R${i}.fastq.gz;

ln -s $ddir/XRM-1230-10_L3_A029.R${i}.fastq.gz $wdir/TE_siMcrs1_K9me3_rep3.R${i}.fastq.gz;
ln -s $ddir/XRM-1230-11_L3_A030.R${i}.fastq.gz $wdir/TE_siMcrs1_K9me3_rep4.R${i}.fastq.gz;
ln -s $ddir/XRM-1230-1_L1_A010.R${i}.fastq.gz $wdir/ICM_control_K9me3_rep3.R${i}.fastq.gz;
ln -s $ddir/XRM-1230-2_L3_A011.R${i}.fastq.gz $wdir/TE_control_K9me3_rep4.R${i}.fastq.gz;
ln -s $ddir/XRM-1230-3_L3_A021.R${i}.fastq.gz $wdir/TE_control_K9me3_rep5.R${i}.fastq.gz;
ln -s $ddir/XRM-1230-4_L3_A022.R${i}.fastq.gz $wdir/ICM_siKat5_K9me3_rep3.R${i}.fastq.gz;
ln -s $ddir/XRM-1230-5_L3_A023.R${i}.fastq.gz $wdir/TE_siKat5_K9me3_rep4.R${i}.fastq.gz;
ln -s $ddir/XRM-1230-6_L3_A024.R${i}.fastq.gz $wdir/TE_siKat5_K9me3_rep5.R${i}.fastq.gz;
ln -s $ddir/XRM-1230-7_L4_A026.R${i}.fastq.gz $wdir/ICM_siMax_K9me3_rep2.R${i}.fastq.gz;
ln -s $ddir/XRM-1230-8_L3_A027.R${i}.fastq.gz $wdir/TE_siMax_K9me3_rep3.R${i}.fastq.gz;
ln -s $ddir/XRM-1230-9_L3_A028.R${i}.fastq.gz $wdir/TE_siMax_K9me3_rep4.R${i}.fastq.gz;
done

#检验
ll *R1*gz | awk 'BEGIN{FS="[ /.]"} {split($8,temp,"_K9me3_");print temp[1],temp[2],$22}'


nohup fastqc -o ${wdir/0.links/1.qc} $wdir/*gz -t 64 &
for i in $wdir/*R1.fastq.gz; do o=${i/0.links/3.cutadapt}
nohup cutadapt -j 8 --trim-n -q 25,25 -m 30 -a ATCAACGCAGAGTAC -a AGATCGGAAGAGC -A ATCAACGCAGAGTAC -A AGATCGGAAGAGC -o $o -p ${o/R1/R2} $i ${i/R1/R2} > ${o/R1.fastq.gz/log} &
done
wait
nohup fastqc -o ${wdir/0.links/4.qc_cutadapt} ${wdir/0.links/3.cutadapt}/*.gz -t 32 &

#
cd ${wdir/0.links/1.qc} 
multiqc --no-data-dir -f -dd 10 ./ -n $wdir/raw.qc.html
cd ${wdir/0.links/3.cutadapt} 
multiqc --no-data-dir -f -dd 10 ./ -n $wdir/cutadapt.html
cd ${wdir/0.links/4.qc_cutadapt}
multiqc --no-data-dir -f -dd 10 ./ -n $wdir/clean.qc.html


#################################################################
#2021年2月23日 siMIT Kat5 Mcrs1 Max NT RNA+cut&run
awk -v FS="," 'NR>1{split($2,temp,"-");if(temp[2]=="C"){temp[2]="control"}else{temp[2]="si"temp[2]}
print "ln -s $ddir/"$1"_*.R${i}.fastq.gz $wdir/"temp[1]"-"temp[3]"_"temp[2]"_RNA_rep"temp[4]".R${i}.fastq.gz;"
}' sample-info.smart.info
awk -v FS="," 'NR>1{split($2,temp,"-");if(temp[2]=="C"){temp[2]="control"}else{temp[2]="si"temp[2]}
print "ln -s $ddir/"$1"_*.R${i}.fastq.gz $wdir/"temp[1]"-"temp[3]"_"temp[2]"_K9me3_rep"temp[4]".R${i}.fastq.gz;"
}' sample-info.cutrun.info

ddir=~/data/9.NT-ChIP/siMIT.NT
wdir=~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links/8.siMIT-NT
mkdir -p $wdir ${wdir/0.links/1.qc} ${wdir/0.links/3.cutadapt} ${wdir/0.links/4.qc_cutadapt}

for i in 1 2; do  
ln -s $ddir/XRM-0114-1_*.R${i}.fastq.gz $wdir/NT-ICM_control_RNA_rep1.R${i}.fastq.gz;
ln -s $ddir/XRM-0114-2_*.R${i}.fastq.gz $wdir/NT-ICM_control_RNA_rep2.R${i}.fastq.gz;
ln -s $ddir/XRM-0114-3_*.R${i}.fastq.gz $wdir/NT-TE_control_RNA_rep1.R${i}.fastq.gz;
ln -s $ddir/XRM-0114-4_*.R${i}.fastq.gz $wdir/NT-TE_control_RNA_rep2.R${i}.fastq.gz;
ln -s $ddir/XRM-0114-5_*.R${i}.fastq.gz $wdir/NT-ICM_oeKat5_RNA_rep1.R${i}.fastq.gz;
ln -s $ddir/XRM-0114-6_*.R${i}.fastq.gz $wdir/NT-ICM_oeKat5_RNA_rep2.R${i}.fastq.gz;
ln -s $ddir/XRM-0114-7_*.R${i}.fastq.gz $wdir/NT-TE_oeKat5_RNA_rep1.R${i}.fastq.gz;
ln -s $ddir/XRM-0114-8_*.R${i}.fastq.gz $wdir/NT-TE_oeKat5_RNA_rep2.R${i}.fastq.gz;
ln -s $ddir/XRM-0114-9_*.R${i}.fastq.gz $wdir/NT-ICM_oeMax_RNA_rep1.R${i}.fastq.gz;
ln -s $ddir/XRM-0114-10_*.R${i}.fastq.gz $wdir/NT-ICM_oeMax_RNA_rep2.R${i}.fastq.gz;
ln -s $ddir/XRM-0114-11_*.R${i}.fastq.gz $wdir/NT-TE_oeMax_RNA_rep1.R${i}.fastq.gz;
ln -s $ddir/XRM-0114-12_*.R${i}.fastq.gz $wdir/NT-TE_oeMax_RNA_rep2.R${i}.fastq.gz;
ln -s $ddir/XRM-0114-13_*.R${i}.fastq.gz $wdir/NT-ICM_oeMcrs1_RNA_rep1.R${i}.fastq.gz;
ln -s $ddir/XRM-0114-14_*.R${i}.fastq.gz $wdir/NT-ICM_oeMcrs1_RNA_rep2.R${i}.fastq.gz;
ln -s $ddir/XRM-0114-15_*.R${i}.fastq.gz $wdir/NT-TE_oeMcrs1_RNA_rep1.R${i}.fastq.gz;
ln -s $ddir/XRM-0114-16_*.R${i}.fastq.gz $wdir/NT-TE_oeMcrs1_RNA_rep2.R${i}.fastq.gz;
ln -s $ddir/XRM-0125-1_*.R${i}.fastq.gz $wdir/NT-TE_control_RNA_rep3.R${i}.fastq.gz;
ln -s $ddir/XRM-0125-2_*.R${i}.fastq.gz $wdir/NT-TE_control_RNA_rep4.R${i}.fastq.gz;
ln -s $ddir/XRM-0125-3_*.R${i}.fastq.gz $wdir/NT-ICM_oeMax_RNA_rep3.R${i}.fastq.gz;
ln -s $ddir/XRM-0125-4_*.R${i}.fastq.gz $wdir/NT-TE_oeMax_RNA_rep3.R${i}.fastq.gz;
ln -s $ddir/XRM-0125-5_*.R${i}.fastq.gz $wdir/NT-TE_oeMax_RNA_rep4.R${i}.fastq.gz;

ln -s $ddir/XRM-0111-1_*.R${i}.fastq.gz $wdir/NT-ICM_control_K9me3_rep1.R${i}.fastq.gz;
ln -s $ddir/XRM-0111-2_*.R${i}.fastq.gz $wdir/NT-ICM_control_K9me3_rep2.R${i}.fastq.gz;
ln -s $ddir/XRM-0111-3_*.R${i}.fastq.gz $wdir/NT-TE_control_K9me3_rep1.R${i}.fastq.gz;
ln -s $ddir/XRM-0111-4_*.R${i}.fastq.gz $wdir/NT-TE_control_K9me3_rep2.R${i}.fastq.gz;
ln -s $ddir/XRM-0111-5_*.R${i}.fastq.gz $wdir/NT-TE_control_K9me3_rep3.R${i}.fastq.gz;
ln -s $ddir/XRM-0111-6_*.R${i}.fastq.gz $wdir/NT-TE_control_K9me3_rep4.R${i}.fastq.gz;
ln -s $ddir/XRM-0111-7_*.R${i}.fastq.gz $wdir/NT-ICM_oeKat5_K9me3_rep1.R${i}.fastq.gz;
ln -s $ddir/XRM-0111-8_*.R${i}.fastq.gz $wdir/NT-TE_oeKat5_K9me3_rep1.R${i}.fastq.gz;
ln -s $ddir/XRM-0111-9_*.R${i}.fastq.gz $wdir/NT-TE_oeKat5_K9me3_rep2.R${i}.fastq.gz;
ln -s $ddir/XRM-0111-10_*.R${i}.fastq.gz $wdir/NT-ICM_oeMax_K9me3_rep1.R${i}.fastq.gz;
ln -s $ddir/XRM-0111-11_*.R${i}.fastq.gz $wdir/NT-TE_oeMax_K9me3_rep1.R${i}.fastq.gz;
ln -s $ddir/XRM-0111-12_*.R${i}.fastq.gz $wdir/NT-TE_oeMax_K9me3_rep2.R${i}.fastq.gz;
ln -s $ddir/XRM-0124-8_*.R${i}.fastq.gz $wdir/NT-ICM_control_K9me3_rep3.R${i}.fastq.gz;
ln -s $ddir/XRM-0124-9_*.R${i}.fastq.gz $wdir/NT-TE_control_K9me3_rep5.R${i}.fastq.gz;
ln -s $ddir/XRM-0124-10_*.R${i}.fastq.gz $wdir/NT-TE_control_K9me3_rep6.R${i}.fastq.gz;
ln -s $ddir/XRM-0124-11_*.R${i}.fastq.gz $wdir/NT-ICM_oeMax_K9me3_rep2.R${i}.fastq.gz;
ln -s $ddir/XRM-0124-12_*.R${i}.fastq.gz $wdir/NT-TE_oeMax_K9me3_rep3.R${i}.fastq.gz;
ln -s $ddir/XRM-0124-13_*.R${i}.fastq.gz $wdir/NT-TE_oeMax_K9me3_rep4.R${i}.fastq.gz;
done
nohup fastqc -o ${wdir/0.links/1.qc} $wdir/*gz -t 64 &
for i in $wdir/*R1.fastq.gz; do o=${i/0.links/3.cutadapt}
nohup cutadapt -j 8 --trim-n -q 25,25 -m 30 -a ATCAACGCAGAGTAC -a AGATCGGAAGAGC -A ATCAACGCAGAGTAC -A AGATCGGAAGAGC -o $o -p ${o/R1/R2} $i ${i/R1/R2} > ${o/R1.fastq.gz/log} &
done
wait
nohup fastqc -o ${wdir/0.links/4.qc_cutadapt} ${wdir/0.links/3.cutadapt}/*.gz -t 32 &

ddir=~/data/9.NT-ChIP/siMIT.NT
wdir=~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links/8.siMIT-NT
rm $wdir/*gz ${wdir/0.links/3.cutadapt}/*gz
for i in 1 2; do 
ln -s $ddir/B0226-XRM-1_*.R${i}.fastq.gz $wdir/NT-ICM_control_RNA_rep3.R${i}.fastq.gz
ln -s $ddir/B0226-XRM-2_*.R${i}.fastq.gz $wdir/NT-ICM_oeMax_RNA_rep4.R${i}.fastq.gz
ln -s $ddir/B0226-XRM-3_*.R${i}.fastq.gz $wdir/NT-ICM_oeMax_RNA_rep5.R${i}.fastq.gz
ln -s $ddir/B0226-XRM-4_*.R${i}.fastq.gz $wdir/NT-ICM_oeMcrs1_RNA_rep3.R${i}.fastq.gz
ln -s $ddir/B0226-XRM-5_*.R${i}.fastq.gz $wdir/NT-TE_control_RNA_rep5.R${i}.fastq.gz
ln -s $ddir/B0226-XRM-6_*.R${i}.fastq.gz $wdir/NT-TE_oeMax_RNA_rep5.R${i}.fastq.gz
ln -s $ddir/B0226-XRM-7_*.R${i}.fastq.gz $wdir/NT-TE_oeMax_RNA_rep6.R${i}.fastq.gz
ln -s $ddir/B0226-XRM-8_*.R${i}.fastq.gz $wdir/NT-TE_oeMcrs1_RNA_rep3.R${i}.fastq.gz

ln -s $ddir/B0224-XRM-1_*.R${i}.fastq.gz $wdir/NT-ICM_control_K9me3_rep4.R${i}.fastq.gz
ln -s $ddir/B0224-XRM-2_*.R${i}.fastq.gz $wdir/NT-ICM_oeMax_K9me3_rep3.R${i}.fastq.gz
ln -s $ddir/B0224-XRM-3_*.R${i}.fastq.gz $wdir/NT-ICM_oeMcrs1_K9me3_rep1.R${i}.fastq.gz
ln -s $ddir/B0224-XRM-4_*.R${i}.fastq.gz $wdir/NT-TE_control_K9me3_rep7.R${i}.fastq.gz
ln -s $ddir/B0224-XRM-5_*.R${i}.fastq.gz $wdir/NT-TE_oeMax_K9me3_rep5.R${i}.fastq.gz
ln -s $ddir/B0224-XRM-6_*.R${i}.fastq.gz $wdir/NT-TE_oeMcrs1_K9me3_rep1.R${i}.fastq.gz
done


ddir=~/data/9.NT-ChIP/siMIT.NF.smart
wdir=~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links/9.siMIT-mph
mkdir -p $wdir ${wdir/0.links/1.qc} ${wdir/0.links/3.cutadapt} ${wdir/0.links/4.qc_cutadapt}
for i in 1 2; do 
ln -s $ddir/XRM-0313-1_*.R${i}.fastq.gz $wdir/NF-ICM_control_mph_rep1.R${i}.fastq.gz
ln -s $ddir/XRM-0313-2_*.R${i}.fastq.gz $wdir/NF-ICM_control_mph_rep2.R${i}.fastq.gz
ln -s $ddir/XRM-0313-3_*.R${i}.fastq.gz $wdir/NF-TE_control_mph_rep1.R${i}.fastq.gz
ln -s $ddir/XRM-0313-4_*.R${i}.fastq.gz $wdir/NF-TE_control_mph_rep2.R${i}.fastq.gz
ln -s $ddir/XRM-0313-5_*.R${i}.fastq.gz $wdir/NF-TE_siMax_mph_rep1.R${i}.fastq.gz
ln -s $ddir/XRM-0313-6_*.R${i}.fastq.gz $wdir/NF-TE_siMax_mph_rep2.R${i}.fastq.gz
done

###################
ddir=~/data/9.NT-ChIP/siMIT.NT/20210420
wdir=~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links/8.siMIT-NT/20210421
mkdir -p $wdir ${wdir/0.links/1.qc} ${wdir/0.links/3.cutadapt} ${wdir/0.links/4.qc_cutadapt}
for i in 1 2; do 
ln -s $ddir/XRM-0409-4*_$i.fq.gz $wdir/NT-ICM_oeKat5_K9me3_rep2.R${i}.fastq.gz; 
ln -s $ddir/XRM-0409-5*_$i.fq.gz $wdir/NT-ICM_oeKat5_IgG_rep1.R${i}.fastq.gz; 
ln -s $ddir/XRM-0409-6*_$i.fq.gz $wdir/NT-ICM_oeMax_IgG_rep1.R${i}.fastq.gz; 
ln -s $ddir/XRM-0409-7*_$i.fq.gz $wdir/NT-TE_oeKat5_IgG_rep1.R${i}.fastq.gz; 
ln -s $ddir/XRM-0409-8*_$i.fq.gz $wdir/NT-TE_oeMax_IgG_rep1.R${i}.fastq.gz; 
done
mv *gz ../


#####################
ddir=~/data/9.NT-ChIP/20210423
wdir=~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links/y.20210425
mkdir -p $wdir ${wdir/0.links/1.qc} ${wdir/0.links/3.cutadapt} ${wdir/0.links/4.qc_cutadapt}
#awk 'NR>1{print "ln -s $ddir/"$1"_*.R$i.fastq.gz $wdir/"$2".R$i.fastq.gz"}' sample-information.tab
for i in 1 2; do 
ln -s $ddir/XRM-0323-1_*.R$i.fastq.gz $wdir/NF-ICM_control_K9me3_rep2.R$i.fastq.gz
ln -s $ddir/XRM-0323-2_*.R$i.fastq.gz $wdir/NF-ICM_control_IgG_rep1.R$i.fastq.gz
ln -s $ddir/XRM-0323-3_*.R$i.fastq.gz $wdir/NF-ICM_mphMcrs1_K9me3_rep1.R$i.fastq.gz
ln -s $ddir/XRM-0323-4_*.R$i.fastq.gz $wdir/NF-TE_control_K9me3_rep2.R$i.fastq.gz
ln -s $ddir/XRM-0323-5_*.R$i.fastq.gz $wdir/NF-TE_control_IgG_rep1.R$i.fastq.gz
ln -s $ddir/XRM-0323-6_*.R$i.fastq.gz $wdir/NF-TE_mphMcrs1_K9me3_rep1.R$i.fastq.gz
ln -s $ddir/XRM-0323-7_*.R$i.fastq.gz $wdir/NF-TE_mphMcrs1_IgG_rep1.R$i.fastq.gz
ln -s $ddir/XRM-0323-8_*.R$i.fastq.gz $wdir/NT-ICM_control_K9me3_rep5.R$i.fastq.gz
ln -s $ddir/XRM-0323-9_*.R$i.fastq.gz $wdir/NT-ICM_control_IgG_rep1.R$i.fastq.gz
ln -s $ddir/XRM-0323-10_*.R$i.fastq.gz $wdir/NT-ICM_oeMcrs1_K9me3_rep2.R$i.fastq.gz
ln -s $ddir/XRM-0323-11_*.R$i.fastq.gz $wdir/NT-ICM_oeMcrs1_IgG_rep1.R$i.fastq.gz
ln -s $ddir/XRM-0323-12_*.R$i.fastq.gz $wdir/NT-TE_control_K9me3_rep8.R$i.fastq.gz
ln -s $ddir/XRM-0323-13_*.R$i.fastq.gz $wdir/NT-TE_control_IgG_rep1.R$i.fastq.gz
ln -s $ddir/XRM-0323-14_*.R$i.fastq.gz $wdir/NT-TE_oeMcrs1_K9me3_rep2.R$i.fastq.gz
ln -s $ddir/XRM-0323-15_*.R$i.fastq.gz $wdir/NT-TE_oeMcrs1_IgG_rep1.R$i.fastq.gz
ln -s $ddir/XRM-0325-1_*.R$i.fastq.gz $wdir/NF-ICM_control_RNA_rep3.R$i.fastq.gz
ln -s $ddir/XRM-0325-2_*.R$i.fastq.gz $wdir/NF-ICM_control_RNA_rep4.R$i.fastq.gz
ln -s $ddir/XRM-0325-3_*.R$i.fastq.gz $wdir/NF-ICM_mphMcrs1_RNA_rep1.R$i.fastq.gz
ln -s $ddir/XRM-0325-4_*.R$i.fastq.gz $wdir/NF-ICM_mphMcrs1_RNA_rep2.R$i.fastq.gz
ln -s $ddir/XRM-0325-5_*.R$i.fastq.gz $wdir/NF-TE_control_RNA_rep3.R$i.fastq.gz
ln -s $ddir/XRM-0325-6_*.R$i.fastq.gz $wdir/NF-TE_control_RNA_rep4.R$i.fastq.gz
ln -s $ddir/XRM-0325-7_*.R$i.fastq.gz $wdir/NF-TE_mphMcrs1_RNA_rep1.R$i.fastq.gz
ln -s $ddir/XRM-0325-8_*.R$i.fastq.gz $wdir/NF-TE_mphMcrs1_RNA_rep2.R$i.fastq.gz
ln -s $ddir/XRM-0325-9_*.R$i.fastq.gz $wdir/NT-TE_oeMcrs1_RNA_rep4.R$i.fastq.gz
ln -s $ddir/XRM-0325-10_*.R$i.fastq.gz $wdir/NT-TE_oeMcrs1_RNA_rep5.R$i.fastq.gz
ln -s $ddir/XRM-0330-1_*.R$i.fastq.gz $wdir/NF-ICM_control_K9me3_rep3.R$i.fastq.gz
ln -s $ddir/XRM-0330-2_*.R$i.fastq.gz $wdir/NF-ICM_mphMcrs1_K9me3_rep2.R$i.fastq.gz
ln -s $ddir/XRM-0330-3_*.R$i.fastq.gz $wdir/NF-TE_mphMax_K9me3_rep2.R$i.fastq.gz
ln -s $ddir/XRM-0330-4_*.R$i.fastq.gz $wdir/NF-TE_mphMax_IgG_rep1.R$i.fastq.gz
ln -s $ddir/XRM-0330-5_*.R$i.fastq.gz $wdir/NF-TE_mphMcrs1_K9me3_rep2.R$i.fastq.gz
ln -s $ddir/XRM-0331-1_*.R$i.fastq.gz $wdir/NF-ICM_mphMax_RNA_rep1.R$i.fastq.gz
ln -s $ddir/XRM-0331-2_*.R$i.fastq.gz $wdir/NF-ICM_mphMax_RNA_rep2.R$i.fastq.gz
done

nohup fastqc -o ${wdir/0.links/1.qc} $wdir/*gz -t 64 &
for i in $wdir/*R1.fastq.gz; do o=${i/0.links/3.cutadapt}
nohup cutadapt -j 8 --trim-n -q 25,25 -m 30 -a ATCAACGCAGAGTAC -a AGATCGGAAGAGC -A ATCAACGCAGAGTAC -A AGATCGGAAGAGC -o $o -p ${o/R1/R2} $i ${i/R1/R2} > ${o/R1.fastq.gz/log} &
done
wait
nohup fastqc -o ${wdir/0.links/4.qc_cutadapt} ${wdir/0.links/3.cutadapt}/*.gz -t 32 &
wait
cd ${wdir/0.links/1.qc}
multiqc --no-data-dir -f -dd 10 ./ -n $wdir/raw.qc.html
cd ${wdir/0.links/3.cutadapt}
multiqc --no-data-dir -f -dd 10 ./ -n $wdir/cutadapt.html
cd ${wdir/0.links/4.qc_cutadapt}
multiqc --no-data-dir -f -dd 10 ./ -n $wdir/clean.qc.html

#############################################################################
linuxnd login -u X101SC19050463-Z01-J118 -p 9bxapmg0
linuxnd cp -d oss://CP2019051000025/H101SC19050463/RSCS0500/X101SC19050463-Z01/X101SC19050463-Z01-J118/1.rawdata/ ~/data/9.NT-ChIP/siMIT.20210507

ddir=~/data/9.NT-ChIP/siMIT.20210507
wdir=~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links/z.20210507
mkdir -p $wdir ${wdir/0.links/1.qc} ${wdir/0.links/3.cutadapt} ${wdir/0.links/4.qc_cutadapt}

for i in 1 2; do 
ln -fs $ddir/XRM-0426-1_*_$i.fq.gz $wdir/NF-ICM_mphMax_K9me3_rep2.R$i.fastq.gz
ln -fs $ddir/XRM-0426-2_*_$i.fq.gz $wdir/NF-ICM_mphMax_IgG_rep1.R$i.fastq.gz
ln -fs $ddir/XRM-0426-3_*_$i.fq.gz $wdir/NF-ICM_mphMcrs1_IgG_rep1.R$i.fastq.gz
done 




#####################
#2021年10月9日 kdm4d NT  L1C L2C NChIP K9me3
ddir=/home/share/DATA/clean_reads/GaoSR_mouse_nuclear_transfer_ChIP_20180301/20211009-4b
wdir=~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links/za.20211009
mkdir -p $wdir ${wdir/0.links/1.qc} ${wdir/0.links/3.cutadapt} ${wdir/0.links/4.qc_cutadapt}

for i in 1 2; do 
ln -s $ddir/XRM-0812-1_*_$i.fq.gz $wdir/NT-14h_kdm4b_input_rep1.R$i.fastq.gz
ln -s $ddir/XRM-0812-2_*_$i.fq.gz $wdir/NT-14h_kdm4b_K9me3_rep1.R$i.fastq.gz
ln -s $ddir/XRM-0812-3_*_$i.fq.gz $wdir/NT-14h_kdm4b_K9me3_rep2.R$i.fastq.gz
ln -s $ddir/XRM-0812-4_*_$i.fq.gz $wdir/NT-l2cell_kdm4b_input_rep1.R$i.fastq.gz
ln -s $ddir/XRM-0812-5_*_$i.fq.gz $wdir/NT-l2cell_kdm4b_K9me3_rep1.R$i.fastq.gz
ln -s $ddir/XRM-0812-6_*_$i.fq.gz $wdir/NT-l2cell_kdm4b_K9me3_rep2.R$i.fastq.gz
ln -s $ddir/XRM-0812-7_*_$i.fq.gz $wdir/NT-l2cell_kdm4b_K9me3_rep3.R$i.fastq.gz
done 


#####################
#2021年3月9日 kdm4d NT  L1C L2C NChIP K9me3
ddir=/home/share/DATA/clean_reads/GaoSR_mouse_nuclear_transfer_ChIP_20180301/20220222-4b
wdir=~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links/zb.20220222
mkdir -p $wdir ${wdir/0.links/1.qc} ${wdir/0.links/3.cutadapt} ${wdir/0.links/4.qc_cutadapt}

for i in 1 2; do 
ln -s $ddir/0209-XRM-1*_$i.fq.gz $wdir/NT-14h_kdm4b_K9me3_rep1.R$i.fastq.gz
ln -s $ddir/0209-XRM-2*_$i.fq.gz $wdir/NT-14h_kdm4b_K9me3_rep2.R$i.fastq.gz
ln -s $ddir/0209-XRM-3*_$i.fq.gz $wdir/NT-l2cell_kdm4b_K9me3_rep1.R$i.fastq.gz
ln -s $ddir/0209-XRM-4*_$i.fq.gz $wdir/NT-l2cell_kdm4b_K9me3_rep2.R$i.fastq.gz
ln -s $ddir/0209-XRM-5*_$i.fq.gz $wdir/NT-l2cell_kdm4b_K9me3_rep3.R$i.fastq.gz
done 

################################
#2022年3月9日 NT oeMcrs1 NChIP
ddir=/home/share/DATA/clean_reads/GaoSR_mouse_nuclear_transfer_ChIP_20180301/20220307-NT-Mcrs1-NChIP
wdir=~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links/zc.20220307
mkdir -p $wdir ${wdir/0.links/1.qc} ${wdir/0.links/3.cutadapt} ${wdir/0.links/4.qc_cutadapt}
for i in 1 2; do 
ln -s $ddir/XRM-0221-1*_$i.fq.gz $wdir/NT-ICM_oeMcrs1_input_rep1.R$i.fastq.gz
ln -s $ddir/XRM-0221-2*_$i.fq.gz $wdir/NT-ICM_oeMcrs1_input_rep2.R$i.fastq.gz
ln -s $ddir/XRM-0221-3*_$i.fq.gz $wdir/NT-ICM_oeMcrs1_K9me3_rep1.R$i.fastq.gz
ln -s $ddir/XRM-0221-4*_$i.fq.gz $wdir/NT-ICM_oeMcrs1_K9me3_rep2.R$i.fastq.gz
ln -s $ddir/XRM-0221-5*_$i.fq.gz $wdir/NT-TE_oeMcrs1_input_rep1.R$i.fastq.gz
ln -s $ddir/XRM-0221-6*_$i.fq.gz $wdir/NT-TE_oeMcrs1_input_rep2.R$i.fastq.gz
ln -s $ddir/XRM-0221-7*_$i.fq.gz $wdir/NT-TE_oeMcrs1_K9me3_rep1.R$i.fastq.gz
ln -s $ddir/XRM-0221-8*_$i.fq.gz $wdir/NT-TE_oeMcrs1_K9me3_rep2.R$i.fastq.gz
done 


########################2023年3月14日 revision data
obsutil config -i=JNXPZKWY32AVAPOYWEW7 -k=E4ssAc5davZzX8yypt26qBRdLlWw8oN2emyUlYYb -e=http://218.94.125.245:8999
obsutil ls obs://tongjiuniversity/20230314
obsutil cp obs://tongjiuniversity/20230313/XRM ~/data/9.NT-ChIP/revision.20230314 -f -r
obsutil cp obs://tongjiuniversity/20230314/XRM ~/data/9.NT-ChIP/revision.20230314 -f -r
obsutil cp obs://tongjiuniversity/20230314/MD ~/data/9.NT-ChIP/revision.20230314 -f -r

mkcd ~/data/9.NT-ChIP/revision.20230412

ddir=~/data/9.NT-ChIP/revision.20230314
wdir=~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links/zd.20230314
mkdir -p $wdir ${wdir/0.links/1.qc} ${wdir/0.links/3.cutadapt} ${wdir/0.links/4.qc_cutadapt}

for i in 1 2; do 
ln -fs $ddir/XRM-0307-1_*R$i*.fastq.gz $wdir/NT-6hpa_sih2_K9me3_rep2.R$i.fastq.gz
ln -fs $ddir/XRM-0307-2_*R$i*.fastq.gz $wdir/Sertoli-ICM_input_rep1.R$i.fastq.gz
ln -fs $ddir/XRM-0307-3_*R$i*.fastq.gz $wdir/Sertoli-ICM_input_rep2.R$i.fastq.gz
ln -fs $ddir/XRM-0307-4_*R$i*.fastq.gz $wdir/Sertoli-ICM_input_rep3.R$i.fastq.gz
ln -fs $ddir/XRM-0307-5_*R$i*.fastq.gz $wdir/Sertoli-TE_input_rep1.R$i.fastq.gz
ln -fs $ddir/XRM-0307-6_*R$i*.fastq.gz $wdir/Sertoli-TE_input_rep2.R$i.fastq.gz
ln -fs $ddir/XRM-0307-7_*R$i*.fastq.gz $wdir/Sertoli-TE_input_rep3.R$i.fastq.gz
ln -fs $ddir/XRM-0307-8_*R$i*.fastq.gz $wdir/Sertoli-ICM_K9me3_rep1.R$i.fastq.gz
ln -fs $ddir/XRM-0307-9_*R$i*.fastq.gz $wdir/Sertoli-ICM_K9me3_rep2.R$i.fastq.gz
ln -fs $ddir/XRM-0307-10_*R$i*.fastq.gz $wdir/Sertoli-ICM_K9me3_rep3.R$i.fastq.gz
ln -fs $ddir/XRM-0307-11_*R$i*.fastq.gz $wdir/Sertoli-TE_K9me3_rep1.R$i.fastq.gz
ln -fs $ddir/XRM-0307-12_*R$i*.fastq.gz $wdir/Sertoli-TE_K9me3_rep2.R$i.fastq.gz
ln -fs $ddir/XRM-0307-13_*R$i*.fastq.gz $wdir/Sertoli-TE_K9me3_rep3.R$i.fastq.gz
ln -fs $ddir/XRM-0307-14_*R$i*.fastq.gz $wdir/Sertoli-ICM_oeMcrs1_K9me3_rep1.R$i.fastq.gz
ln -fs $ddir/XRM-0307-15_*R$i*.fastq.gz $wdir/Sertoli-ICM_oeMcrs1_K9me3_rep2.R$i.fastq.gz
ln -fs $ddir/XRM-0307-16_*R$i*.fastq.gz $wdir/Sertoli-TE_oeMcrs1_K9me3_rep1.R$i.fastq.gz
ln -fs $ddir/XRM-0307-17_*R$i*.fastq.gz $wdir/Sertoli-TE_oeMcrs1_K9me3_rep2.R$i.fastq.gz
ln -fs $ddir/XRM-0307-18_*R$i*.fastq.gz $wdir/Sertoli-ICM_oeMcrs1_input_rep1.R$i.fastq.gz
ln -fs $ddir/XRM-0307-19_*R$i*.fastq.gz $wdir/Sertoli-ICM_oeMcrs1_input_rep2.R$i.fastq.gz
ln -fs $ddir/XRM-0307-20_*R$i*.fastq.gz $wdir/Sertoli-TE_oeMcrs1_input_rep1.R$i.fastq.gz
ln -fs $ddir/XRM-0307-21_*R$i*.fastq.gz $wdir/Sertoli-TE_oeMcrs1_input_rep2.R$i.fastq.gz
done 


nohup fastqc -o ${wdir/0.links/1.qc} $wdir/*gz -t 64 &
for i in $wdir/*R1.fastq.gz; do o=${i/0.links/3.cutadapt}
nohup cutadapt -j 8 --trim-n -q 25,25 -m 30 -a ATCAACGCAGAGTAC -a AGATCGGAAGAGC -A ATCAACGCAGAGTAC -A AGATCGGAAGAGC -o $o -p ${o/R1/R2} $i ${i/R1/R2} > ${o/R1.fastq.gz/log} &
done
wait
nohup fastqc -o ${wdir/0.links/4.qc_cutadapt} ${wdir/0.links/3.cutadapt}/*.gz -t 32 &
wait
cd ${wdir/0.links/1.qc}
multiqc --no-data-dir -f -dd 10 ./ -n $wdir/raw.qc.html
cd ${wdir/0.links/3.cutadapt}
multiqc --no-data-dir -f -dd 10 ./ -n $wdir/cutadapt.html
cd ${wdir/0.links/4.qc_cutadapt}
multiqc --no-data-dir -f -dd 10 ./ -n $wdir/clean.qc.html
wait-m

obsutil cp obs://tongjiuniversity/20230412/XRM ~/data/9.NT-ChIP/revision.20230412 -f -r
ddir=~/data/9.NT-ChIP/revision.20230412
wdir=~/workspace/9.NT-ChIP/0.cutadapt-qc/0.links/ze.20230412
mkdir -p $wdir ${wdir/0.links/1.qc} ${wdir/0.links/3.cutadapt} ${wdir/0.links/4.qc_cutadapt}
for i in 1 2; do 
ln -fs $ddir/XRM-0403-1_*R$i*.fastq.gz $wdir/6hpa_control_K9me3_rep3.R$i.fastq.gz
ln -fs $ddir/XRM-0403-2_*R$i*.fastq.gz $wdir/6hpa_control_K9me3_rep4.R$i.fastq.gz
ln -fs $ddir/XRM-0403-3_*R$i*.fastq.gz $wdir/6hpa_sih2_K9me3_rep3.R$i.fastq.gz
ln -fs $ddir/XRM-0403-4_*R$i*.fastq.gz $wdir/6hpa_sih2_K9me3_rep4.R$i.fastq.gz
done
