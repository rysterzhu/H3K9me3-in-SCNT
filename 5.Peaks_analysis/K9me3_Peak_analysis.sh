
mkd P3;for i in *broadPeak; do awk '$8>3' $i > P3/$i;done;cd P3
for i in *_peaks.broadPeak; do k=$(basename $i _peaks.broadPeak)
awk '{a=int($2/1000+0.5);b=int($3/1000+0.5)
for(i=a;i<b;i++){print $1,i*1000,i*1000+1000}
}' $i | sort -k1,1 -k2n,2 | uniq > $k.fraction &
done
for i in *_peaks.broadPeak; do k=$(basename $i _peaks.broadPeak)
awk '{a=int($2/500+0.5);b=int($3/500+0.5)
for(i=a;i<b;i++){print $1,i*500,i*500+500}
}' $i | sort -k1,1 -k2n,2 | uniq > 2.500bp-frations/$k.fraction &
done

echo -e "sample\ttotalLength\tpeakNumber\tcoverage\tmeanLength" > peaks.stat;for i in ${NTF[@]}; do awk -v i=$i '{a+=1;b+=$3-$2} END{print i,b,a,b/27e8,b/a}' ${i}_peaks.broadPeak >> peaks.stat;done

(echo -e "sample\ttotalPrevious\ttotalNext\tboth\testablish\tdisappear" > est-dis.NT.stat
for n in CC:6h 6h:14h 14h:NT-2cell NT-2cell:NT-4cell NT-4cell:NT-8cell NT-8cell:NT-Morula NT-Morula:NT-ICM NT-Morula:NT-TE; do P=${n%%:*}.fraction;N=${n##*:}.fraction
awk -v name=${n##*:} 'NR==FNR{a[$1"-"$2];totala+=1} NR>FNR{b[$1"-"$2];totalb+=1}
END{for(i in a){if(i in b){both+=1}};print name,totala,totalb,both,(totalb-both)/totalb,(both-totala)/totala}' $P $N >> est-dis.NT.stat
done)

(echo -e "sample\ttotalPrevious\ttotalNext\tboth\testablish\tdisappear"
for n in Oocyte:PN3 Sperm:PN3 PN3:PN5 PN5:NF-2cell NF-2cell:NF-4cell NF-4cell:NF-8cell NF-8cell:NF-Morula NF-Morula:NF-ICM NF-Morula:NF-TE; do P=${n%%:*}.fraction;N=${n##*:}.fraction
awk -v name=${n##*:} 'NR==FNR{a[$1"-"$2];totala+=1} NR>FNR{b[$1"-"$2];totalb+=1}
END{for(i in a){if(i in b){both+=1}};print name,totala,totalb,both,(totalb-both)/totalb,(both-totala)/totala}' $P $N
done)


#直接以bp单位计算
for n in cc:6h 6h:14h 14h:2cell 2cell:4cell 4cell:8cell 8cell:morula morula:icm morula:te morula:blast; do P=${n%%:*}_peaks.broadPeak;N=${n##*:}_peaks.broadPeak
totalN=`awk '{t+=$3-$2} END{print t}' $N`
totalP=`awk '{t+=$3-$2} END{print t}' $P`
establish=`bedtools intersect -a $N -b $P -v | awk '{t+=$3-$2} END{print t}'`
disappear=`bedtools intersect -a $P -b $N -v | awk '{t+=$3-$2} END{print t}'`
echo -e "${n##*:}\t$totalN\t$establish\t$totalP\t$disappear"
done

##################Number of LTR marked by K9me3 domains
#LTR 数量972593 平均长度328.848
LTR=~/ann/mm10.LTR.bed;type=LTR
total=`cat $LTR | wc -l `
echo -e "sample\tnumber\tratio\ttype" > marked.stat.tab
for j in ${NT2[@]}; do i=${j##*:}
g=`intersectBed -a $LTR -b ${i}_peaks.broadPeak -r -e -f 0.5 -u -wa | wc -l`
awk -v i=$i -v g=$g -v total=$total -v t=$type 'BEGIN{printf("%s\t%d\t%g\t%s\n",i,g,g/total,t)}' >> marked.stat.tab
done

##################Number of gene promoter(-1kb~1kb) marked by K9me3 domains
LTR=~/ann/mm10_promoter_1000-1000.bed;type=promoter-1kb
total=`cat $LTR | wc -l `
for j in ${NT2[@]}; do i=${j##*:}
g=`intersectBed -a $LTR -b ${i}_peaks.broadPeak -r -e -f 0.5 -u -wa | wc -l`
awk -v i=$i -v g=$g -v total=$total -v t=$type 'BEGIN{printf("%s\t%d\t%g\t%s\n",i,g,g/total,t)}' >> marked.stat.tab
done

LTR=~/ann/mm10.CGI.bed;type=CpG-Island
total=`cat $LTR | wc -l `
for j in ${NT2[@]}; do i=${j##*:}
g=`intersectBed -a $LTR -b ${i}_peaks.broadPeak -r -e -f 0.5 -u -wa | wc -l`
awk -v i=$i -v g=$g -v total=$total -v t=$type 'BEGIN{printf("%s\t%d\t%g\t%s\n",i,g,g/total,t)}' >> marked.stat.tab
done

##############NT NF specific fractions
echo -e "stage\tboth\tNT\tNF" >specific.tab
for n in "CC:Oocyte" "CC:Sperm" "6h:PN3" "14h:PN5";do P=${n%%:*};N=${n##*:}
awk -v stage=$N 'NR==FNR{a[$1"\t"$2];t1+=1} NR>FNR{b[$1"\t"$2];t2+=1} END{for(i in a){if(i in b)both+=1};print stage,both,t1-both,t2-both}' $P.fraction $N.fraction  >> specific.tab
done
for n in "2cell" "4cell" "8cell" "Morula" "ICM" "TE";do
awk -v stage=$n 'NR==FNR{a[$1"\t"$2];t1+=1} NR>FNR{b[$1"\t"$2];t2+=1} END{for(i in a){if(i in b)both+=1};print stage,both,t1-both,t2-both}' NT-$n.fraction NF-$n.fraction  >> specific.tab
done

###NT specific peaks bp
echo -e "stage\tBoth\tNT-specific" > specific.NT-peaks.tab
for n in "CC:Oocyte" "6h:PN3" "14h:PN5";do P=${n%%:*};N=${n##*:}
intersectBed -a ${P}_peaks.broadPeak -b ${N}_peaks.broadPeak -wao | awk -v stage=$P 'NR==FNR{a+=$NF} NR>FNR{b+=$3-$2} END{print stage,a,b-a}' - ${P}_peaks.broadPeak >> specific.NT-peaks.tab
done
for n in "2cell" "4cell" "8cell" "Morula" "ICM" "TE";do
intersectBed -a NT-${n}_peaks.broadPeak -b NF-${n}_peaks.broadPeak -wao | awk -v stage=$n 'NR==FNR{a+=$NF} NR>FNR{b+=$3-$2} END{print stage,a,b-a}' - NT-${n}_peaks.broadPeak >> specific.NT-peaks.tab
done

##NT inherit CC peaks bp
echo -e "stage\toverlapCC\tspecific" > overlapCC.NT-peaks.tab
NTC=(CC 6h 14h NT-2cell NT-4cell NT-8cell NT-Morula NT-ICM NT-TE)
for i in ${NTC[@]}; do
intersectBed -a ${i}_peaks.broadPeak -b CC_peaks.broadPeak -wao | awk -v stage=$i 'NR==FNR{a+=$NF} NR>FNR{b+=$3-$2} END{print stage,a,b-a}' - ${i}_peaks.broadPeak >> overlapCC.NT-peaks.tab
done

#NT specific with NF peaks, inherit from CC
echo -e "stage\toverlapCC\tspecific" > overlapCC.NT-specific-peaks.tab
for n in "CC:Oocyte" "6h:PN3" "14h:PN5" "NT-2cell:NF-2cell" "NT-4cell:NF-4cell" "NT-8cell:NF-8cell" "NT-2cell:NF-2cell" "NT-Morula:NF-Morula" "NT-ICM:NF-ICM" "NT-TE:NF-TE";do P=${n%%:*};N=${n##*:}
subtractBed -a ${P}_peaks.broadPeak -b ${N}_peaks.broadPeak > temp
intersectBed -a temp -b CC_peaks.broadPeak -wao | \
awk -v stage=$P 'NR==FNR{a+=$NF} NR>FNR{b+=$3-$2} END{print stage,a,b-a}' - temp >> overlapCC.NT-specific-peaks.tab
done

###NF and NT peaks , overlap with CC
echo -e "stage\toverlapCC\tspecific" > overlapCC.NF-NT-peaks.tab
NTC=(6h 14h NT-2cell NT-4cell NT-8cell NT-Morula NT-ICM NT-TE)
NFC=(PN3 PN5 NF-2cell NF-4cell NF-8cell NF-Morula NF-ICM NF-TE)
for i in ${NTC[@]} ${NFC[@]}; do
intersectBed -a ${i}_peaks.broadPeak -b CC_peaks.broadPeak -wao | awk -v stage=$i 'NR==FNR{a+=$NF} NR>FNR{b+=$3-$2} END{print stage,a,b-a}' - ${i}_peaks.broadPeak >> overlapCC.NF-NT-peaks.tab
done
