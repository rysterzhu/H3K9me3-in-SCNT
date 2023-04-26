NT3=(CC 6h 14h NT-2cell NT-4cell NT-8cell NT-Morula NT-ICM NT-TE)
awk 'FNR==1{gsub("^.*/|.fraction$","",FILENAME)} {print $0,FILENAME}' $ddir/[C61]*fraction $ddir/NT*fraction | datamash -s groupby 1,2,3 collapse 4 > $wdir/fraction.NT.tab
awk -v samples="${NT3[*]}" 'BEGIN{len=split(samples,S);printf("chr\tstart\tend");for(i=1;i<=len;i++)printf("\t%s",S[i]);printf("\n")}
{m=($1"\t"$2"\t"$3);for(i=1;i<=len;i++){match($0,S[i]);m=m"\t"(RSTART==0?"unmarked":"marked")};print m}' $wdir/fraction.NT.tab > $wdir/fraction.NT.matrix
sed '1d' $wdir/fraction.NT.matrix | cut -f 4-7 | sort | uniq -c | awk -v samples="CC 6h 14h NT-2cell" 'BEGIN{$0=samples" count";$1=$1;print $0}
{for(i=2;i<=NF;i++)printf("%s\t",$i);print $1}' > $wdir/before-2cell.tab
sed '1d' $wdir/fraction.NT.matrix | cut -f 4-6 | sort | uniq -c | awk -v samples="CC 6h 14h" 'BEGIN{$0=samples" count";$1=$1;print $0}
{for(i=2;i<=NF;i++)printf("%s\t",$i);print $1}' > $wdir/before-14h.tab

