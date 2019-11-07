cat pval-wt.txt | sort | awk '{idx=$1}{a[idx]=(idx in a)? a[idx]","$2:$2}END{for(i in a) print i,a[i]}' | sort > collapsed.txt

cat pval-wt.txt | sort > tmp
join tmp collapsed.txt > pval-all.txt
