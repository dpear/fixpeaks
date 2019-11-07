cd /data/daniela/Projects/fixpeaks/h99genome

cat revised_h99_genome.fasta |  tr -d "\n" | tr ">" "\n" | tr -d 'chr[0-9]' | grep -i CCCT -obna > CCCT-matches.txt
cat revised_h99_genome.fasta |  tr -d "\n" | tr ">" "\n" | tr -d 'chr[0-9]' | grep -i AGGG -obna > AGGG-matches.txt
mkdir chrfiles
cat CCCT-matches.txt | sed -r 's/:/\t/g' | awk '{print "chr" $1-1 "\t" $2 "\t" $3}' | sed 's/chr15/chrM/' | sed -r 's/\t/|/g' | awk -F\| '{print>"chrfiles/"$1}'

cat revised_h99_genome.fasta | grep -P '\>' | sed 's/|//g' | awk '{print $1 "\t" $4}' | sed 's/>//' | sed 's/length=//' > chrlengths.txt

echo 0 > chrcumulative.txt
cat chrlengths.txt | awk '{print $2}' | awk '{total += $0; $0 = total}1' >> chrcumulative.txt
paste chrlengths.txt chrcumulative.txt | awk '{print $1 "|" $3}' > chrsubtract.txt
cat chrsubtract.txt |  awk -F\| '{print>>"chrfiles/"$1}'
