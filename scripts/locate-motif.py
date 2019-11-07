import sys
from Bio import SeqIO

print('python3 locate-motif.py FASTA OUTFILE')

FASTA = sys.argv[1]
OUT = sys.argv[2]
WIN = int(sys.argv[3])
motif = 'CCCT'
m = []

def findMotif(motif,strand):
	with open(FASTA,"r") as handle:
		for record in SeqIO.parse(handle,"fasta"):
			# 0-based numbering (2 is AGCCCT)	
			start=0
			while record.seq.find(motif,start) != -1:
				start = record.seq.find(motif,start)+1
				m.append([record.id,max(0,start-WIN),min(len(record.seq),start+WIN+len(motif)),strand])
findMotif('CCCT','+')
findMotif('AGGG','-')

with open(OUT,"w") as out_f:
	for match in m:
		out_f.write(match[0]+"\t"+str(match[1])+"\t"+str(match[2])+"\t"+match[3]+"\n")


		
		
