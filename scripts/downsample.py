import sys
import random

n = int(sys.argv[0])
in_fname = sys.argv[1]
out_fname = sys.argv[2]

print ('n=', n)
print ('Reading', in_fname)
with open(in_fname) as in_f:
	lines = in_f.readlines()

print ('Writing out to', out_fname)
with open(out_fname, 'w') as out_f:
	out_f.write(random.sample(lines, n))


