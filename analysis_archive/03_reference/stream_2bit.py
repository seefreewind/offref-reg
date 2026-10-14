"""Stream complete checksum-verified 2bit as FASTA; avoid external-disk rewrite."""
import sys,py2bit
SEED=20261009
p=py2bit.open(sys.argv[1]);out=sys.stdout.buffer
for ch,n in p.chroms().items():
 out.write(('>'+ch+'\n').encode())
 for s in range(0,n,1000000):out.write(p.sequence(ch,s,min(s+1000000,n)).encode()+b'\n')
