import sys
import optparse
import re
from Bio import SeqIO


parser=optparse.OptionParser()
parser.add_option('-i', '--infile', help='', type='str')
parser.add_option('-f', '--fasta', help='', type='str')
parser.add_option('-o', '--outfile', help='', type='str')


options, args=parser.parse_args()


amb_d = {('A','C'):'M',
		 ('A','G'):'R',
		 ('A','T'):'W',
		 ('C','G'):'S',
		 ('C','T'):'Y',
		 ('G','T'):'K',
		 ('A','C','G'):'V',
		 ('A','C','T'):'H',
		 ('A','G','T'):'D',
		 ('C','G','T'):'B',
		 ('A','C','G','T'):'N'
		 }


def parse_seq(seq_f):
	seq = ''
	for record in SeqIO.parse(seq_f, "fasta"):
		name = record.description
		seq = str(record.seq)
	return(name, seq)


def parse_vcf(fname):
	d = {}
	with open(fname) as f:
		for line in f:
			arr = line.strip('\n').split('\t')
			if not line.startswith('#'):
				if arr[1] not in d:
					d[arr[1]] = []
				d[arr[1]].append(arr)
	return(d)


def analyze_freqs(d, amb_d):
	changes = {}
	for i in d:
		alt_freqs = 0
		alt_letters = []
		for j in d[i]:
			ref = j[3]
			alt = j[4]
			alt_freq = float(j[7].split(';')[1].replace('AF=', ''))
			if alt_freq >= 0.2:
				alt_freqs += alt_freq
				alt_letters.append(alt)
		if alt_freqs <= 0.8:
			alt_letters.append(ref)
		new_letter = amb_d[tuple(sorted(alt_letters))]
		changes[int(i)] = new_letter
	return(changes)


def apply_changes(seq, changes):
	seq_list = list(seq)
	for i in changes:
		seq_list[i-1] = changes[i]
	return(''.join(seq_list))


name, seq = parse_seq(options.fasta)
d = parse_vcf(options.infile)
changes = analyze_freqs(d, amb_d)
new_seq = apply_changes(seq, changes)

out = open(options.outfile, 'w')
out.write('>' + name + '\n' + new_seq + '\n')
out.close()
