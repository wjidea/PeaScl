#! /usr/bin/env python

import argparse
from Bio import SeqIO

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('fasta_file')
	parser.add_argument('list_file')
	parser.add_argument('new_fasta')

	args = parser.parse_args()

	#Recover sequences
	ID_set =[line.strip().split('\t')[0] for line in open(args.list_file)]
	
	#Parse fasta 
	seqiter = SeqIO.parse(open(args.fasta_file), 'fasta')

	#Find match and write new fasta file
	SeqIO.write((seq for seq in seqiter if seq.id in ID_set), args.new_fasta, "fasta")
	

if __name__ == '__main__':
	main()
