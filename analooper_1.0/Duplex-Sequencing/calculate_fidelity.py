import sys
import re
from argparse import ArgumentParser
from Bio import SeqIO
from Bio.Seq import Seq
from itertools import groupby

def fasta_iter(fasta_name):
    """
    given a fasta file. yield tuples of header, sequence
    """
    fh = open(fasta_name, 'rU')
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        header = header.next()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.next())
        yield header, seq


def distance(a, b):
    return sum(map(lambda (x, y): 0 if x == y else 1, zip(a, b)))

def distance_less_than(k, a, b):
    return distance(a, b) < k

def revcomp(seq):
    
    #create a sequence object
    my_seq = Seq(seq)
    
    return str(my_seq.reverse_complement())

def compare_motifs(motif_size, seq1, seq2):
    regex_pattern = '[AGTC]{%s}' % str(motif_size)
    sList1 = re.findall(regex_pattern, seq1)
    sList2 = re.findall(regex_pattern, seq2)
    
    distances = [distance(s1,s2) for s1, s2 in zip(sList1, sList2)]
    #badMotifs = [(s1,s2,1) for s1, s2, d in zip(sList1, sList2, distances) if d > 0]
    #goodMotifs = [(s1,s2,0) for s1, s2, d in zip(sList1, sList2, distances) if d == 0]
    return zip(sList1, sList2, distances)


def main():
    parser = ArgumentParser()
    parser.add_argument('--infile', dest = 'infile', help = 'Input fasta file.  ', required=True)
    parser.add_argument('--outfile', dest = 'outfile', help = 'Output file.  ', required=True)
    parser.add_argument('--motif_length', type = int, default = 5, dest = 'mlength', help = 'Length of the template motif sequence. [5]')
    o=parser.parse_args()
    
    in1 = fasta_iter(o.infile)
    out = open(o.outfile, 'w')
    motifDict = {}

    while True:
        try:
            (h1, seq1) = in1.next()
            (h2, seq2) = in1.next()

            motifList = compare_motifs(o.mlength,seq1,seq2)
            for ts, ds, d in motifList:
                rcTs = revcomp(ts)
                if not rcTs in motifDict:
                    motifDict[rcTs] = {}
                    motifDict[rcTs]['bad'] = []
                    motifDict[rcTs]['count'] = 0
                    motifDict[rcTs]['badCount'] = 0
                    motifDict[rcTs]['ts'] = ts
                if d > 0:
                    motifDict[rcTs]['bad'].append(ds)
                    motifDict[rcTs]['badCount'] = motifDict[rcTs]['badCount'] + 1
                motifDict[rcTs]['count'] = motifDict[rcTs]['count'] + 1
                    
        except StopIteration:
            break

    for motif in motifDict.iterkeys():
        out.write("%s\t%s\t%s\t%s\t%s\n" % (motif, motifDict[motif]['ts'], motifDict[motif]['count'], motifDict[motif]['badCount'], "\t".join(motifDict[motif]['bad'])))
        
    out.close()


if __name__ == "__main__":
    main()
