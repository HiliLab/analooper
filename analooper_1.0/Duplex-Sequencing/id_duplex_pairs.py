import sys
import re
from argparse import ArgumentParser
from Bio import SeqIO
from Bio.Seq import Seq


def distance(a, b):
    return sum(map(lambda (x, y): 0 if x == y else 1, zip(a, b)))

def distance_less_than(k, a, b):
    return distance(a, b) < k

def revcomp(seq):
    
    #create a sequence object
    my_seq = Seq(seq)
    
    return str(my_seq.reverse_complement())

def main():
    parser = ArgumentParser()
    parser.add_argument('--infile', dest = 'infile', help = 'Input fasta file.  ', required=True)
    parser.add_argument('--outfile', dest = 'outfile', help = 'Output file.  ', required=True)
    parser.add_argument('--barcode_length', type = int, default = 12, dest = 'blength', help = 'Length of the duplex tag sequence. [12]')
    o=parser.parse_args()
    
    in1 = open(o.infile, 'rU')
    out = open(o.outfile, 'w')
    readDict = {}
    readTags = []
    for seqRecord in SeqIO.parse(in1,"fasta"):
        readDict[seqRecord.id.split('-')[0]] = {'seq' : seqRecord.seq, 'tseq' : seqRecord.id.split('-')[1]}
    in1.close()
    readTags = readDict.keys()

##    seqDict = SeqIO.index(o.infile, "fasta")

    for tag in readTags:
        switchTag = tag[o.blength:] + tag[:o.blength]
        tseq = readDict[tag]['tseq']
        if not distance_less_than(2,tseq,'GATTCG'):
            temp = tag
            tag = switchTag
            switchTag = temp
        
        if readDict.get(tag, None) and readDict.get(switchTag, None):
            out.write(">%s\n%s\n>%s\n%s\n" %(tag, readDict[tag]['seq'], switchTag, revcomp(str(readDict[switchTag]['seq']))))
            sys.stderr.write("Found duplex pair for: %s\t%s\n" % (tag, switchTag))
            readTags.remove(tag)
            readTags.remove(switchTag)
        else:
            sys.stderr.write("Missing duplex pair for: %s\t%s\n" % (tag, switchTag))
    out.close()

if __name__ == "__main__":
    main()
