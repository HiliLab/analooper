
'''
Tag To Header
Version 2.0.1
By Joe Hiatt, Scott Kennedy(1), Brendan Kohrn and Mike Schmitt(1)
(1) Department of Pathology, University of Washington School of Medicine, Seattle, WA 98195
March 24, 2014

Isolate duplex tags, move them from within the sequenced read to the header region, and remove the spacer region.  

usage: tag_to_header.py [-h] [--infile1 INFILE1] [--infile2 INFILE2]
                        [--outfile1 OUTFILE1] [--outfile2 OUTFILE2]
                        [--barcode_length BLENGTH] [--spacer_length SLENGTH]
                        [--read_out ROUT] [--adapter ADAPTERSEQ]

optional arguments:
  -h, --help            show this help message and exit
  --infile1 INFILE1     First input raw fastq file.
  --infile2 INFILE2     Second input raw fastq file.
  --outfile1 OUTFILE1   Output file for first fastq reads.
  --outfile2 OUTFILE2   Output file for second fastq reads.
  --barcode_length BLENGTH
                        Length of the duplex tag sequence. [12]
  --spacer_length SLENGTH
                        Length of the spacer sequences used. [5]
  --read_out ROUT       How often you want to be told what the program is
                        doing. [1000000]
  --adapter ADAPTERSEQ  Optional: Spacer sequence for filtering on the
                        presence of the spacer. This could be thrown off by
                        low quality scores.

'''

import sys
import re
from argparse import ArgumentParser
from Bio import SeqIO
from Bio.Seq import Seq

class fastQRead:
    def __init__(self, in1, in2, in3, in4):
        '''This class is meant to hold a single fastQ read.
        '''
        self.name=in1.strip().strip("@").replace(' ', '_')
        self.seq=in2.strip()
        self.spacer="+"
        self.qual=in4.strip()
        if len(self.seq)!=len(self.qual):
            raise ValueError("Sequence and quality scores of different lengths!/n%s/n%s/n%s/n%s" % (in1, in2, "+", in4))

    def __getitem__(self, index):
        '''This should allow slicing of the read to proceed properly.
        '''
        if isinstance(index, int):
            return self.seq[index]
        elif isinstance(index, slice):
            answer = self.__class__(self.name, self.seq[index], self.spacer, self.qual[index])
            return answer
        raise ValueError("Invalid index")


class fastQItterator:
    def __init__(self, inFile):
        '''This class will go through a fastQ file one line at a time.
        '''
        self.source=inFile
        self.eof=False
    
    def next(self):
        new=[]
        for j in xrange(4):
            try:
                tmp=self.source.next()
            except StopIteration:
                self.eof=True
                return("EOF") 
            new.append(tmp)
        newRead=fastQRead(new[0],new[1],new[2],new[3])
        return(newRead)

    def close(self):
        self.source.close()
        return(True)


class fastqWriter:
    def __init__(self, outFile):
        self.file=outFile
        self.firstLine=True
    
    def write(self, read):
        if self.firstLine==True:
            self.file.write("@" + read.name)
            self.firstLine=False
        else:
            self.file.write("\n@" + read.name)
        self.file.write("\n" + read.seq)
        self.file.write("\n" + read.spacer)
        self.file.write("\n" + read.qual)
        return(True)
    
    def close(self):
        self.file.close()
        return(True)

class fastaWriter:
    def __init__(self, outFile):
        self.file=outFile
        self.firstLine=True
    
    def write(self, name, seq):
        if self.firstLine==True:
            self.file.write(">" + name)
            self.firstLine=False
        else:
            self.file.write("\n>" + name)
        self.file.write("\n" + seq)
        return(True)
    
    def close(self):
        self.file.close()
        return(True)

def revcomp(seq):
    
    #create a sequence object
    my_seq = Seq(seq)
    
    return str(my_seq.reverse_complement())

def tagExtractFxn(x, blen):
    '''this is the function that extracts the UID tags from both the 
    forward and reverse read.  Assigns read1 the sequence from some 
    position to the end, then read2 from some position to the end, 
    then assigns tag1 from the 5'-end to length of the UID tag for 
    read1 and then read 2.
    '''
    ###return(x[0][:blen], x[1][:blen])
    return(tagMatch(x))

### Edit the value in {14} if you want to change size of barcode ###
def tagMatch(seq):
    pattern = '(?P<stag>[ATCG]{14})CAGTA(?P<templatestart>[ATCG]{18})(?P<motifs>[ATCG]{40})(?P<templateend>[ATCG]{18})TACTG(?P<rtag>[ATCG]{14})'
    match = re.search(pattern, seq)
    if match:
        return(match.group(0), match.groupdict())
    else:
        return None
    
    ##return(match.group('stag'), match.group('templatestart'), match.group('motifs'), match.group('templateend'), match.group('rtag')
    
def hdrRenameFxn(x, y, z):
    '''this function renames the header with the formatting of 
    *header coordinates,etc*, *index seq*, *tag from read1*, *tag from read2*, *spacer from this read*
    *read designation from original header*
    '''
    return("%s%s%s" % (x[:-1], y, z))

def main():
    parser = ArgumentParser()
    parser.add_argument('--infile1', dest = 'infile1', help = 'First input raw fastq file.  ', required=True)
    parser.add_argument('--outfile1', dest = 'outfile1', help = 'Output file for first fastq reads.  ', required=True)
    parser.add_argument('--barcode_length', type = int, default = 12, dest = 'blength', help = 'Length of the duplex tag sequence. [12]')
    parser.add_argument('--spacer_length', type = int, default = 5, dest = 'slength', help = 'Length of the spacer sequences used. [5]')
    parser.add_argument('--read_out', type = int, default = 1000000, dest = 'rOut', help = 'How often you want to be told what the program is doing. [1000000]')
    parser.add_argument('--adapter',  default = None,  dest = 'adapterSeq', help = 'Optional: Spacer sequence for filtering on the presence of the spacer.  This could be thrown off by low quality scores.')
    o=parser.parse_args()


    in1=fastQItterator(open(o.infile1, 'rU'))
    out1=fastaWriter(open(o.outfile1, 'w'))

    ctr=0
    nospacer = 0
    goodreads = 0
    badtag = 0
    oldBad = 0
    isEOF=False

    while isEOF==False:
        read1 = in1.next()

        if read1 == "EOF":
            isEOF = True
        else:
            
            ctr += 1
            if o.adapterSeq != None and (read1.seq[o.blength:o.blength + o.slength] != o.adapterSeq):
                nospacer += 1
            else:
                #extract tags
                r1parts = tagExtractFxn(read1.seq,o.blength)
                if not r1parts:
                    continue

                ftag = r1parts[1]['stag']
                rtag = revcomp(r1parts[1]['rtag'])
                
                #header reconstruction
                #read1.name = hdrRenameFxn(read1.name, ftag, rtag)
                readName = read1.name
                read1.name = ftag + rtag + "-" + r1parts[1]['templatestart'][:6]

                tag1 = ftag

                #fastq reconstruction
                if (tag1.isalpha() and tag1.count('N') == 0) :
                    rOut1 = r1parts[1]['motifs']
                    out1.write(read1.name, rOut1)
                    encodedStr = "\t".join((readName, read1.name, r1parts[1]['motifs']))
                    print encodedStr
                    goodreads += 1
                else:
                    badtag += 1
            if ctr%o.rOut==0:
                sys.stderr.write("Total sequences processed: %s\n" % (ctr))
                sys.stderr.write("Sequences passing filter: %s\n" % (goodreads))
                sys.stderr.write("Missing spacers: %s\n" % (nospacer))
                sys.stderr.write("Bad tags: %s\n\n" % (badtag))
                if badtag == oldBad+o.rOut:
                    sys.stderr.write("Warning!  Potential file error between lines %s and %s.  " % ((ctr-o.rOut)*4,(ctr)*4))
                oldBad = badtag

    in1.close()
    out1.close()


    sys.stderr.write("Summary statistics:\n")
    sys.stderr.write("Total sequences processed: %s\n" % (ctr))
    sys.stderr.write("Good sequences: %s\n" % (goodreads))
    sys.stderr.write("Missing spacers: %s\n" % (nospacer))
    sys.stderr.write("Bad tags: %s\n\n" % (badtag))

if __name__ == "__main__":
    main()
