import sys
import numpy as np

def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break

if __name__ == "__main__":
    # first argument is fasta file
    # second argument is the rng seed
    np.random.seed(int(sys.argv[2]))
    shape1, shape2 = 0.4, 0.6
    rulow, ruhigh = 0.2, 0.6
    includethr = 0.8
    # percentgapatends = 0.4
    locshape1, locshape2= 0.4, 0.4
    beginendtrim = 20
    res = []
    firsttime = True
    startindex = []
    y = []
    seqindex = 0

    for name, seq, qual in readfq(open(sys.argv[1])):
        seq = seq.strip()
        if firsttime:
            firsttime = False
            ru = np.random.uniform(rulow, ruhigh)
            y = np.floor(np.random.beta(ru, 1-ru, size=10000)*len(seq))
            startindex = np.floor(np.random.beta(locshape1, locshape2, size=10000)*(len(seq)-y))
            startindex[startindex < beginendtrim] = 0
            endindex = startindex + y
            endindex[endindex > (len(seq)-beginendtrim)] = len(seq)
            startindex = endindex - y
            # starts_ends = np.random.binomial(1, percentgapatends, 10000)
            # startindex = startindex * starts_ends
            # ends = starts_ends * np.random.binomial(1, 0.5, 10000)*(len(seq)-y)
            # startindex += ends

        if y[seqindex]/len(seq) <= includethr:
            seqnew = seq[0:int(startindex[seqindex])] + '-'*int(y[seqindex]) + seq[int(startindex[seqindex]+y[seqindex]):]
            if len(seqnew) != len(seq):
                print("zaa")
            assert len(seqnew) == len(seq)
            res.append(">" + name)
            res.append(seqnew)
        seqindex += 1
    print("\n".join(res))        
