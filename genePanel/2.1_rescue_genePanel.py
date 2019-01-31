
from __future__ import division
import vcf
import pysam
import argparse
import numpy

def arguments():
    '''Parses the arguments from the program invocation'''
    
    #Call the argument parse
    parser = argparse.ArgumentParser() 
    
    #Specify arguments
    #parser.add_argument('--example', nargs='?', const=1, type=int, default=1)

    parser.add_argument("--vcf", type = str,
                        help="Path to SNV VCF file")
    parser.add_argument('--bamDir', default="", type = str,
                        help="Path to dir with BAMs"
                        )
    
    
    #parser.add_argument('--callers', type=int, default=1,
                        #help = 'Filter SNV variants called by this number of variant callers [1]')
    
    
    
    args = parser.parse_args()       
    return args



def checkBAMfreq(chrm, position, ref, nucl, bamfile, _type):
    bam = pysam.AlignmentFile(bamfile, 'rb') 
    nucList = []
    
    #Careful with pileup and positions, it follows python conventions. Truncate works to report only the column asked.
    if _type == 'snp':
        
        for pileupcol in bam.pileup(str(chrm), int(position)-1, int(position), truncate = True):
            for pileupread in pileupcol.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    nucList.append(pileupread.alignment.query_sequence[pileupread.query_position])
        count = nucList.count(nucl)
        cov = len(nucList)
    elif _type == 'del':
        coverages = []
        counts = []
        for pileupcol in bam.pileup(str(chrm), int(position), int(position)+len(ref), truncate = True):
            tmpCount = 0
            tmpCov = 0
            for pileupread in pileupcol.pileups:
                if pileupread.is_del:
                    tmpCount += 1
                tmpCov += 1
            coverages.append(tmpCov)
            counts.append(tmpCount)
        count = numpy.mean(counts)
        cov = numpy.mean(coverages)
    
    return count, cov

def vcfParser(vcfFile, bamdir):
    
    with open(vcfFile, 'r') as infile:
        for line in infile:
            line = line.rstrip()
            if line.startswith('#'):
                print line
                if line.startswith('#CHROM'):
                    tabline = line.split('\t')
                    samples = {}
                    for i, s in enumerate(tabline):                        
                        if i < 9:
                            continue
                        else:
                            samples[i] = s
                continue
            else:
                line = line.split('\t')
                for i, s in enumerate(line):
                    if i < 9:
                        continue
                    gt = s.split(':')[0]
                    if gt == "0/1" or gt == "1/1":
                        continue
                    elif gt == "./." or gt == "0/0":
                        sample = samples[i]
                        chrom = line[0]
                        pos = line[1]
                        ref = line[3]
                        alt = line[4]
                        
                        if len(ref) == len(alt) and len(ref) == 1 and len(alt) == 1:
                            t = 'snp'
                        elif len(ref) > len(alt):
                            t = 'del'
                        elif len (ref) < len(alt):
                            t = 'ins'
                            continue
                        altBam = bamdir + '/' + sample + '.bam'
                        
                        altCount, altDepth = checkBAMfreq(chrom, pos, ref, alt, altBam, t)
                        if altCount / altDepth >= 0.05:
                            line[i] = "0/1:%i,%i:%i" % (altDepth-altCount, altCount, altDepth)
                        else:
                            continue
                    else:
                        raise ValueError("Don't know what to do with GT "+gt)
                print '\t'.join(line)
                    


def main():
    args=arguments()
    vcfParser(args.vcf, args.bamDir)
    
   
###Execute main body

if __name__ == '__main__': main()
