#! /usr/bin/python

from __future__ import division
import vcf
import pysam
import argparse
import matplotlib.pyplot as plt

def arguments():
    '''Parses the arguments from the program invocation'''
    
    #Call the argument parse
    parser = argparse.ArgumentParser() 
    
    #Specify arguments
    #parser.add_argument('--example', nargs='?', const=1, type=int, default=1)

    parser.add_argument('-f', "--readfile", type = str,
                        help="Path to VCF file", required=True)
    parser.add_argument('-b', '--bamdir', 
                        help="Directory with bam files",required=True )
    parser.add_argument('-a', '--afthreshold', default=0.05, type = float,
                        help="Allele frequency threshold (per 1) [0.05]")
    parser.add_argument('-r', '--readthreshold', default='0', type = float,
                        help="# of reads supporting the alt allele for rescue [0]")
    parser.add_argument('--report', action="store_true",
                        help="Report the rescued variants ")
    parser.add_argument('--nonsyn', action="store_true",
                        help="Use only non-synonimous variants ")
    parser.add_argument('--noref', action="store_true",
                        help="No reference available ")
    args = parser.parse_args()       
    return args



def blood_dict():
    '''Establish the relation of each sample with the corresponding blood sample to use in the BAM rescue'''
    bdict = {'T1':'B1', 'O1':'B1',
             'T2':'B1', 'O2':'B1',
             'T4':'B-IN', 'O4':'B-IN',
             'T5':'B-IN', 'O5':'B-IN',
             'T7':'B-IN', 'O7':'B-IN',
             'T8':'B-IN', 'O8':'B-IN',
             'T9':'B9', 'O9':'B9', 'O9-S':'B9',
             'T10':'B9', 'O10':'B9', 'O10-S':'B9',
             'T11':'B11', 'O11':'B11',
             'T12':'B12', 'O12':'B12',
             'T13':'B13', 'O13':'B13', 'O13-S':'B13',
             'T14':'B14', 'O14':'B14',
             'T18':'B14', 'O18':'B14',
             'T15':'B15', 'O15':'B15', 'O15-S':'B15',
             'T16':'B16', 'O16':'B16',
             'T17':'B16', 'O17':'B16',
             'T21':'B16', 'O21':'B16',
             'T25':'B25', 'O25':'B25', 'O25-S':'B25',
             'T37':'B25', 'O37':'B25',
             'T47':'B25', 
             'T61':'B25', 'O61':'B25', 'O61-S':'B25',
             'T74':'B25', 'O74':'B25',
             'T29':'B29', 'O29':'B29',
             'T40':'B40', 'O40':'B40',
             'T41':'B40', 'O41':'B40',
             'T42':'B42', 'O42':'B42',
             'T27':'B16', 'O27':'B27',
             'T49':'B49', 'O49':'B49',
             'T50':'B49', 'O50':'B49',
             'T52':'B52', 'O52':'B52',
             'T53':'B52', 'O53':'B52',
             'T54':'B54', 'O54':'B54',
             'T55':'B55', 'O55':'B55',
             'T56':'B55', 'O56':'B55',
             'T58':'B58', 'O58':'B58',
             'T59':'B58', 'O59':'B58',
             'T65':'B65', 'O65':'B65',
             'T66':'B65', 'O66':'B65', 'T66-B':'B65',
             'T68':'B67', 'O68':'B67',
             'T70':'B67', 'O70':'B67'
             
             
             }
    return bdict


def vcfReader(vcfFilename, nonsynonimous):
    
    with open(vcfFilename, 'r') as vcfFile:
        vcfReader = vcf.Reader(vcfFile)
        samples = vcfReader.samples        
        
        
        annotationRank = ["chromosome_number_variation", 'exon_loss_variant', 'frameshift_variant', 'stop_gained', 'stop_lost', 'start_lost', 
             'splice_acceptor_variant', 'splice_donor_variant', 'rare_amino_acid_variant', 'missense_variant', 'inframe_insertion',
             'disruptive_inframe_insertion', 'inframe_deletion', 'disruptive_inframe_deletion', '5_prime_UTR_truncation+exon_loss_variant',
             '3_prime_UTR_truncation+exon_loss', 'splice_branch_variant', 'splice_region_variant', 'splice_branch_variant', 
             'stop_retained_variant', 'initiator_codon_variant', '5_prime_UTR_premature_start_codon_gain_variant', ]
        vcfDict = {}
        #vcfDict = {Gene:[{chrom:1, pos:183483843, id:rs822883, ref:T, alt:TTC, qual:678, filter:string, info:string, format:string}, {dict2}]}
        #vcfDict Keys are genes, values are a list with dictionaries for each variants
        
        
        for record in vcfReader:
            
            if nonsynonimous:
                try:
                    ann = []
                    flag = False
                    for a in record.INFO['ANN']:
                        tmp = a.strip().split('|')[1]
                        for a in annotationRank:
                            #SOMETIMES THE ANNOTATION ARE COMPOSED
                            if a in tmp:
                                flag = True
                        break
                    if not flag:
                        continue
                    
                except KeyError:
                    continue
                           
            tmpDict = {'pos':record.POS, 'id':record.ID, 'ref':record.REF, 
                        'alt':record.ALT, 'qual':record.QUAL, 'filter':record.FILTER, 'info':record.INFO, 'format':record.FORMAT, 'samples':record.samples}
            try:
                vcfDict[record.CHROM].append(tmpDict)
            except KeyError:
                vcfDict[record.CHROM] = [tmpDict]
            
    
    
    return vcfDict, samples

def bamCheck(bamfile, chrm, position, allele):
    
    bam = pysam.AlignmentFile(bamfile, 'rb') 
    nucList = []
    
    #Careful with pileup and positions, it follows python conventions. Truncate works to report only the column asked.
    for pileupcol in bam.pileup(str(chrm), int(position)-1, int(position), truncate = True):
        for pileupread in pileupcol.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                nucList.append(pileupread.alignment.query_sequence[pileupread.query_position])
    if len(nucList) == 0:
        freq = 0
    else:
        freq = nucList.count(allele)/len(nucList)
    return freq
                    




def bamCheckMaster(varDict, bamDir, afThreshold, refFlag):
    
    varComp = {}
    count = 1
    #test = []
    bloodDict = blood_dict()
    for chrm, recordList in varDict.items():
        #if chrm != '21':
            #continue        
        for record in recordList:
            varComp[count] = {'chrm':chrm, 'pos':record['pos'], 'ref':record['ref'], 'alt':record['alt']}
            checkList = []
            for call in record['samples']:
                #0 for not there
                #1 for there
                #2 for rescued
                
                sample = call.sample
                gt = call['GT']
                
                if gt != './.' and gt != '0/0':
                    varComp[count][sample] = 1
                else:
                    varComp[count][sample] = 0
                    #if len(record['alt']) == 1 and len(record['alt'][0]) == 1 and len(record['ref']) == 1:
                    bamfile = bamDir + '/' + sample + '.bam'
                    altfreq = bamCheck(bamfile, chrm, record['pos'], record['alt'][0])
                    
                    if not refFlag:
                        reference = bamDir + '/' + bloodDict[sample] + '.bam'
                        bloodfreq = bamCheck(reference, chrm, record['pos'], record['alt'][0])
                    else: 
                        bloodfreq = 0
                    if altfreq >= afThreshold and (bloodfreq == 0 or refFlag):
                        varComp[count][sample] = 2
                        #test.append(altfreq*100)
                varComp[count]['site'] = call.site
            count += 1
     
    #plt.hist(test) 
    return varComp




def printMaster(bamCheckedDict, samples):
    printlist = []
    header = ['variant'] + samples
    printlist.append(header)
    for count, var in bamCheckedDict.items():
        tmp = [str(count)]
        for sample in samples:
            tmp.append(str(var[sample]))
        printlist.append(tmp)
    for line in printlist:
        print '\t'.join(line)
    return
            
            
def printMaster_report(bamCheckedDict, samples):
    for count, var in bamCheckedDict.items():
        rescueFlag = False
        for sample in samples:
            if var[sample] == 2:
                rescueFlag = True
        if rescueFlag:
            print var['site']
    return 




###MAIN BODY###

def main():
    
    args=arguments()
    varDict, sampleList = vcfReader(args.readfile, args.nonsyn)
    bamCheckedDict = bamCheckMaster(varDict, args.bamdir, args.afthreshold, args.noref)    

    
    if args.report:
        printMaster_report(bamCheckedDict, sampleList)
    else:
        printMaster(bamCheckedDict, sampleList)
    
    
###Execute main body

if __name__ == '__main__': main()
    
    
      