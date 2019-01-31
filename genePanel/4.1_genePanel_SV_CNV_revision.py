#! /usr/bin/python

from __future__ import division
import vcf
import argparse
from patientDictionary import patientDict
import geneDicts

def arguments():
    '''Parses the arguments from the program invocation'''
    
    #Call the argument parse
    parser = argparse.ArgumentParser() 
    
    #Specify arguments
    #parser.add_argument('--example', nargs='?', const=1, type=int, default=1)

    parser.add_argument("--snv", type = str,
                        help="Path to SNV VCF file")
    parser.add_argument("--cnv", type = str,
                        help="Path to FREEC CNV file")
    parser.add_argument("--sv", type = str,
                        help="Path to SV VCF file")      
    parser.add_argument('-p', '--panel', default="december2018",
                        help="Gene panel to use",                        
                        )
    
    parser.add_argument('-o', '--outfile', default="",
                        help="Output file instead of STDOUT [STDOUT]",
                        )
    parser.add_argument('--patient', default="",
                        help="Patient label []",
                        )
    parser.add_argument('--sample', default="",
                        help="Sample type",
                        )
    
    #parser.add_argument('--callers', type=int, default=1,
                        #help = 'Filter SNV variants called by this number of variant callers [1]')
    
    
    
    args = parser.parse_args()       
    return args


def genepanel():
    '''Gene Panel to use, in the form of a dictionary. It is a bit redundant since I filtered the VCFs already with SnpSift,
    but it is the result of legacy'''
    {
    'CXCR1':("2",219027568,219031725),
    'CDK12':("17",37617739,37721160),
    'CCNE1':("19",30302805,30315215),
    'E2F1':("20",32263489,32274210),
    'MYC':("8",128747680,128753674),
    'RB1':("13",48877883,49056122),
    'ARID1A':("1",27022524,27108595),
    'BRAF':("7",140424943,140624564),
    'KRAS':("12",25357723,25403870),
    'NRAS':("1",115247090,115259515),
    'PPP2R1A':("19",52693274,52730687),
    'PIK3CA':("3",178865902,178957881),
    'PTEN':("10",89622870,89731687),
    'TP53':("17",7565097,7590856),
    'MDC1':("6",30667584,30685666),
    'PARP1':("1",224102741,226595780),
    'FGFR2':("10",123237848,123357972),
    'SLC19A1':("21",46913486,46964325),
    'APOB':("2",21224301,21266945),
    'BCORL1':("X",129115083,129192058),
    'BRCA1':("17",41196312,41322290),
    'BRCA2':("13",32889611,32973805),
    'CREBBP':("16",3775055,3930727),
    'FLG2':("1",152321213,152332482),
    'MYH1':("17",10395629,10421859),
    'MYH2':("17",10424465,10453274),
    'PPP1R3A':("7",113516832,113715975),
    'SYNE1':("6",152442819,152958936),
    'TRIOBP':("22",38093011,38172563),
    'USH2A':("1",215796236,216596738),
    'MMP16':("8",89044237,89340254),
    'MMP26':("11",4726157,5013659),
    'DAB2':("5",39371780,39462402),
    'TP53BP1':("15",43699407,43802926),
    'CDKN2A':("9",21967751,21995300),
    'CDKN2B':("9",22002902,22009280),
    'RPTOR':("17",78518619,78940171),
    'AURKA':("20",54944445,54967393),
    'FGFR1':("8",38268656,38326352),
    'ATR':("3",142168077,142297668),
    'CIC':("19",42772689,42799948),
    'FAT3':("11",92085262,92629618),
    'FAT4':("4",126237554,126414087),
    'CSMD3':("8",113235157,114449328),
    'NDRG1':('8',134249414,134314265),
    'PLEC':("8",144989321,145050902),    
    'RECQL4':('8',145736667,145743229),
    'MECOM':("3",168801287,169381406),
    'PTK2':("8",141667999,142012315),
    'AGO2':('8',141541264,141645718),
    'PRKCI':("3",169940153,170023769),
    'EXT1':('8',118806729,119124092),
    'NF1':("17",29421945,29708905),
    'GTSE1':('22',46692638,46726707),
    'WWOX':("16",78133310,79246564),
    'LRP1B':("2",140988992,142889270),
    'GABRA6':("5",160974069,161129599),
    'PARK2':('6',161768452,163148803),
    'RASA1':('5',86563705,86687748),     
    'MAP3K1':("5",56111401,56191979),
    'MAP3K4':("6",161412759,161551917)
   }




    
    
def snvReader(vcfFilename, geneDict):

    with open(vcfFilename, 'r') as vcfFile:
        vcfReader = vcf.Reader(vcfFile)
       
                
        
        vcfDict = {}
        #vcfDict = {Gene:[{chrom:1, pos:183483843, id:rs822883, ref:T, alt:TTC, qual:678, filter:string, info:string, format:string}, {dict2}]}
        #vcfDict Keys are genes, values are a list with dictionaries for each variants
        s = vcfReader.samples[0]
        for gene, coords in geneDict.items():
            chrm, start, end = coords
            vcfDict[gene] = []
            try:
                for record in vcfReader.fetch(chrm, start-1, end):
                    if record.FILTER != []:
                        continue
                                       
                    genotypeRaw = record.genotype(s)['GT']
                    if genotypeRaw == '0/1':
                        genotype = 'HET'
                    elif genotypeRaw == '1/1':
                        genotype = 'HOM'
                    else:
                        continue
                    tmpDict = {'chrom':record.CHROM, 'pos':record.POS, 'id':record.ID, 'ref':record.REF, 
                            'alt':record.ALT, 'qual':record.QUAL, 'filter':record.FILTER, 'info':record.INFO, 'format':record.FORMAT, 'genotype':genotype}
                    vcfDict[gene].append(tmpDict)
            except ValueError: 
                continue
            
                                                                                                              
                #print '\t'.join([str(record.CHROM), str(record.POS), record.ID, 
                                #record.REF, record.ALT, str(record.QUAL), 
                                #record.FILTER, record.INFO, record.FORMAT])
    return vcfDict



def cnvReader(cnvFilename, geneDict):
    cnvWhole = {}
    #Chr:[(start, end, cn, type), (start, end, cn, type)]
    with open(cnvFilename, 'r') as infile:
        for line in infile:
            chrm, start, end, cn, cntype = line.strip().split('\t')
            if chrm == 'X':
                chrm = 0
            elif chrm == 'Y':
                chrm = 42
            chrm = int(chrm)
            start = int(start)
            end = int(end)
            cn = int(cn)
            if not chrm in cnvWhole:
                cnvWhole[chrm] = [(start, end, cn, cntype)]
            else:
                cnvWhole[chrm].append((start, end, cn, cntype))
    geneCNV = {}
    for gene, coords in geneDict.items():
        geneChr, geneStart, geneEnd = (coords[0], coords[1], coords[2])
        if geneChr == 'X':
            geneChr = 0
        try:
            for cnv in cnvWhole[int(geneChr)]:
                cnStart, cnEnd, cnCN, cnType = cnv
                if geneEnd > cnStart and geneStart < cnEnd:
                    geneCNV[gene] = cnv
                    break
        except KeyError:
            geneCNV[gene] = (0,0,2,'normal')
        if not gene in geneCNV:
            geneCNV[gene] = (0,0,2,'normal')
    return geneCNV



def svReader(svFilename, geneDict):
    '''This info is not used anymore, still it is here... I know...'''
    ##NOT THAT MANY SVs, so parse manually
    
    svWhole = {}
    with open(svFilename, 'r') as infile:
        for line in infile:
            if line.startswith('#'):
                continue
            try:
                chrm, pos, ident, ref, alt, qual, filt, info, form, refGen, sampleGen = line.strip().split('\t')
            except ValueError:
                chrm, pos, ident, ref, alt, qual, filt, info, form, sampleGen = line.strip().split('\t')
            
            sampleGenList = sampleGen.split(':')
            if len(sampleGenList) == 1:
                sampleGensplit = sampleGenList[0].split(',')
                refCount = int(sampleGensplit[0])
                altCount = int(sampleGensplit[1])
                totalCount = refCount + altCount
                
                if altCount / totalCount >= 0.95:
                    genotype = 'HOM'
                else:
                    genotype = 'HET'
            elif len(sampleGenList) == 2:
                sampleGensplit = sampleGenList[0].split(',')
                refCount1 = int(sampleGensplit[0])
                altCount1 = int(sampleGensplit[1])
                totalCount1 = refCount1 + altCount1
                if totalCount1 == 0:
                    genotype1 = '0'
                elif altCount1 / totalCount1 >= 0.95:
                    genotype1 = 'HOM'
                else:
                    genotype1 = 'HET'
                sampleGensplit = sampleGenList[1].split(',')
                refCount2 = int(sampleGensplit[0])
                altCount2 = int(sampleGensplit[1])
                totalCount2 = refCount2 + altCount2
                if totalCount2 == 0:
                    genotype2 = '0'
                elif altCount2 / totalCount2 >= 0.95:
                    genotype2 = 'HOM'
                else:
                    genotype2 = 'HET'
                genotype = 'HET'
                if genotype1 == 'HOM':
                    if genotype2 != 'HET':
                        genotype = 'HOM'
                    
                elif genotype2 == 'HOM':
                    if genotype1 != 'HET':
                        genotype = 'HOM'
            elif len(sampleGenList) > 2:
                #Leiden
                sampleGensplit = sampleGenList[-2].split(',')
                refCount1 = int(sampleGensplit[0])
                altCount1 = int(sampleGensplit[1])
                totalCount1 = refCount1 + altCount1
                if totalCount1 == 0:
                    genotype1 = '0'
                elif altCount1 / totalCount1 >= 0.95:
                    genotype1 = 'HOM'
                else:
                    genotype1 = 'HET'
                sampleGensplit = sampleGenList[-1].split(',')
                refCount2 = int(sampleGensplit[0])
                altCount2 = int(sampleGensplit[1])
                totalCount2 = refCount2 + altCount2
                if totalCount2 == 0:
                    genotype2 = '0'
                elif altCount2 / totalCount2 >= 0.95:
                    genotype2 = 'HOM'
                else:
                    genotype2 = 'HET'
                genotype = 'HET'
                if genotype1 == 'HOM':
                    if genotype2 != 'HET':
                        genotype = 'HOM'
                    
                elif genotype2 == 'HOM':
                    if genotype1 != 'HET':
                        genotype = 'HOM'
                   
            
            if chrm == 'X':
                chrm = 0
            elif chrm == 'Y':
                chrm = 42
            chrm = int(chrm)
            pos = int(pos)
            realId = ident.split(':')[0].replace('Manta', '')
            if realId == 'BND':
                start = int(pos)
                end = '.'               
                
            else:
                start = pos
                infoList = info.split(';')
                for el in infoList:
                    if 'END=' in el:
                        end = int(el.replace('END=', ''))
                        
                        break
            if not chrm in svWhole:
                svWhole[chrm] = [(start, end, realId, alt, genotype)]
            else:
                svWhole[chrm].append((start, end, realId, alt, genotype))
                
    geneSV = {}
    for gene, coords in geneDict.items():
        geneChr, geneStart, geneEnd = (coords[0], coords[1], coords[2])
        if geneChr not in svWhole:
            geneSV[gene] = ''
            continue
        for sv in svWhole[geneChr]:            
            svStart, svEnd, svID, svAlt, svGen = sv
            if svID == 'BND':
                ##DO SOMETHING
                continue
            else:
                if geneEnd > svStart and geneStart < svEnd:
                    geneSV[gene] = sv
                    break
        if not gene in geneSV:
            geneSV[gene] = ''
    return geneSV
        
          
          
def snvFilter(snvDict):  
    #filtered version
    annotationRank = ["chromosome_number_variation", 'exon_loss_variant', 'frameshift_variant', 'stop_gained', 'stop_lost', 'start_lost', 
             'splice_acceptor_variant', 'splice_donor_variant', 'rare_amino_acid_variant', 'missense_variant', 'inframe_insertion',
             'disruptive_inframe_insertion', 'inframe_deletion', 'disruptive_inframe_deletion', '5_prime_UTR_truncation+exon_loss_variant',
             '3_prime_UTR_truncation+exon_loss', 'splice_branch_variant', 'splice_region_variant', 'splice_branch_variant', 
             'stop_retained_variant', 'initiator_codon_variant', 'synonymous_variant', 'initiator_codon_variant+non_canonical_start_codon', 
             'coding_sequence_variant', '5_prime_UTR_variant', '3_prime_UTR_variant', '5_prime_UTR_premature_start_codon_gain_variant',
             'upstream_gene_variant', 'downstream_gene_variant', 'TF_binding_site_variant', 'regulatory_region_variant', 'miRNA']
    
    filteredGeneDict = {}    
    for gene, snvList in snvDict.items():
        annotations = []
        for snv in snvList:
            try:
                aList = []                
                a = snv['info']['ANN'][0]
                aList.append((a.strip().split('|')[1], snv['genotype']))
                annotations += aList
            except KeyError:
                continue
        
        setAnnotations = annotations       
        sortedAnnotations = []
        for a in setAnnotations:
            if a[0] not in annotationRank:
                continue
            else:
                sortedAnnotations.append((a, annotationRank.index(a[0])))
        
        if sortedAnnotations == []:            
            continue
        
        sortedAnnotations.sort(key=lambda x: x[1])
        
        if len(sortedAnnotations) > 1:
            snv1 = sortedAnnotations[0][0]
            snv2 = sortedAnnotations[1][0]
        elif len(sortedAnnotations) == 1:
            snv1 = sortedAnnotations[0][0]
            snv2 = ('.', '.')
        filteredGeneDict[gene] = (snv1, snv2)
    return filteredGeneDict
          
          
          
def printVar(geneDict, snvDict, svDict, cnvDict, sample, patient):
    ##SAMPLE GENE CHR snv1 snv2 sv cn cntype
    
    filteredSnv = snvFilter(snvDict)
    
    for gene, coords in geneDict.items():
        chrm = coords[0]
        try:
            snv1 = filteredSnv[gene][0][0]
            gen1 = filteredSnv[gene][0][1]
            snv2 = filteredSnv[gene][1][0]
            gen2 = filteredSnv[gene][1][1]
        except KeyError:
            snv1 = '.'
            gen1 = '.'
            snv2 = '.'
            gen2 = '.'
        if svDict[gene] != '':
            sv = svDict[gene][2]
            svgen = svDict[gene][4]
        else:
            sv = '.'
            svgen = '.'
        cn, cnType = cnvDict[gene][2], cnvDict[gene][3]
        
    
        print '\t'.join([patient, sample, gene, str(chrm), snv1, gen1, snv2, gen2, sv, svgen, str(cn), cnType])    
            

                
            
        


###MAIN BODY###

def main():
    args=arguments()
  
       
    geneDict = genepanel()
          

    snvDict = snvReader(args.snv, geneDict)

    svDict = svReader(args.sv, geneDict)
    cnvDict = cnvReader(args.cnv, geneDict)
    sample = args.sample
    patient = args.patient
    printVar(geneDict, snvDict, svDict, cnvDict, sample, patient)
    
    

      
   
###Execute main body

if __name__ == '__main__': main()