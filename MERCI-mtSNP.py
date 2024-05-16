
# coding: utf-8

# In[ ]:

import os
import pysam
import pandas as pd
import numpy as np
from functools import partial
import time
from optparse import OptionParser
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


# In[ ]:

##########################################################################################
# load MT reference genome
##########################################################################################
def load_MTgenome(path_fa):
    ref_fa = pysam.FastaFile(path_fa)
    if 'MT' in ref_fa.references:
        MT_ref = ref_fa.fetch("MT")
    if 'chrM' in ref_fa.references:
        MT_ref = ref_fa.fetch("chrM")
    if 'MT' not in ref_fa.references and 'chrM'not in ref_fa.references:
        print('error can not load MT reference, because no \'MT\' or \'chrM\' reference name are found!')
    MT_ref = MT_ref.upper()
    return(MT_ref)


# In[ ]:

##########################################################################################
# load selected cell barcodes for 10x data
##########################################################################################
def load_CBs(path_barcodes):
    os.chdir(path_barcodes)
    if('barcodes.tsv' not in os.listdir() ):
        os.system('gunzip -dc barcodes.tsv.gz>barcodes.tsv')
    filename = path_barcodes +  "/" + 'barcodes.tsv'
    with open(filename, "r") as f:
        CBs = f.read()
    CBs = CBs.strip('').split('\n') 
    return(CBs)


# In[ ]:

##########################################################################################
#write the reads of selected cells mapping to MT genome (time～20 mins), 10x scRNA-seq
##########################################################################################
def WrteSbamfile_10xRNA(path_out, bamfile, CBs, sampleID='S0010', Mquality=255, MTref_name='MT'):  
    filename = path_out + "/" + sampleID + '.MT.bam'
    Sfile = pysam.AlignmentFile(filename, 'wb', template=bamfile)
    i = 0
    s_CBs = set()
    for read in bamfile.fetch(MTref_name):
        if (not read.is_duplicate) and (not read.is_secondary) and (not read.is_supplementary) and (not read.is_unmapped):
            if read.mapping_quality>=Mquality:
                if read.has_tag('CB'):
                    t_CB = read.get_tag(tag='CB')
                    if t_CB in CBs:
                        Sfile.write(read)
                        i += 1
                        s_CBs.add(t_CB)                    
    Sfile.close()
    pysam.index(filename)
    CB_number = len(s_CBs)
    print('A total of %s reads were written out'% i)
    print('%s cell barcodes have unqiue reads mapping to MT'% CB_number)
   
##########################################################################################
#write the reads of selected cells mapping to MT genome, BD Rhapsody
##########################################################################################
def WrteSbamfile_BDRh(path_out, bamfile, CBs, sampleID='S0010', Mquality=255, MTref_name='MT'):  
    filename = path_out + "/" + sampleID + '.MT.bam'
    Sfile = pysam.AlignmentFile(filename, 'wb', template=bamfile)
    i = 0
    s_CBs = set()
    for read in bamfile.fetch(MTref_name):
        if (not read.is_duplicate) and (not read.is_secondary) and (not read.is_supplementary) and (not read.is_unmapped):
            if read.mapping_quality>=Mquality:
                if read.has_tag('CB'):
                    t_CB = read.get_tag(tag='CB')
                    t_CB=t_CB.replace('_', '')
                    if t_CB in CBs:
                        Sfile.write(read)
                        i += 1
                        s_CBs.add(t_CB)                    
    Sfile.close()
    pysam.index(filename)
    CB_number = len(s_CBs)
    print('A total of %s reads were written out'% i)
    print('%s cell barcodes have unqiue reads mapping to MT'% CB_number)

##########################################################################################
#write the high quality reads mapping to MT genome, 10x scATAC-seq
##########################################################################################
def WrteSbamfile_10xATAC(path_out, bamfile, CBs, sampleID='S0010', Mquality=5, MTref_name='ChrM'):  
    filename = path_out + "/" + sampleID + ".MT.bam"
    Sfile = pysam.AlignmentFile(filename, 'wb', template=bamfile)
    i = 0
    s_CBs = set()
    for read in bamfile.fetch(MTref_name):
        if (not read.is_duplicate) and (not read.is_secondary) and (not read.is_supplementary) and (not read.is_unmapped):
            if read.mapping_quality >= Mquality:
                if read.has_tag('CB'):
                    t_CB = read.get_tag(tag='CB')
                    if t_CB in CBs:
                        Sfile.write(read)
                        i += 1
                        s_CBs.add(t_CB)     
    Sfile.close()
    pysam.index(filename)
    CB_number = len(s_CBs)
    print('A total of %s reads were written out'% i)
    print('%s cell barcodes have unqiue reads mapping to MT'% CB_number)
    
#######################################################################################
#write the high quality reads mapping to MT genome, 10x scATAC-seq for mouse genome, retain part of reads with MQ=0
#######################################################################################
def WrteSbamfile_10xATAC_mouse(path_out, bamfile, CBs, sampleID='S0010', Mquality=5, MTref_name='ChrM'):  
    filename = path_out + "/" + sampleID + ".MT.bam"
    Sfile = pysam.AlignmentFile(filename, 'wb', template=bamfile)
    i = 0
    s_CBs = set()
    for read in bamfile.fetch(MTref_name):
        if (not read.is_duplicate) and (not read.is_secondary) and (not read.is_supplementary) and (not read.is_unmapped):
            read_XA = False
            if read.has_tag('XA'):
                tagXA = read.get_tag(tag='XA')
                if ('chr1,+246' in tagXA )| ('chr1,-246' in tagXA ):
                     read_XA = True
            if (read.mapping_quality >= Mquality) | read_XA:
                if read.has_tag('CB'):
                    t_CB = read.get_tag(tag='CB')
                    if t_CB in CBs:
                        Sfile.write(read)
                        i += 1
                        s_CBs.add(t_CB)     
    Sfile.close()
    pysam.index(filename)
    CB_number = len(s_CBs)
    print('A total of %s reads were written out'% i)
    print('%s cell barcodes have unqiue reads mapping to MT'% CB_number)

##########################################################################################
#write the high quality reads mapping to MT genome, smart-seq2 or RNA-seq
##########################################################################################
def WrteSbamfile_RNAseq(path_out, bamfile, sampleID='S0010', Mquality=255, MTref_name='MT'):  
    filename = path_out + "/" + sampleID + ".MT.bam"
    Sfile = pysam.AlignmentFile(filename, 'wb', template=bamfile)
    i = 0
    for read in bamfile.fetch(MTref_name):
        if (not read.is_duplicate) and (not read.is_secondary) and (not read.is_supplementary) and (not read.is_unmapped):
            if read.mapping_quality >= Mquality:
                Sfile.write(read)
                i += 1
    Sfile.close()
    pysam.index(filename)
    print('A total of %s reads were written out'% i)
    
##########################################################################################
#write the high quality reads mapping to MT genome, ATATseq or scATAC-seq
##########################################################################################
def WrteSbamfile_ATAC(path_out, bamfile, sampleID='S0010', Mquality=5, MTref_name='chrM'):  
    filename = path_out + "/" + sampleID + ".MT.bam"
    Sfile = pysam.AlignmentFile(filename, 'wb', template=bamfile)
    i = 0
    for read in bamfile.fetch(MTref_name):
        if (not read.is_duplicate) and (not read.is_secondary) and (not read.is_supplementary) and (not read.is_unmapped):
            if read.mapping_quality >= Mquality:
                Sfile.write(read)
                i += 1
    Sfile.close()
    pysam.index(filename)
    print('A total of %s reads were written out'% i)
    
#############################################################################################
#WrteSbamfile function
#############################################ATAC#############################################
def WrteSbamfile(path_out, bamfile, sampleID, Mquality, MTref_name, CBs=None, dataType='10x_scRNA-seq', sp='human'):
    if dataType=='10x_scRNA-seq':
        WrteSbamfile_10xRNA(path_out=path_out, bamfile=bamfile, CBs=CBs, sampleID=sampleID, Mquality=Mquality, MTref_name=MTref_name)
    if (dataType=='10x_mtscATAC-seq') and (sp=='human'):
        WrteSbamfile_10xATAC(path_out=path_out, bamfile=bamfile, CBs=CBs, sampleID=sampleID, Mquality=Mquality, MTref_name=MTref_name)
    if (dataType=='10x_mtscATAC-seq') and (sp=='mouse'):
        WrteSbamfile_10xATAC_mouse(path_out=path_out, bamfile=bamfile, CBs=CBs, sampleID=sampleID, Mquality=Mquality, MTref_name=MTref_name)
    if dataType in ['smart-seq2', 'bulk_RNA-seq']:
        WrteSbamfile_RNAseq(path_out=path_out, bamfile=bamfile, sampleID=sampleID, Mquality=Mquality, MTref_name=MTref_name)
    if dataType in ['bulk_ATAC-seq', 'scATAC-seq']:
        WrteSbamfile_ATAC(path_out=path_out, bamfile=bamfile, sampleID=sampleID, Mquality=Mquality, MTref_name=MTref_name)
    if dataType=='BD-Rhapsody_scRNA-seq':
        WrteSbamfile_BDRh(path_out=path_out, bamfile=bamfile,  CBs=CBs, sampleID=sampleID, Mquality=Mquality, MTref_name=MTref_name)
# In[ ]:

##########################################################################################
#Return the read List in each cell, 10x and RD-Rhapsody
##########################################################################################
def ReadCounts_cell(Sfile, MTref_name='MT', dataType='10x_scRNA-seq'):
    Cell_reads = {}
    for read in Sfile.fetch(MTref_name):
        t_CB = read.get_tag(tag='CB')
        if dataType=='BD-Rhapsody_scRNA-seq':
            t_CB=t_CB.replace('_', '')
        if t_CB not in Cell_reads:
            Cell_reads[t_CB] = [read]
        else:
            Cell_reads[t_CB].append(read) 
    return(Cell_reads)

# In[ ]:

##########################################################################################
#Return the read Number in each cell, 10x and RD-Rhapsody
##########################################################################################
def ReadNumber_cell(Cell_reads):
    Read_Num = pd.Series(0, index=list(Cell_reads.keys()))
    for CB_, reads in Cell_reads.items():
        tmp = len(reads)
        Read_Num[CB_] = tmp
    return(Read_Num)


# In[ ]:

##########################################################################################
#Return the read Coverage in each cell, 10x and RD-Rhapsody
##########################################################################################
def ReadCoverage_cell(Cell_reads, MT_ref, Qcutoff=15):
    Coverage_Cell = pd.DataFrame(columns=list(Cell_reads.keys()))
    for CB_, readList_ in Cell_reads.items():
        Coverage_ = np.zeros(len(MT_ref), dtype='i')
        for read in readList_:
            tmp = [(i, j) for i, j in read.get_aligned_pairs() if j in read.get_reference_positions() ]
            overlap_pos =[j for i, j in tmp if read.query_qualities[i]>=Qcutoff]
            Coverage_[overlap_pos] += 1
        Coverage_Cell[CB_] = Coverage_
    return(Coverage_Cell)
##########################################################################################
#Return the read Coverage for each MT site, smart-seq2 and bulk RNA-seq
##########################################################################################
def writeCoverage_RNA(Sfile, path_out, sampleID, Qcutoff=20, MTref_name='MT'):
    MT_total_coverage = Sfile.count_coverage(MTref_name, quality_threshold=Qcutoff)
    MT_Coverage = pd.DataFrame(columns=["A", "C", "G", "T"], dtype='int64')
    MT_Coverage["A"]=MT_total_coverage[0]
    MT_Coverage["C"]=MT_total_coverage[1]
    MT_Coverage["G"]=MT_total_coverage[2]
    MT_Coverage["T"]=MT_total_coverage[3]
    os.chdir(path_out)
    outname = sampleID + '.MT_Coverage.csv'
    MT_Coverage.to_csv(outname)
##########################################################################################
#Return the read Coverage for each MT site, ATAC-seq and scATAC-seq
##########################################################################################
def writeCoverage_ATAC(Sfile, path_out, sampleID, Qcutoff=25, MTref_name='chrM'):
    MT_total_coverage = Sfile.count_coverage(MTref_name, quality_threshold=Qcutoff)
    MT_Coverage = pd.DataFrame(columns=["A", "C", "G", "T"], dtype='int64')
    MT_Coverage["A"]=MT_total_coverage[0]
    MT_Coverage["C"]=MT_total_coverage[1]
    MT_Coverage["G"]=MT_total_coverage[2]
    MT_Coverage["T"]=MT_total_coverage[3]
    os.chdir(path_out)
    outname = sampleID + '.MT_Coverage.csv'
    MT_Coverage.to_csv(outname)


# In[ ]:

##########################################################################################
#Return the covered genome fraction, 10x and RD-Rhapsody
##########################################################################################
def funct(x, length, cutoff=1):
    y = sum(pd.Series(x)>=cutoff)
    frac = y/length
    return(frac)

def CovergeF(Coverage_Cell, path_out, sampleID, lengthMT, minC=1):
    funct2 = partial(funct, length=lengthMT, cutoff=minC)
    coverageF_stat = Coverage_Cell.apply(funct2, axis=0)
    plt.hist(coverageF_stat, range = (0, 1), bins = 10, edgecolor='black')
    plt.title("Reads covered MT fraction")
    xlable = "The fraction of MT genome covered by at least " + str(minC) + "x"
    plt.xlabel(xlable)
    plt.ylabel("Cell counts")
    os.chdir(path_out)
    plt.savefig(sampleID + ".Coverage.pdf")
    print(np.mean(coverageF_stat))
    return(coverageF_stat)
##########################################################################################
#Total Coverage, smart-seq2 and bulk RNA-seq
##########################################################################################
def Coverage_stat_RNA(Sfile, path_out, sampleID, Qcutoff=25, minC=1, MTref_name='MT'):
    MT_total_coverage = Sfile.count_coverage(MTref_name, quality_threshold=Qcutoff) 
    A_coverage = MT_total_coverage[0]
    C_coverage = MT_total_coverage[1]
    G_coverage = MT_total_coverage[2]
    T_coverage = MT_total_coverage[3]
    b = np.array(A_coverage) + np.array(C_coverage) + np.array(G_coverage) + np.array(T_coverage)
    x = [i for i in range(len(b))]
    plt.plot(x, b, color='blue', linestyle='solid') ;
    plt.xlabel('MT genome coordinate')
    plt.ylabel('Site coverage')
    frac = np.count_nonzero(b>=minC)/len(b)
    Title_text = 'Coverage\nmedian coverage'+ str(np.median(b))+ ';\n'+ str(frac*100) + '% MT genome were covered by at least' + str(minC) + ' reads' 
    plt.title(Title_text)
    os.chdir(path_out)
    plt.savefig(sampleID + ".Coverage.pdf")
##########################################################################################
#Total Coverage, scATAC-seq and ATAC-seq
##########################################################################################
def Coverage_stat_ATAC(Sfile, path_out, sampleID, Qcutoff=25, minC=1, MTref_name='chrM'):
    MT_total_coverage = Sfile.count_coverage(MTref_name, quality_threshold=Qcutoff) 
    A_coverage = MT_total_coverage[0]
    C_coverage = MT_total_coverage[1]
    G_coverage = MT_total_coverage[2]
    T_coverage = MT_total_coverage[3]
    b = np.array(A_coverage) + np.array(C_coverage) + np.array(G_coverage) + np.array(T_coverage)
    x = [i for i in range(len(b))]
    plt.plot(x, b, color='blue', linestyle='solid') ;
    plt.xlabel('MT genome coordinate')
    plt.ylabel('Site coverage')
    frac = np.count_nonzero(b>=minC)/len(b)
    Title_text = 'Coverage\nmedian coverage'+ str(np.median(b))+ ';\n'+ str(frac*100) + '% MT genome were covered by at least' + str(minC) + ' reads' 
    plt.title(Title_text)
    os.chdir(path_out)
    plt.savefig(sampleID + ".Coverage.pdf")


# In[ ]:

##########################################################################################
#Return the allele frequency of SNPs for each cell, 10x and RD-Rhapsody
##########################################################################################
def AF_cal(Cell_mutations):
    for t_Cell, SNP_record_ in Cell_mutations.items():
        SNP_records = SNP_record_
        for mutID, info in SNP_record_.items():
            record = info
            AF = info[1]/info[2]
            record.append(AF)
            SNP_records[mutID] = record
        Cell_mutations[t_Cell] = SNP_records
    return(Cell_mutations)
##########################################################################################
#Return the SNPs for each cell, 10x and RD-Rhapsody
##########################################################################################
def SNP_caller_10x(Cell_reads, Coverage_Cell, MT_ref, outpath, sampleID, Qcutoff=15):
    os.chdir(outpath)
    Cell_mutations = {}
    for t_Cell, t_readList in Cell_reads.items():
        SNP_records = {}
        for read in t_readList:
            aligned_pairs = read.get_aligned_pairs()
            for i, j in aligned_pairs:
                if isinstance(i, int) and isinstance(j,int):
                    if read.query_sequence[i]==MT_ref[j] or read.query_qualities[i] < Qcutoff :
                        pass
                    else:
                        SNP = "MT_" + str(j) + "_" + MT_ref[j] + "-" + read.query_sequence[i]
                        if SNP not in SNP_records:
                            t_coverage = Coverage_Cell[t_Cell][j]
                            record = [j, 1, t_coverage]  #[pos, s_reads, t_reads]
                            SNP_records[SNP] = record
                        else:
                            SNP_records[SNP][1] += 1
        Cell_mutations[t_Cell] = SNP_records
    Cell_mutations = AF_cal(Cell_mutations)
    filename = sampleID + '.Cell_mutations'
    np.save(filename, Cell_mutations)
    return(Cell_mutations)
##########################################################################################
#SNP clusters fiteration
#filter the some MTmutations for each cell, those SNPs of only 1 s_read and co-locate with other SNPs in ln length will be removed
##########################################################################################
def MTmutations_filters(Cell_mutations, sampleID, ln=5):
    New_Cell_mutations = Cell_mutations
    for Cell_, SNP_records_ in Cell_mutations.items():
        tmp_variants = list(SNP_records_.keys())
        tmp_records = list(SNP_records_.values())
        positions = np.array([x[0] for x in tmp_records])
        t_index = []
        for pos in positions:
            t_max = pos + ln
            t_min = pos - ln
            func = lambda x: x <= t_max and x >= t_min
            t_stat = sum( list(map(func, positions)) )
            t_index.append(t_stat > 1)
        if sum(t_index) == 0:
            continue
        c_variants = np.array(tmp_variants)[t_index]
        c_records = np.array(tmp_records)[t_index]
        s_reads = c_records[:, 1] 
        t_index2 = s_reads == 1
        if sum(t_index2) == 0:
            continue
        r_variants = c_variants[t_index2]
        New_SNP_records_ = SNP_records_
        for var in r_variants:
            del New_SNP_records_[var]
        New_Cell_mutations[Cell_] = New_SNP_records_
    filename = sampleID + '.Filterd.Cell_mutations'
    np.save(filename, New_Cell_mutations)
    return New_Cell_mutations
##########################################################################################
#Return the MT variants table, smart-seq2 and bulk RNA-seq
##########################################################################################
def SNP_caller_RNA(Sfile, MT_ref, sampleID,  path_out, Qcutoff=25, MTref_name='MT'):
    os.chdir(path_out)
    Site_mutations = {}
    for pileupColum in Sfile.pileup(MTref_name, max_depth=100000, min_base_quality=Qcutoff, stepper='nofilter', flag_filter=0, ignore_overlaps=False):
        t_pos = pileupColum.reference_pos
        t_coverage = 0
        Ref_pos = pileupColum.reference_pos
        Ref_base = MT_ref[Ref_pos]
        SNP_IDs = []
        baseQ = {}
        for pileupread in pileupColum.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                Read = pileupread.alignment
                Read_pos = pileupread.query_position
                Read_base =  Read.query_sequence[Read_pos]
                if Read_base == Ref_base:
                    t_coverage += 1
                else:
                    ID = "MT_" + str(t_pos) + "_" + Ref_base + "-" + Read_base
                    SNP_IDs.append(ID)
                    if ID not in baseQ:
                        baseQ[ID] = [Read.query_qualities[Read_pos]]
                    else:
                        baseQ[ID].append(Read.query_qualities[Read_pos]) 
                    t_coverage += 1
        SNP_sat = pd.value_counts(SNP_IDs)
        IDs = list(SNP_sat.index)
        for t_ID in IDs:
            pos= t_pos
            s_reads = SNP_sat[t_ID]
            avg_BQ = np.mean(baseQ[t_ID])
            AF = s_reads/t_coverage
            record = [pos, s_reads, avg_BQ, t_coverage, AF]
            Site_mutations[t_ID] = record
    filename = sampleID +'.Site_mutations'        
    np.save(filename, Site_mutations)
    return(Site_mutations)
##########################################################################################
#Return the MT variants table, scATAC-seq and bulk ATAC-seq
##########################################################################################
def SNP_caller_ATAC(Sfile, MT_ref, sampleID, path_out, Qcutoff=25, MTref_name='chrM'):
    os.chdir(path_out)
    Site_mutations = {}
    for pileupColum in Sfile.pileup(MTref_name, max_depth=100000, min_base_quality=Qcutoff, stepper='nofilter', flag_filter=0, ignore_overlaps=False):
        t_pos = pileupColum.reference_pos
        t_coverage = 0
        Ref_pos = pileupColum.reference_pos
        Ref_base = MT_ref[Ref_pos]
        SNP_IDs = []
        baseQ = {}
        for pileupread in pileupColum.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                Read = pileupread.alignment
                Read_pos = pileupread.query_position
                Read_base =  Read.query_sequence[Read_pos]
                if Read_base == Ref_base:
                    t_coverage += 1
                else:
                    ID = "MT_" + str(t_pos) + "_" + Ref_base + "-" + Read_base
                    SNP_IDs.append(ID)
                    if ID not in baseQ:
                        baseQ[ID] = [Read.query_qualities[Read_pos]]
                    else:
                        baseQ[ID].append(Read.query_qualities[Read_pos]) 
                    t_coverage += 1
        SNP_sat = pd.value_counts(SNP_IDs)
        IDs = list(SNP_sat.index)
        for t_ID in IDs:
            pos= t_pos
            s_reads = SNP_sat[t_ID]
            avg_BQ = np.mean(baseQ[t_ID])
            AF = s_reads/t_coverage
            record = [pos, s_reads, avg_BQ, t_coverage, AF]
            Site_mutations[t_ID] = record
    filename = sampleID +'.Site_mutations'        
    np.save(filename, Site_mutations)
    return(Site_mutations)


# In[ ]:

##########################################################################################
#Write SNPs into txt， 10x and RD-Rhapsody
##########################################################################################
def SNP_to_txt_10x(Cell_mutations, ReadNumber, coverageF_stat1X, sampleID, outpath):
    os.chdir(outpath)
    filename = sampleID+'.MT_variants.txt'
    fileObject = open(filename, 'w')
    fileObject.write("Cell\t" + 'ID\t'+'pos\t'+'s_reads\t'+'t_reads\t'+'AF\t'+'Cell_reads\t'+ 'CovFraction_minCx' +'\n')
    for Cell_, SNP_records_ in Cell_mutations.items():
        CB = Cell_
        CellReads = ReadNumber[Cell_]
        Covered_MT = coverageF_stat1X[Cell_]
        for SNP, record in SNP_records_.items():
            mutID = SNP
            pos = record[0]
            support_reads = record[1]
            total_reads = record[2]
            AF = record[3]
            fileObject.write(CB + '\t')
            fileObject.write(mutID + '\t')
            fileObject.write(str(pos) + '\t')
            fileObject.write(str(support_reads) + '\t')
            fileObject.write(str(total_reads) + '\t')
            fileObject.write(str(AF) + '\t')
            fileObject.write(str(CellReads) + '\t')
            fileObject.write(str(Covered_MT))
            fileObject.write('\n')
    fileObject.close()
##########################################################################################
#Write SNPs into txt, smart-seq2/bulk RNA-seq/scATAC-seq/ATAC-seq
##########################################################################################
def SNP_to_txt(Site_mutations, sampleID, path_out):
    os.chdir(path_out)
    filename = sampleID +'.MT_variants.txt'
    fileObject = open(filename, 'w')
    fileObject.write('ID\t'+'pos\t'+'s_reads\t'+'avg_BQ\t' +'t_reads\t'+'AF'+'\n')
    for ID_, SNP_records_ in Site_mutations.items():
        mutID = ID_
        pos = SNP_records_[0]
        support_reads = SNP_records_[1]
        avg_BQ = SNP_records_[2]
        total_reads = SNP_records_[3]
        AF = SNP_records_[4]
        fileObject.write(mutID + '\t')
        fileObject.write(str(pos) + '\t')
        fileObject.write(str(support_reads) + '\t')
        fileObject.write(str(avg_BQ) + '\t')
        fileObject.write(str(total_reads) + '\t')
        fileObject.write(str(AF))
        fileObject.write('\n')
    fileObject.close()


# In[ ]:

##########################################################################################
#parameters setting
##########################################################################################
def CommandLineParser():
    usage = "python MERCI-mtSNP.py [-D <dataType>] [-o <Directory>] [-S <sampleID>] [-b <path_bam>] [-f <path_fa>] [-c <path_barcodes>] [-M <Mquality>] [-B <Qcutoff>] [-r <Species>] [-l <ln>] [-m <minC>]"
    parser=OptionParser(usage)
    print ('''
    MERCI-mtSNP version 1.5.0!
    mtSNP calling start...
    parameters setting...
    ''')
    parser.add_option("-D","--dataType", action="store", dest="dataType", default='10x_scRNA-seq', help="The data type of your sequencing data. One of '10x_scRNA-seq'(default), 'smart-seq2', 'bulk_ATAC-seq', 'scATAC-seq', 'bulk_RNA-seq', 'BD-Rhapsody_scRNA-seq', or '10x_mtscATAC-seq'")
    parser.add_option("-o","--output", action="store", dest="Directory", default='./', help="Output directory for intermediate and final outputs.")
    parser.add_option("-S","--sampleID", action="store", dest="sampleID", default='sampleX', help="the sample name, also serve as the name of output file. if not given, the names of all intermeidate or final output files will be automatically set as sampleX")
    parser.add_option("-b","--Bamfile", action="store", dest="path_bam", default='', help="Input bam file for MT mutation calling")
    parser.add_option("-f","--fastafile", action="store", dest="path_fa", default='', help="The genome reference sequence in fasta format, usually named as genome.fa")
    parser.add_option("-c","--CellBarcode", action="store", dest="path_barcodes", default='None', help="This parameter only work for dataTypes with 10x_scRNA-seq, BD-Rhapsody_scRNA-seq or 10x_mtscATAC-seq, the directory where cell barcodes file (barcodes.tsv.gz or barcodes.tsv) generated by cellranger exists")
    parser.add_option("-M","--MQcutoff", action="store", type="int", dest="Mquality", default='255', help="The lowest alignment quality that are accepted, the reads with alignment scores below the given value will be discarded, default=5 for scATAC-seq, 10x_mtscATAC-seq or bulk_ATAC-seq, default=255 for other dataTypes")    
    parser.add_option("-B","--BQcutoff", action="store", type="int", dest="Qcutoff", default='15', help="The base qulaity cutoff, only alleles with BQ higher than this value will be retained, default=15 for 10x_scRNA-seq and BD-Rhapsody_scRNA-seq, default=25 for other dataTypes")
    parser.add_option("-r","--ref", action="store", dest="Species", default='human', help="This parameter only works for 10x_mtscATAC-seq dataType, user can set 'human' or 'mouse' depending on what species the sequencing data is, default=mouse")
    parser.add_option("-l","--ln", action="store", type="int", dest="ln", default='5', help="This parameter only works for 10x_scRNA-seq and BD-Rhapsody_scRNA-seq dataTypes, the maximum sequence range of snp clusters, reads supporting multiple variants within a small genomic region (ln bp) will be reomved, default=5")
    parser.add_option("-m","--minC", action="store", type="int", dest="minC", default='1', help="This parameter works for all data types expcept those of 10x and BD-Rhapsody platforms, A threshold for coverage, the faction of MT genome that was covered by read counts no less than than this value will be recorded on the generated coverage figure, default=1")
    return parser.parse_args()


# In[ ]:

def main():
    opt, args = CommandLineParser()
    dataType = opt.dataType
    sampleID = opt.sampleID
    path_bam = opt.path_bam
    path_fa = opt.path_fa
    path_barcodes = opt.path_barcodes
    path_out = opt.Directory
    Mquality = opt.Mquality
    Qcutoff = opt.Qcutoff
    sp = opt.Species
    ln = opt.ln
    minC = opt.minC
    
    if dataType not in ['10x_scRNA-seq', 'smart-seq2', 'bulk_ATAC-seq', 'scATAC-seq', 'bulk_RNA-seq', '10x_mtscATAC-seq', 'BD-Rhapsody_scRNA-seq']:
        print('Warning: dataType is not assigned or not within supported category')
        
    if dataType in ['scATAC-seq', 'bulk_ATAC-seq', '10x_mtscATAC-seq'] and Mquality==255:
        Mquality=5
    if dataType not in ['10x_scRNA-seq', 'BD-Rhapsody_scRNA-seq'] and Qcutoff==15:
        Qcutoff=25
    
    bamfile = pysam.AlignmentFile(path_bam, "rb")
    #loading the MT reference genome sequence, 16,569 for human, 16,299 for mouse
    MT_ref = load_MTgenome(path_fa)
    if dataType in ['10x_scRNA-seq', '10x_mtscATAC-seq', 'BD-Rhapsody_scRNA-seq']:
        CBs = load_CBs(path_barcodes)
    else:
        CBs = None
    
    allref_names = bamfile.references
    if 'MT' in allref_names:
        MTref_name = "MT"
    if 'chrM' in allref_names:
        MTref_name = "chrM"
    if 'MT' not in allref_names and 'chrM'not in allref_names:
        print('error can not load MT reference, because no \'MT\' or \'chrM\' reference name are found!')
    
    #write out the qualified MT reads in to new bam file
    WrteSbamfile(path_out, bamfile, sampleID,  Mquality, MTref_name, CBs=CBs, dataType=dataType, sp=sp)

    path_bam2 =  path_out + "/" + sampleID + ".MT.bam"
    Sfile = pysam.AlignmentFile(path_bam2, 'rb')
    
    #the read List in each cell, 10x rna-seq, BD-Rhapsody_scRNA-seq or 10x mtscATAC-seq
    if dataType in ['10x_scRNA-seq', '10x_mtscATAC-seq', 'BD-Rhapsody_scRNA-seq']:
        Cell_reads = ReadCounts_cell(Sfile, MTref_name=MTref_name, dataType=dataType)
        Read_Num = ReadNumber_cell(Cell_reads)
        print('The median MT read count is %s per cell'% str( np.median(Read_Num) ) )
    #read Coverage in each cell
    if dataType in ['10x_scRNA-seq', '10x_mtscATAC-seq', 'BD-Rhapsody_scRNA-seq']:
        Coverage_Cell = ReadCoverage_cell(Cell_reads, MT_ref, Qcutoff=Qcutoff)
        os.chdir(path_out)
        filename= sampleID + '.Coverage_Cell'
        Coverage_Cell.to_csv( filename + ".csv")
    if dataType in ['smart-seq2', 'bulk_RNA-seq']:
        writeCoverage_RNA(Sfile, path_out, sampleID, Qcutoff=Qcutoff, MTref_name=MTref_name)
    if dataType in ['scATAC-seq', 'bulk_ATAC-seq']:
        writeCoverage_ATAC(Sfile, path_out, sampleID, Qcutoff=Qcutoff, MTref_name=MTref_name)
 
    #the covered genome fraction at 1x (or minC x)
    if dataType in ['10x_scRNA-seq', '10x_mtscATAC-seq', 'BD-Rhapsody_scRNA-seq']:
        lengthMT = len(MT_ref)
        coverageF_stat = CovergeF(Coverage_Cell, path_out, sampleID, lengthMT, minC=minC)
    if dataType in ['smart-seq2', 'bulk_RNA-seq']:
        Coverage_stat_RNA(Sfile, path_out, sampleID, Qcutoff, minC, MTref_name=MTref_name)
    if dataType in ['scATAC-seq', 'bulk_ATAC-seq']:
        Coverage_stat_ATAC(Sfile, path_out, sampleID, Qcutoff, minC, MTref_name=MTref_name)

    #the variant calling for each cell/sample
    if dataType in ['10x_scRNA-seq', 'BD-Rhapsody_scRNA-seq']:
        Cell_mutations = SNP_caller_10x(Cell_reads, Coverage_Cell, MT_ref, path_out, sampleID, Qcutoff=Qcutoff)
        New_Cell_mutations = MTmutations_filters(Cell_mutations, sampleID, ln=ln)
    if dataType == '10x_mtscATAC-seq':
        Cell_mutations = SNP_caller_10x(Cell_reads, Coverage_Cell, MT_ref, path_out, sampleID, Qcutoff=Qcutoff)
    if dataType in ['smart-seq2', 'bulk_RNA-seq']:
        Site_mutations = SNP_caller_RNA(Sfile, MT_ref, sampleID, path_out, Qcutoff=Qcutoff, MTref_name=MTref_name)
    if dataType in ['scATAC-seq', 'bulk_ATAC-seq']:    
        Site_mutations = SNP_caller_ATAC(Sfile, MT_ref, sampleID, path_out, Qcutoff=Qcutoff, MTref_name=MTref_name) 
    
    #Write results into txt
    if dataType in ['10x_scRNA-seq', 'BD-Rhapsody_scRNA-seq']:
        SNP_to_txt_10x(New_Cell_mutations, Read_Num, coverageF_stat, sampleID, path_out)
    elif dataType == '10x_mtscATAC-seq':
        SNP_to_txt_10x(Cell_mutations, Read_Num, coverageF_stat, sampleID, path_out)
    else:
        SNP_to_txt(Site_mutations, sampleID, path_out)


# In[ ]:

if __name__ == "__main__":
    t0=time.time()
    main()
    print ("Total time elapsed: %f" %(time.time()-t0))

