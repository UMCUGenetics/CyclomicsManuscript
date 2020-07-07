import sys, math, re, glob, os, random, gzip, warnings, pprint

from itertools import combinations_with_replacement
from statistics import mean
from collections import Counter, defaultdict
from subprocess import run, check_output

import numpy as np
import pandas as pd
from scipy import stats

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from matplotlib import rcParams

from Bio import pairwise2, Entrez, SeqIO
from Bio.SubsMat import MatrixInfo as matlist
from Bio import pairwise2



#####################################################
##################    BACKBONES    ##################
#####################################################
BB22 = 'GGGCATGCACAGATGTACACGTGACGCAACGANTGATGTTAGCTATTTGTTCAATGACATATNCTGGTATGATCAATACNAGATCTGATATTGATATNCTGATACTCATATATGTAGAATATCACATTATATTTATTATAATACATCGTCGAACATATACACAATGCATCTTATCTATACGTATCGGGATAGCGTTGGCATAGCACTGGATGGCATGACCCTCATTAGATGCTGCATGACATAGCCC'
BB24 = 'GGGCATGCACAGATGTACACGAATCCCGAAGANTGTTGTCCATTCATTGAATATGAGATCTCNATGGTATGATCAATATNCGGATGCGATATTGATANCTGATAAATCATATATGCATAATCTCACATTATATTTATTATAATAAATCATCGTAGATATACACAATGTGAATTGTATACAATGGATAGTATAACTATCCAATTTCTTTGAGCATTGGCCTTGGTGTAGATGCTGCATGACATAGCCC'
BB25 = 'GGGCATGCACAGATGTACACGAATCCGTGAGANTGAAGATCTTATTTGTGACATTCATCGATNCTGGATATGATCAATANCCATGCGATATTGATTANCTGATAAATCATATATGTAGAATATCACATTATATTAATTATAATAAATCGTCGTACATATACATCCACAATTAGCTATGTATACTATCTATAGAGATGGTGCATCATCGTACTCCACCATTCCCACTAGATGCTGCATGACATAGCCC'
BBCR = 'GGGCGGTATGTCATGCACACGAATCCCGAAGANTGTTGTCCATTCATTGAATATGAGATCTCNATGGTATGATCAATATNCGGATGCGATATTGATANCTGATAAATCATATATGCATAATCTCACATTATATTTATTATAATAAATCATCGTAGATATACACAATGTGAATTGTATACAATGGATAGTATAACTATCCAATTTCTTTGAGCATTGGCCTTGGTGTAGATTGCATGACATACCGCCC'
PJET = 'ATCTTTCTAGAAGATCTCCTACAATATTCTCAGCTGCCATGGAAAATCGATGTTCTTCTTTTATTCTCTCAAGATTTTCAGGCTGTATATTAAAACTTATATTAAGAACTATGCTAACCACCTCATCAGGAACCGTTGTAGGTGGCGTGGGTTTTCTTGGCAATCGACTCTCATGAAAACTACGAGCTAAATATTCAATATGTTCCTCTTGACCAACTTTATTCTGCATTTTTTTTGAACGAGGTTTAGAGCAAGCTTCAGGAAACTGAGACAGGAATTTTATTAAAAATTTAAATTTTGAAGAAAGTTCAGGGTTAATAGCATCCATTTTTTGCTTTGCAAGTTCCTCAGCATTCTTAACAAAAGACGTCTCTTTTGACATGTTTAAAGTTTAAACCTCCTGTGTGAAATTGTTATCCGCTCACAATTCCACACATTATACGAGCCGGAAGCATAAAGTGTAAAGCCTGGGGTGCCTAATGAGTGAGCTAACTCACATTAATTGCGTTGCGCTCACTGCCAATTGCTTTCCAGTCGGGAAACCTGTCGTGCCAGCTGCATTAATGAATCGGCCAACGCGCGGGGAGAGGCGGTTTGCGTATTGGGCGCTCTTCCGCTTCCTCGCTCACTGACTCGCTGCGCTCGGTCGTTCGGCTGCGGCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTTATCCACAGAATCAGGGGATAACGCAGGAAAGAACATGTGAGCAAAAGGCCAGCAAAAGGCCAGGAACCGTAAAAAGGCCGCGTTGCTGGCGTTTTTCCATAGGCTCCGCCCCCCTGACGAGCATCACAAAAATCGACGCTCAAGTCAGAGGTGGCGAAACCCGACAGGACTATAAAGATACCAGGCGTTTCCCCCTGGAAGCTCCCTCGTGCGCTCTCCTGTTCCGACCCTGCCGCTTACCGGATACCTGTCCGCCTTTCTCCCTTCGGGAAGCGTGGCGCTTTCTCATAGCTCACGCTGTAGGTATCTCAGTTCGGTGTAGGTCGTTCGCTCCAAGCTGGGCTGTGTGCACGAACCCCCCGTTCAGCCCGACCGCTGCGCCTTATCCGGTAACTATCGTCTTGAGTCCAACCCGGTAAGACACGACTTATCGCCACTGGCAGCAGCCACTGGTAACAGGATTAGCAGAGCGAGGTATGTAGGCGGTGCTACAGAGTTCTTGAAGTGGTGGCCTAACTACGGCTACACTAGAAGGACAGTATTTGGTATCTGCGCTCTGCTGAAGCCAGTTACCTTCGGAAAAAGAGTTGGTAGCTCTTGATCCGGCAAACAAACCACCGCTGGTAGCGGTGGTTTTTTTGTTTGCAAGCAGCAGATTACGCGCAGAAAAAAAGGATCTCAAGAAGATCCTTTGATCTTTTCTACGGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGGATTTTGGTCATGAGATTATCAAAAAGGATCTTCACCTAGATCCTTTTAAATTAAAAATGAAGTTTTAAATCAATCTAAAGTATATATGAGTAAACTTGGTCTGACAGTTACCAATGCTTAATCAGTGAGGCACCTATCTCAGCGATCTGTCTATTTCGTTCATCCATAGTTGCCTGACTCCCCGTCGTGTAGATAACTACGATACGGGAGGGCTTACCATCTGGCCCCAGTGCTGCAATGATACCGCGAGACCCACGCTCACCGGCTCCAGATTTATCAGCAATAAACCAGCCAGCCGGAAGGGCCGAGCGCAGAAGTGGTCCTGCAACTTTATCCGCCTCCATCCAGTCTATTAATTGTTGCCGGGAAGCTAGAGTAAGTAGTTCGCCAGTTAATAGTTTGCGCAACGTTGTTGCCATTGCTACAGGCATCGTGGTGTCACGCTCGTCGTTTGGTATGGCTTCATTCAGCTCCGGTTCCCAACGATCAAGGCGAGTTACATGATCCCCCATGTTGTGCAAAAAAGCGGTTAGCTCCTTCGGTCCTCCGATCGTTGTCAGAAGTAAGTTGGCCGCAGTGTTATCACTCATGGTTATGGCAGCACTGCATAATTCTCTTACTGTCATGCCATCCGTAAGATGCTTTTCTGTGACTGGTGAGTACTCAACCAAGTCATTCTGAGAATAGTGTATGCGGCGACCGAGTTGCTCTTGCCCGGCGTCAATACGGGATAATACCGCGCCACATAGCAGAACTTTAAAAGTGCTCATCATTGGAAAACGTTCTTCGGGGCGAAAACTCTCAAGGATCTTACCGCTGTTGAGATCCAGTTCGATGTAACCCACTCGTGCACCCAACTGATCTTCAGCATCTTTTACTTTCACCAGCGTTTCTGGGTGAGCAAAAACAGGAAGGCAAAATGCCGCAAAAAAGGGAATAAGGGCGACACGGAAATGTTGAATACTCATACTCTTCCTTTTTCAATATTATTGAAGCATTTATCAGGGTTATTGTCTCATGAGCGGATACATATTTGAATGTATTTAGAAAAATAAACAAATAGGGGTTCCGCGCACATTTCCCCGAAAAGTGCCACCTGACGTCTAAGAAACCATTATTATCATGACATTAACCTATAAAAATAGGCGTATCACGAGGCCGCCCCTGCAGCCGAATTATATTATTTTTGCCAAATAATTTTTAACAAAAGCTCTGAAGTCTTCTTCATTTAAATTCTTAGATGATACTTCATCTGGAAAATTGTCCCAATTAGTAGCATCACGCTGTGAGTAAGTTCTAAACCATTTTTTTATTGTTGTATTATCTCTAATCTTACTACTCGATGAGTTTTCGGTATTATCTCTATTTTTAACTTGGAGCAGGTTCCATTCATTGTTTTTTTCATCATAGTGAATAAAATCAACTGCTTTAACACTTGTGCCTGAACACCATATCCATCCGGCGTAATACGACTCACTATAGGGAGAGCGGCCGCCAGATCTTCCGGATGGCTCGAGTTTTTCAGCAAGAT'
#####################################################
#####################################################
#####################################################


#####################################################
###################    INSERTS    ###################
#####################################################
A12_WT = 'CTTGCTTACCTCGCTTAGTGCTCCCTGGGGGCAGCTCGTGGTGAGGCTCCCCTTTCTTGCGGAGATTCTCTTCCTCTGTGCGCCGGTCTCTCCCAGGACAGGCACAAACACGCACCTCAAAGCTGTTCCGTCCCAGTAGAT' #Synthetic

A12_MU = 'CTTGCTTACCTCGCTTAGTGCTCCCTGGGGGCAGCTCGTGGTGAGGCTCCCCTTTCTTGCGGAGATTCTCTTCCTCTGTGCGCCAGTCTCTCCCAGGACAGGCACAAACATACACCTCAAAGCTGTTCCGTCCCAGTAGAT' #Synthetic

S1_WT = 'GGTGGCAAGTGGCTCCTGACCTGGAGTCTTCCAGTGTGATGATGGTGAGGATGGGCCTCCGGTTCATGCCGCCCATGCAGGAACTGTTACACATGTAGTTGTAGTGGATGGTGGTACAGTCAGAGCCAACCTAGGAGATAACACAGGCCCAAGATGA' #Synthetic

S1_MU = 'GGTGGCAAGTGGCTCCTGACCTGGAGTCTTCCAGTGTGAAGATGGTGAGGATGGGCCTCCGGTTCATGCCGCCCATGCAGGAACTGTTACACATGTAGTTGTCGTGGATGGTGGTACAGTCAGAGCCAACCTAGGAGATAACACAGGCCCAAGATGA' #Synthetic

TP53_5 = 'TCCTTCCCAGCCTGGGCATCCTTGAGTTCCAAGGCCTCATTCAGCTCTCGGAACATCTCGAAGCGCTCACGCCCACGGATCTGCAGCAACAGAGGAGGGGGAGAAGTAAGTATATACA'

TP53_9 = 'ACTGGAAACTTTCCACTTGATAAGAGGTCCCAAGACTTAGTACCTGAAGGGTGAAATATTCTCCATCCAGTGGTTTCTTCTTTGGCTGGGGAGAGGAGCTGGTGTTGTTGGGCAGTGCTAGGAAAGAGGCAAGGAAAGGTGATAAA'

TP53_12 = 'CTTGCTTACCTCGCTTAGTGCTCCCTGGGGGCAGCTCGTGGTGAGGCTCCCCTTTCTTGCGGAGATTCTCTTCCTCTGTGCGCCGGTCTCTCCCAGGACAGGCACAAACACGCACCTCAAAGCTGTTCCGTCCCAGTAGAT'

TP53_19 = 'TGTCGTCTCTCCAGCCCCAGCTGCTCACCATCGCTATCTGAGCAGCGCTCATGGTGGGGGCAGCGCCTCACAACCTCCGTCATGTGCTGTGACTGCTTGTAGATGGCCAT'

TP53_20 = 'CTCACAACCTCCGTCATGTGCTGTGACTGCTTGTAGATGGCCATGGCGCGGACGCGGGTGCCGGGCGGGGGTGTGGAATCAACCCACAGCTGCACAGGGCAGGTCTTGGCCAGTTGGCAAAACATCTTGTTGAGGGC'

TP53_21 = 'GGGGTGTGGAATCAACCCACAGCTGCACAGGGCAGGTCTTGGCCAGTTGGCAAAACATCTTGTTGAGGGCAGGGGAGTACTGTAGGAAGAGGAAGGAGACAGAGTTGAAAGTCAGGG'
#####################################################
#####################################################
#####################################################


#####################################################
###############    Other variables    ###############
#####################################################
pp = pprint.PrettyPrinter(indent=4,depth=6)

palette = {'purple':"#9b59b6",
           'blue':"#3498db",
           'gray':"#95a5a6",
           'red':"#e74c3c",
           'black':"#34495e",
           'green':"#2ecc71"}

style = 'seaborn-whitegrid' #matplotlib style
#####################################################
#####################################################
#####################################################



#####################################################
##################    Functions    ##################
#####################################################
def fp_percent(row):
    try:
        ref = row['REF']
        return percent_of(row['COV'], row['COV']-row[ref])
    except:
        return np.nan


def fp_percent_no_del(row):
    try:
        ref = row['REF']
        cov = row['COV'] - row['DEL']
        return percent_of(cov, cov-row[ref])
    except:
        return np.nan

def fp_cosmic(row):
    try:
        assert not row['REF'] in row['COSMIC']
    except AssertionError as e:
        #print(row['REF'], row['COSMIC'])
        return np.nan
    muts = row['COSMIC']
    fps = 0
    for m in muts:
        try:
            fps += row[m]
        except Exception as e:
            #print(e)
            pass
    try:
        if 'DEL' in muts:
            cov = row['COV']
        else:
            #correct for DELs
            cov = row['COV'] - row['DEL']
        return percent_of(cov, fps)
    except ZeroDivisionError:
        pass
    return np.nan

    
def get_seq(row):
    seq = sequence_from_coordinates(
        chromosome = row['REF'],
        strand=1,
        start=row['POS']+1,
        end=row['POS']+1,
        ref_genome=37
    )
    return seq


def get_ref(row):
    return row['GRCh37 coord. (chr:base)'].split(':')[0]


def get_pos(row):
    return int(row['GRCh37 coord. (chr:base)'].split(':')[1])


def get_mut_alleles(row):
    '''Parse SNP and single base deletions only'''
    if 'Substitution' in row['Type']:
        m = row['c.Mutation'][-1]
        m = complement(m) #TP53 is on the minus strand
        return m
    elif 'Deletion' in row['Type']:
        return 'DEL'
    else:
        return np.nan
    
    
def get_wt_alleles(row):
    '''Parse SNP and single base deletions only'''
    if 'Substitution' in row['Type']:
        m = row['c.Mutation'][-3]
        m = complement(m) #TP53 is on the minus strand
        return m
    elif 'Deletion' in row['Type']:
        try:
            m = row['c.Mutation'].split('del')[-1]
            m = reverse(complement(m))
            return m
        except:
            #print(row['c.Mutation'])
            pass
    else:
        return np.nan

    
def match_references(row):
    cy_ref = row['REF']
    cy_pos = row['POS']
    cosmic_ref = set(db[db.POS == cy_pos]['g.Ref_Allele'])
    if cy_ref in cosmic_ref:
        return True
    return f'{cy_pos} {cy_ref} {cosmic_ref}'


def delta_percent(before, after, debug=False):
    np.seterr(divide='ignore', invalid='ignore')
    
    before = float(before)
    after  = float(after)
    
    try:
        x = ((after-before)/abs(before))*100
        if x == 0.0:
            return 0.000000001 #avoid -inf
        else:
            return x
    except Exception as e:
        if debug: print('Exception raised by delta_percent():',e)
        return 0.000000001 #avoid -inf


def vcf_to_df(file):
    skiprows = False
    with open(file, 'r') as f:
        for i,line in enumerate(f):
            if line.startswith('#CHROM'):
                skiprows = i
                break
    return pd.read_csv(file, sep='\t', skiprows=skiprows) 


def compute_Qscore(row, ref, with_deletions=True):
    '''
    Compute Phred Quality Score of CyclomicsSeq reads grouped by number of repeats.
    Works with df generated from cont files (e.g. consensus_8_full_consensus_sambamba_ouput_BB200_4.txt)
    Args:
        row: is a Pandas.Series as obtained by df.iterrows()
        ref: a reference sequence
        with_deletions: Include/Exclude DEL calls from the stats
    '''
    base = ref[row['POS']]
    if with_deletions:
        errors = sum([row[x] for x in ['A','C','G','T','DEL'] if x != base]) or 0.1 #to avoid 0/something => inf
        coverage = row['COV']
        
    else:
        errors = sum([row[x] for x in ['A','C','G','T'] if x != base]) or 0.1
        coverage = row['COV']-row['DEL']

    Q = -10*np.log10(errors/coverage)
    return Q, coverage, errors


def quick_Qscore(coverage, errors):
    if coverage == 0:
        return np.nan
    if errors == 0:
        errors = 0.1
    
    Q = -10*np.log10(errors/coverage)
    return Q


def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    https://stackoverflow.com/questions/5967500/how-to-correctly-sort-a-string-with-a-number-inside
    '''
    return [atoi(c) for c in re.split(r'(\d+)', text)]


def jitter(n, mu=0, sigma=0.1):
    '''Return a jittered version of n'''
    return n + np.random.normal(mu, sigma, 1)

    
def plot_calls(
    df, #dataframe
    indels = False,
    xs=np.array([0,100]), #starting boundaries of scratter line
    dx=np.array([10,-10]), #borders between positions
    xt=50, #xticks starting position
    line_width=100, #lenght of scratter line
    savefig = False
):
    '''...'''

    _xticks = [xt]
    plt.style.use('ggplot')
    plt.figure(figsize=(16,9))
    rcParams.update({'font.size': 22})
    plt.title(f"{df.iloc[0]['SAMPLE']}\n")
    palette = {'purple':"#9b59b6",
               'blue':"#3498db",
               'gray':"#95a5a6",
               'red':"#e74c3c",
               'black':"#34495e",
               'green':"#2ecc71"}

    for i, row in df.iterrows():
        #print(_xticks)
        #print(xs)
        A = plt.scatter(
            [random.randint(*(xs+dx)) for n in range(row['A'])],
            [jitter(0, 0, 1) for n in range(row['A'])],
            c=palette['red'], alpha=0.3, s=30 
        )
        C = plt.scatter(
            [random.randint(*(xs+dx)) for n in range(row['C'])],
            [jitter(5, 0, 1) for n in range(row['C'])],
            c=palette['purple'], alpha=0.3, s=30
        )
        G = plt.scatter(
            [random.randint(*(xs+dx)) for n in range(row['G'])],
            [jitter(10, 0, 1) for n in range(row['G'])],
            c=palette['blue'], alpha=0.3, s=30
        )
        T = plt.scatter(
            [random.randint(*(xs+dx)) for n in range(row['T'])],
            [jitter(15, 0, 1) for n in range(row['T'])],
            c=palette['green'], alpha=0.3, s=30
        )
        if indels:
            DEL = plt.scatter(
                [random.randint(*(xs+dx)) for n in range(row['DEL'])],
                [jitter(20, 0, 1) for n in range(row['DEL'])],
                c=palette['black'], alpha=0.3, s=30
            )
            
        xs[0] += line_width
        xs[1] += line_width
        _xticks.append(_xticks[-1]+line_width)
        
    plt.xticks(_xticks, df['POS'], rotation=90)
    plt.ylabel('base called')
    plt.xlabel('reference position')
    
    if indels:
        leg = plt.legend([A,C,G,T, DEL],['A','C','G','T','DEL'])
        plt.yticks([0, 5, 10, 15, 20], ['A','C','G','T','DEL'])
        plt.ylim(-5, 25)
            
    else:
        leg = plt.legend([A,C,G,T],['A','C','G','T'])
        plt.yticks([0, 5, 10, 15], ['A','C','G','T'])
        plt.ylim(-5, 20)
    
    for l in leg.legendHandles:
        l.set_alpha(1)

    if savefig:
        plt.savefig(savefig)
    plt.show()
    plt.close('all')

def simple_consensus(aligned_sequences_file, ignore_indels=True):
    '''file => string
    Return the consensus of a series of fasta sequences aligned with muscle.
    '''
    # Generate consensus from Muscle alignment
    sequences = []
    seq = False
    ref_seq = ''
    with open(aligned_sequences_file,'r') as f:
        for line in f:
            if line.startswith('\n'):
                continue
            if line.startswith('>'):
                if seq:
                    sequences.append(seq)
                try:
                    int(line[1])
                except ValueError:
                    #it's the ref sequence
                    ref_seq += seq
                    
                seq = ''
            else:
                seq += line.strip()

    #check if all sequenced have the same length
    for seq in sequences:
        assert len(seq) == len(sequences[0])
    
    #compute consensus by majority vote
    consensus = ''
    try:
        for i in range(len(sequences[0])):
            char_count = Counter()
            for seq in sequences:
                char_count.update(seq[i])
                top_candidate = char_count.most_common()[0][0]
                #Ignore deletions
                try:
                    if top_candidate == '-' and ignore_indels:
                        top_candidate = char_count.most_common()[1][0]
                except IndexError:
                    top_candidate = ''
            consensus += top_candidate
    except IndexError:
        print('Error in',aligned_sequences_file)
        pass

    if ignore_indels:
        for i,char in enumerate(ref_seq):
            if char == '-':
                consensus = consensus[:i]+consensus[i+1:]
        
    return consensus#.replace('-','')


def array_percent(array):
    '''
    Compute the percentage weight of each item in the array.

    :param list array: a list/array of numbers

    :return: list of weigths from 0 to 100 (percent)
    :rtype: list
    
    Example:
        >>> array_percent([30,30,15])
        [40.0, 40.0, 20.0]
    '''
    tot = sum(array)
    r = []
    for i in array:
        r.append(percent_of(tot,i))
    return r
            
    
def parse_fasta(fasta_file):
    '''file_path => dict
    Return a dict of id:sequences.
    '''
    d = {}
    _id = False
    seq = ''
    with open(fasta_file,'r') as f:
        for line in f:
            if line.startswith('\n'):
                continue
            if line.startswith('>'):
                if not _id:
                    _id = line[1:].strip()
                elif _id and seq:
                    d.update({_id:seq})
                    _id = line[1:].strip()
                    seq = ''
            else:
                seq += line.strip()
        d.update({_id:seq})
    return d


def quick_align(reference, sample, matrix=matlist.blosum62, gap_open=-10, gap_extend=-0.5):
    '''
    Return a binary score matrix for a pairwise alignment.
    '''

    alns = pairwise2.align.globalds(reference, sample, matrix, gap_open, gap_extend)

    top_aln = alns[0]
    aln_reference, aln_sample, score, begin, end = top_aln

    score = []
    for i, base in enumerate(aln_reference):
        if aln_sample[i] == base:
            score.append(1)
        else:
            score.append(0)

    return score


def reverse(sequence):
    r = ''
    for i in range(len(sequence),0,-1):
        r += sequence[i-1]
    return r


def complement(sequence):
    d = {'A':'T','a':'t',
         'T':'A','t':'a',
         'C':'G','c':'g',
         'G':'C','g':'c'}
    r = ''
    for b in sequence.upper():
        r += d[b]
    return r

def dice_coefficient(sequence_a, sequence_b):
    '''(str, str) => float
    Return the dice cofficient of two sequences.
    '''
    a = sequence_a
    b = sequence_b
    if not len(a) or not len(b): return 0.0
    # quick case for true duplicates
    if a == b: return 1.0
    # if a != b, and a or b are single chars, then they can't possibly match
    if len(a) == 1 or len(b) == 1: return 0.0
    
    # list comprehension, preferred over list.append() '''
    a_bigram_list = [a[i:i+2] for i in range(len(a)-1)]
    b_bigram_list = [b[i:i+2] for i in range(len(b)-1)]
    
    a_bigram_list.sort()
    b_bigram_list.sort()
    
    # assignments to save function calls
    lena = len(a_bigram_list)
    lenb = len(b_bigram_list)
    # initialize match counters
    matches = i = j = 0
    while (i < lena and j < lenb):
        if a_bigram_list[i] == b_bigram_list[j]:
            matches += 2
            i += 1
            j += 1
        elif a_bigram_list[i] < b_bigram_list[j]:
            i += 1
        else:
            j += 1
    
    score = float(matches)/float(lena + lenb)
    return score


def extract_U(sequence, structure):
    segments = structure.split(',')
    m = []
    for s in segments:
        l,t = s.split(':')[:2]
        m.append((t,int(l)))
    
    unmapped_sequences = []
    i = 0
    for t,l in m:
        if t == 'U':
            unmapped_sequences.append(sequence[i:i+l])
        i += l
        
    return unmapped_sequences


def list_of_files(path, extension, recursive=False):
    '''
    Return a list of filepath for each file into path with the target extension.
    If recursive, it will loop over subfolders as well.
    '''
    if not recursive:
        for file_path in glob.iglob(path + '/*.' + extension):
            yield file_path
    else:
        for root, dirs, files in os.walk(path):
            for file_path in glob.iglob(root + '/*.' + extension):
                yield file_path

                
                
                
def split_overlap(iterable,size,overlap):
    '''(list,int,int) => [[...],[...],...]
    Split an iterable into chunks of a specific size and overlap.
    Works also on strings! 

    Examples:
        >>> split_overlap(iterable=list(range(10)),size=3,overlap=2)
        [[0, 1, 2, 3], [2, 3, 4, 5], [4, 5, 6, 7], [6, 7, 8, 9]]

        >>> split_overlap(iterable=range(10),size=3,overlap=2)
        [range(0, 3), range(1, 4), range(2, 5), range(3, 6), range(4, 7), range(5, 8), range(6, 9), range(7, 10)]
    '''
    if size < 1 or overlap < 0:
        raise ValueError('"size" must be an integer with >= 1 while "overlap" must be >= 0')
    result = []
    while True:
        if len(iterable) <= size:
            result.append(iterable)
            return result
        else:
            result.append(iterable[:size])
            iterable = iterable[size-overlap:]
            
            

def is_concatemer(olist, n=2, limit=0.75, verbose=False):
    '''
    Return True if a list is made of n elements that regularly alternates (ordinated list).
    Else, return False.
    '''
    #implementation for n==2
    r = 0
    if n == 2:
        try:
            e0,e1 = set(olist)
        except ValueError as e:
            if verbose:
                print(e)
            return False
        blocks = split_overlap(olist,n,1)
        for i in blocks:
            if i == [e0,e1] or i == [e1,e0]:
                r += 1
                
        correct = r/len(blocks)
    if verbose:
        print(f'correct: {r} out of {len(blocks)} = {correct}')
        
    return True if correct >= limit else False


def entropy(sequence):
    '''(string, bool) => float
    Return the Shannon Entropy of a string.
    Calculated as the minimum average number of bits per symbol
    required for encoding the string.
    The theoretical limit for data compression:
    Shannon Entropy of the string * string length.
    '''
    
    alphabet = set(sequence) # list of symbols in the string

    # calculate the frequency of each symbol in the string
    frequencies = []
    for symbol in alphabet:
        frequencies.append(sequence.count(symbol) / len(sequence))

    # Shannon entropy
    ent = 0.0
    for freq in frequencies:
        ent -= freq * math.log(freq, 2)
        
    return ent



def get_reps(sequence):
    '''str => [(str, str, float), ...]
    Return a list of tuple (generator), each containing:
    [0] the sequence of a concatemer
    [1] its building block
    [2] the number of repetition of the building block in the concatemer sequence
    
    Since the regex is computationally expensive,
    collapsed strings are the preferred inputs. ##See collapse(structure)
    
    Example:
            >>> s = 'U0U1U0U0U0U1U0U1U0U0U1U0U1U0U1U0U1U0U1U0U1U0U1U0U1U0U1U0U1U0U1U0U'
            >>> list(get_reps(s))
            [('U0U0U0', 'U0', 3.0),
             ('U1U0U1U0', 'U1U0', 2.0),
             ('U0U1U0U1U0U1U0U1U0U1U0U1U0U1U0U1U0U1U0U1U0U1', 'U0U1', 11.0)]
    
    '''
    return ((s, b,len(s)/len(b)) for s, b in re.findall(r'((\w+?)\2+)', sequence))



def get_reps2(sequence):
    """Find the most repetitive sequence in a string.

    :param str sequence: string for search
    :param int rep_min_len: minimal length of repetitive substring
    :return the most repetitive substring or None
    """
    greedy, non_greedy = re.compile(r'((\w+)\2+)'), re.compile(r'((\w+?)\2+)')

    all_rep_seach = lambda regex: \
        (regex.search(sequence[shift:]) for shift in range(len(sequence)))

    return (res.groups() for res in chain(all_rep_seach(greedy), all_rep_seach(non_greedy)) if res)




def print_sbar(n,m,s='|#.|',size=30,message=''):
    '''(int,int,string,int) => None
    Print a progress bar using the simbols in 's'.
    Example:
        range_limit = 1000
        for n in range(range_limit):
            print_sbar(n+1,m=range_limit)
            time.sleep(0.1)
    '''
    #adjust to bar size
    if m != size:
        n =(n*size)/m
        m = size
    #calculate ticks
    _a = int(n)*s[1]+(int(m)-int(n))*s[2]
    _b = round(n/(int(m))*100,1)
    #adjust overflow
    if _b >= 100:
        _b = 100.0
    #to stdout    
    sys.stdout.write(f'\r{message}{s[0]}{_a}{s[3]} {_b}%     ')
    sys.stdout.flush()        
        

def percent_of(total, fraction):
    '''(int_or_float,int_or_float) => float
    Return the percentage of 'fraction' in 'total'.
    
    Examples:
        >>> percent_of(150, 75)
        50.0
        
        >>> percent_of(30, 90)
        300.0
    '''
    return (100*fraction)/total


def add_Uid_to_structure(structure, bb_len, ins_len):
    '''
    Add a Uid to each 'U' segment of a structure.
    Return a new structure.

    :param str structure: a structure
    :param int bb_len: len of BB in bp
    :param int ins_len: len of I in bp

    :return: a new structure
    :rtype: str

    '''
    #The ID is given based on the expected correct length using solve_U().
    #First we need to get a map of possible U regions,
    #if we know the length of the insert and the one of the backbones 
    #we get the map by using get_Ubins()
    u_bins = gen_Ubins(bb_len, ins_len)

    #Now we can add this info to the U segments
    new_s = ''
    for segment in structure.split(','):
        if segment[-1] == 'U':
            u_len = int(segment.split(':')[0])
            _ = list(solve_U(u_len, u_bins).keys())
            assert len(_) == 1
            Uid = _[0]
            segment += f':{Uid}'
        new_s += segment + ','

    return new_s[:-1] #strip last ','
            

def collapse(structure, unmapped_char=None, Uid=False, remove_dust=0, trim_edges=False):
    '''(str, None/str, Bool/int, Bool => (str, dict)
    Maps each block to a number.
    Return the collapsed string plus the a his mapping dict.
    
    In the mapping dict the value:
        'B' indicates a backbone
        'I' indicates an insert
        'F'/'R' indicates the orientation
        All other info are discarded.
        
    Each unique combination of mapping values will be encoded into a an integer
    that will used as a key in the mapping dict.
    The same integer is then used to encode the represented value into the returned string.
    
    ##OLD
    'U' stays always the same,
    so it is directly encoded into the string but id doesn't need to be mapped.
    
    ##NEW
    'U' can have an id as well, in that case, pass True to the Uid param.
    
    Example:
        s     = '28:U,195:BB:R:BB200_2:1:201:0,1:U,153:I:F:17:7577479:158:0,
        242:BB:R:BB200_2:1:241:0,5:U,170:I:F:17:7577484:153:0,239:BB:R:BB200_2:6:238:0,
        164:I:F:17:7577471:155:0,14:U,253:BB:R:BB200_2:0:243:0,182:I:F:17:7577478:159:0,
        191:BB:R:BB200_2:50:192:0,1407:U,217:BB:R:BB200_2:2:228:0,6:U,
        121:I:F:17:7577495:140:0,3:U,200:BB:R:BB200_2:30:208:0,37:U,138:BB:F:BB200_2:0:151:0,
        192:U,141:BB:F:BB200_2:15:142:0,142:U,257:BB:F:BB200_2:0:244:0,
        161:I:R:17:7577478:157:0,47:U,162:BB:F:BB200_2:67:177:0,22:U,131:I:R:17:7577481:138:0,
        3:U,230:BB:F:BB200_2:1:231:0,13:U,152:I:R:17:7577483:154:0,4:U,
        241:BB:F:BB200_2:0:237:0,6:U,152:I:R:17:7577481:156:0,74:U,
        182:BB:R:BB200_2:1:183:0,709:U,155:I:F:17:7577478:159:0,
        10:U,224:BB:R:BB200_2:0:234:0,25:U,109:I:F:17:7577506:131:0,229:BB:R:BB200_2:1:242:0,
        39:U,91:BB:F:BB200_2:30:90:0,559:U'
        
        new_s = add_Uid_to_structure(s, bb_len=245, ins_len=150)
        
        print(collapse(s))
        ('U2U32U323U232U2U3U2U0U0U01U0U1U0U1U0U1U2U3U2U32U0U',
        {0: 'B0F', 1: 'I0R', 2: 'B0R', 3: 'I0F'})
        
        print(collapse(new_s, Uid=True))
        ('32372372732726237323010105303530353035328732372304',
        {0: 'B0F', 1: 'U150', 2: 'B0R', 3: 'U0', 4: 'U545', 5: 'I0R', 6: 'U980', 7: 'I0F', 8: 'U695'})
        
    '''
    string = ''
    for s in structure.split(','):
            f = s.split(':') #fields
            
            #remove dust, i.e. very small mapped or unmapped regions
            size = int(f[0])
            if size <= remove_dust:
                continue
            
            s_type = f[1]
            string += s_type
            
            if (not Uid and s_type != 'U') or Uid: #'U' has only two fields
                orientation    = f[2] if s_type != 'U' else ''
                id_number      = int(f[-1])
                string += str(id_number) + orientation
            
                
    c_str = string.replace('BB','B')
    #c_str = string.replace('U','*')
    #print(c_str)
    r = [x for x in set(split_overlap(c_str,3,2)) if x[-1] in 'FR']
    
    if Uid:
        u = []
        for x in set(split_overlap(c_str,4,3)): #supposing and ID of 3 chars i.e. U145
            if x[0] == 'U':
                if x[1] == '0': #U0 is still an option
                                #check add_Uid_to_structure() for more info
                    u.append(x[:2])
                else:
                    u.append(x)
        r += u
    
    r = set(r)
    r = list(zip(range(len(r)),r))
    
    for n, seq in r:
        c_str = c_str.replace(seq,str(n))
        
    if unmapped_char:
        c_str = c_str.replace('U', unmapped_char)
    
    
    #remove edges if they are unmapped
    if trim_edges:
        if unmapped_char:
            pass
        else:
            unmapped_char = 'U'
        if c_str[0] == unmapped_char:
            c_str = c_str[1:]
        if c_str[-1] == unmapped_char:
            c_str = c_str[:-1]

            
    return c_str, dict(r)


def plot_structure(structure, bb_colors, insert,
                   expected_chr,
                   vmin, vmax,
                   title='untitled',
                   show=True, savefig=False):
    '''
    Generate read-plots from a structure string.
    'vmin' and 'vmax' specify the coordinates of the region of interest
    (for example TP53) on the 'expected_chr'.
    Each read is rendered using 3 different color maps
    in order to highlight any chimerism. (different inserts = different color)
    To save the plot on a file, pass the file path (../file.png or .pdf) to 'savefig')
    '''
    
    ###init
    ##Colormap normalization based on TP53 coordinates 17:7565097-7590856 (GRCh37)
    #https://matplotlib.org/users/colormaps.html
    cmap_0 = cm.gist_rainbow
    cmap_1 = cm.prism #sharper difference between close numbers but not unique
    cmap_2 = cm.tab20
    cmap_chr = cm.hsv #for chromosomes
    #normalize coordinates
    norm = Normalize(vmin, vmax)
    norm_chr = Normalize(1, 24) #for chromosomes

    
    segments = structure.split(',')
    xes = [] #keep track of the coordinates of each segment
    y  = 0 #y value for the middle plot
    dy = 0 #y distance for arrow coordinate (it should be 0 for orizontal arrows)
    read_len = 0
    plt.style.use('ggplot')
    plt.figure(figsize=(20,4), dpi=120) #figuresize must be inside the for-loop

    for n, s in enumerate(segments):

        f = s.split(':') #segment fields
        start_pos_read = read_len
        segment_len = int(f[0])
        orientation = False
        s_type = f[1]
        
        if s_type != 'U': #'U' has only two fields
            orientation    = f[2]
            chromosome     = f[3] #in case of BB, it is the name of the BB (BB200_2)
            start_pos_ref  = int(f[4])
            map_len        = int(f[5])
            id_number      = int(f[6])

            ##Plot arrow
            #BB
            bb_color = False
            if s_type == 'BB':
                if chromosome in bb_colors:
                    bb_color = bb_colors[chromosome]
                else:
                    bb_color = 'black' #wrong BB

                #BB will be the same color on all the plots
                color_0 = bb_color
                color_1 = bb_color
                color_2 = bb_color

            #I
            elif s_type == 'I':
                is_other = False
                if chromosome == expected_chr:
                    if vmin <= start_pos_ref <= vmax:
                        #assign color based on coordinates
                        color_0 = cmap_0(norm(start_pos_ref))#'blue' #map to row[1]['start_position_ref']
                        color_1 = cmap_1(norm(start_pos_ref))
                        color_2 = cmap_2(norm(start_pos_ref))
                    else:
                        is_other = True
                else:
                    is_other = True

            #if I is not the correct one assign a single color for each chromosome       
            if s_type == 'I' and is_other:
                if chromosome.upper() == 'X':
                    chromosome = 23
                elif chromosome.upper() == 'Y':
                    chromosome = 24
                else:
                    try:
                        chromosome = int(chromosome)
                    except ValueError:
                        pass
                
                if type(chromosome) is int:
                    color_0 = cmap_chr(norm_chr(chromosome))
                    color_1 = cmap_chr(norm_chr(chromosome))
                    color_2 = cmap_chr(norm_chr(chromosome))
                else:
                    color_0 = 'black'
                    color_1 = 'black'
                    color_2 = 'black'
        else:
            #if s_type == 'U':
            color_0 = 'white'
            color_1 = 'white'
            color_2 = 'white'


        #get plotting coordinates for each fragment
        x  = start_pos_read
        xes.append(x) #for scale setting
        dx = segment_len

        #Update read_len
        read_len += segment_len

        #invert coordinates in case of R orientation
        if orientation == 'R':
            x += dx
            dx *= -1
        
        ##Plot!!
        #1
        plt.arrow(x, y+0.05, dx, dy,
                  shape='full',
                  lw=3,
                  length_includes_head=True,
                  head_width=.01 if s_type != 'U' else 0,
                  color=color_0,
                  alpha=1)
        #0
        plt.arrow(x, y, dx, dy,
                  shape='full',
                  lw=3,
                  length_includes_head=True,
                  head_width=.01 if s_type != 'U' else 0,
                  color=color_1,
                  alpha=1)
        #-1
        plt.arrow(x, y-0.05, dx, dy,
                  shape='full',
                  lw=3,
                  length_includes_head=True,
                  head_width=.01 if s_type != 'U' else 0,
                  color=color_2,
                  alpha=1)
               
    plt.title(title, fontsize=12)
    plt.ylim(-0.1, 0.1)
    plt.xlim(-int(mean(xes)/10), read_len+int(mean(xes)/10))

    if savefig:
        plt.savefig(savefig)

    if show:
        plt.show()
        
    plt.close('all') #to avoid RuntimeWarning: More than 20 figures have been opened.
    
    
def gen_Ubins(bb_len, ins_len, max_combinations_len=4):
    '''(int, int, int) => dict
    Return a dict of all of the possible lenghts (up to max_combinations_len)
    associated with unmapped occurrences of BB and I.
    '''
    d_len = {'B'    :bb_len,
             'I'    :ins_len}
    u_bins = {v:k for k,v in d_len.items()}
    for n in range(2, max_combinations_len+1):
        for i in combinations_with_replacement(list(d_len.keys()), n):
            v = sum([d_len[k] for k in i])
            assert v not in u_bins
            u_bins[v] = i
    u_bins[0] = 'E' #it's probably just a little mapping error]
    
    for k in u_bins:
        if k == 0:
            pass
        else:
            assert 100 <= k <= 999 #the unmapped region should be smaller than 999bp,
                                   #since the Uid can be 0 or a 3digits int
    return u_bins


def solve_U(u_len, u_bins):
    '''(int, dict) => dict
    Return the expected correct lenght of an 'U' segment
    and it's most lickley composition.
    '''
    d = {}
    for n in u_bins.keys():
        d[n] = abs(n-u_len)

    m = min(d.values())
    return {k:Counter(v) for k,v in u_bins.items() if d[k] == m}


def opposite(x):
    assert x in 'FR'
    return 'R' if x == 'F' else 'F'


def check_dna_shift(c_str, c_blocks):
    '''
    Return the type of template shift in the read.

    :param list c_str: a collapsed string (see cyclomics.collapse())
    :param dict c_blocks: a block table (see cyclomics.collapse())

    :return: one of **['None', '??', 'Unexpected', 'Regular']**
    :rtype: str
    '''

    combined = {''.join(sorted(k[1])) for k in c_blocks if 'U' not in k[1]}
    compare = {}
    for k in combined:
        compare[k] = {}
        for _id in k:
            try:
                v = c_str[-1][int(_id)]
                compare[k][v[:-1]] = v[-1]
            except ValueError:
                pass
    try:
        k1,k2 = compare.keys()
    except ValueError:
        return '??' #todo: detect complex shifts.

    passed = []
    for k,v in compare[k1].items():
        if k in compare[k2] and compare[k2][k] == opposite(v):
            passed.append(k1)
    
    if not passed:
        return 'None'
    if len(passed) == 1:
        return 'Unexpected'
    return 'Regular'


def segments_len(structure):
    d = {'BB':[],'I':[],'U':[]}
    segments = structure.split(',')
    for s in segments:
        v,k = s.split(':')[:2]
        d[k].append(int(v))
    for k,v in d.items():
        d[k] = stats.trim_mean(v, 0.1) #10% trimmed mean
    return d


#Downaload sequence from ensembl
def sequence_from_coordinates(chromosome, strand, start, end, ref_genome=37):
    '''
    Download the nucleotide sequence from the gene_name.
    '''
    Entrez.email = "a.marcozzi@umcutrecht.nl" # Always tell NCBI who you are
    
    if int(ref_genome) == 37:
        #GRCh37 from http://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.25/#/def_asm_Primary_Assembly
        NCBI_IDS = {'1':'NC_000001.10','2':'NC_000002.11','3':'NC_000003.11','4':'NC_000004.11',
                    '5':'NC_000005.9','6':'NC_000006.11','7':'NC_000007.13','8':'NC_000008.10',
                    '9':'NC_000009.11','10':'NC_000010.10','11':'NC_000011.9','12':'NC_000012.11',
                    '13':'NC_000013.10','14':'NC_000014.8','15':'NC_000015.9','16':'NC_000016.9',
                    '17':'NC_000017.10','18':'NC_000018.9','19':'NC_000019.9','20':'NC_000020.10',
                    '21':'NC_000021.8','22':'NC_000022.10','X':'NC_000023.10','Y':'NC_000024.9'}
    elif int(ref_genome) == 38:
        #GRCh38 from https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.38
        NCBI_IDS = {'1':'NC_000001.11','2':'NC_000002.12','3':'NC_000003.12','4':'NC_000004.12',
                    '5':'NC_000005.10','6':'NC_000006.12','7':'NC_000007.14','8':'NC_000008.11',
                    '9':'NC_000009.12','10':'NC_000010.11','11':'NC_000011.10','12':'NC_000012.12',
                    '13':'NC_000013.11','14':'NC_000014.9','15':'NC_000015.10','16':'NC_000016.10',
                    '17':'NC_000017.11','18':'NC_000018.10','19':'NC_000019.10','20':'NC_000020.11',
                    '21':'NC_000021.9','22':'NC_000022.11','X':'NC_000023.11','Y':'NC_000024.10'}
        

    try:        
        handle = Entrez.efetch(db="nucleotide", 
                               id=NCBI_IDS[str(chromosome)], 
                               rettype="fasta", 
                               strand=strand, #"1" for the plus strand and "2" for the minus strand.
                               seq_start=start,
                               seq_stop=end)
        record = SeqIO.read(handle, "fasta")
        handle.close()
        sequence = str(record.seq)
        return sequence
    except ValueError:
        print('ValueError: no sequence found in NCBI')
        return False
#####################################################
#####################################################
#####################################################