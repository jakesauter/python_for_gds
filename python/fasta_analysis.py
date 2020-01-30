def read_fasta(filename):
  
  d = {}
  for line in open(filename):
    line = line.rstrip()
    if line.startswith('>'):
      key = line[1:]
      d[key] = ''
    else: 
      d[key] = d[key] + line
      
  return d


def length_of_seqs(d):
  return [len(x) for x in d]

def describe_longest(d, return_val = False):
  lens = [len(x) for x in d.values()]
  longest = max(lens)
  num_longest = sum([len(x) == longest for x in d.values()])
  
  print('There are %s sequences of length %s' % (num_longest, longest))
  
  longest_d = {}
  
  for k,v in d.items():
    if len(v) == longest: 
      print('gene %s has the longest length of %s' % 
      (k, longest))
      longest_d[k] = v 
      
  if return_val: 
    return longest_d
  
  
def describe_shortest(d, return_val = False): 
  lens = [len(x) for x in d.values()]
  shortest = min(lens)
  num_shortest = sum([len(x) == shortest for x in d.values()])
  
  print('There are %s sequences of length %s' % (num_shortest, shortest))
  
  shortest_d = {}
  
  for k,v in d.items():
    if len(v) == shortest: 
      print('gene %s has the shortest length of %s' % 
      (k, shortest))
      shortest_d[k] = v
  
  if return_val: 
    return shortest_d
    
    
# Given an input reading frame on the forward strand (1, 2, or 3) 
# your program should be able to identify all ORFs present in each 
# sequence of the FASTA file, and answer the following questions: 
#   
# what is the length of the longest ORF in the file? 
#   
# What is the identifier of the sequence containing the longest ORF? 
#   
# For a given sequence identifier, what is the longest ORF contained 
# in the sequence represented by that identifier? 
# 
# What is the starting position of the longest ORF in the sequence that contains it? 
# The position should indicate the character number in the sequence. For instance, \
# the following ORF in reading frame 1:

import sys

def chop_seq(seq, orf = 1): 
  if orf not in (1,2,3):
    sys.stderr.write('orf parameter must be one of (1,2,3), not ' + str(orf) + '!!')
    return()
  
  return [seq[i:i+3] for i in range(0, len(seq[orf-1:])-2, 3)]
  

# test: 
#
# 'ATGCGGTAGGACATGGCGCAGTAAATGCGCTACTGA'

def find_orfs(seq_list, 
              start_codon = 'ATG',
              stop_codons = ('TAA', 'TAG', 'TGA')):
  
  orfs = []
  
  starts = [i for i, x in enumerate(seq_list) if x == start_codon]
  stops = [i for i,x in enumerate(seq_list) if x in stop_codons]
  
  start_idx = 0
  stop_idx = 0
  
  while start_idx < len(starts):
    while stop_idx < len(stops):
      if stops[stop_idx] > starts[start_idx]:
        orfs.append({'start': start_idx, 'len': stops[stop_idx]-starts[start_idx]})
        bigger_starts = [i for i,x in enumerate(starts) if x > stops[stop_idx]]
        if len(bigger_starts) < 1:
          start_idx = len(starts)
          break
        start_idx = bigger_starts[0]
        stop_idx = stop_idx + 1
        break
      else: 
        stop_idx = stop_idx + 1
    start_idx = start_idx + 1
  
  
  return orfs
      

# from concurrent.futures import ThreadPoolExecutor 
# from multiprocessing.dummy import Pool as ThreadPool


def build_orf_table(d, orf = 1):

  orfs = {}

  for k,v in d.items():
    print(k)
    orfs[k] = find_orfs(chop_seq(v, orf = orf))
    
  return(orfs)


# def describe_orfs(orf = 1, filename):
# 
#   if orf not in (1,2,3):
#     sys.stderr.write('orf parameter must be one of (1,2,3), not ' + str(orf) + '!!')
#     return()
#     
#   # Now that we have a valid orf we need to gather all orfs from the file
  # while maintaining which sequences they come from 
  
  # In order to do this we will make a dictionary in which the keys are the fasta identifiers, 
  # and the entries are lists of orfs. The orf entries will be dictionaries themselves containing
  # a 'len' and 'start' key
  #
  #     { seq1: [{start:, len:}, ..., {start:, len:}]
  #
  #                  ...
  #
  #       seqn: [{start:, len:}, ..., {start:, len:}]
  #     }
  
  
  # fasta = read_fasta(filename)
  
  # orfs = {}
  
  # for k in fasta:
    
# 
# 
# def describe_fasta_file(filename): 
#   
#   fasta = read_fasta(filename)
#   
#   print('Number of fasta records identified in the', 
#   'provided file %s: %s' % (filename, len(fasta)))
