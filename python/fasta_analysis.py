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
  
  seq_list = [seq[i:i+3] for i in range(orf-1, len(seq)-1, 3)]
  
  if len(seq_list[-1]) < 3:
    seq_list = seq_list[:-1]
    
  return seq_list
  

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
  
  if len(starts) == 0 or len(stops) == 0:
    return []
  
  while stop_idx < len(stops):
    if stops[stop_idx] > starts[start_idx]:
      # print(starts[start_idx], stops[stop_idx])
      orfs.append({'start': start_idx*3+1, 'len': (stops[stop_idx]-starts[start_idx]+1)*3})
      bigger_starts = [i for i,x in enumerate(starts) if x > stops[stop_idx]]
      if len(bigger_starts) < 1:
        break
      else:
        start_idx = bigger_starts[0]
    else:
      stop_idx += 1
  
  return orfs


# TODO: explore mulit-thread processing 
# from concurrent.futures import ThreadPoolExecutor 
# from multiprocessing.dummy import Pool as ThreadPool


def build_orf_table(d, orf = 1):
  
  orfs = {}
  
  for k,v in d.items():
    seq_list = chop_seq(v, orf = orf)
    orfs[k] = find_orfs(seq_list)
  
  return(orfs)

def test_fun(x): 
  
  print(x)

  
def describe_orfs(filename = 'data/dna.example.fasta', 
                  orf = 1):
  
  d = read_fasta(filename)
  t = build_orf_table(d, orf = orf)
  
  # what is the length of the longest ORF in the file? 
  max_len = 0
  for k,v in t.items():
    if len(v) >= 1:
      lens = [x['len'] for x in v]
      cur_max = max(lens)
      if cur_max > max_len:
        max_len = cur_max
  
  # What is the identifier of the sequence containing the longest ORF? 
  for k,v in t.items():
    if max_len in [x['len'] for x in v]:
      print('Sequence', k, 'has longest orf at length:', max_len)
  
  
  # For a given sequence identifier, what is the longest ORF contained 
  # in the sequence represented by that identifier? 
  seq_id = 'gi|142022655|gb|EQ086233.1|16'
  
  for k,v in t.items():
    if seq_id in k:
      lens = [x['len'] for x in v]
      print('longest sequence in gene id:', seq_id, 'is:', max(lens))
  
  
  # What is the starting position of the longest ORF in the sequence that contains it? 
  # The position should indicate the character number in the sequence. For instance, \
  # the following ORF in reading frame 1: -- slightly modified
  
  for k,v in t.items():
    orfs = v
    for orf in orfs:
      if orf['len'] == max_len:
        print(k)
        print('The starting possition of the longest orf in ', 
              'the file is:', orf['start'])




#==========================================================================================
# Repeats

import re

def describe_repeats(filename = 'data/dna.example.fasta', n = 3):
   
  d = read_fasta(filename)
  
  seq = ''.join(d.values())
  
  unique_seqs = set([seq[i:i+n] for i in range(len(seq)-n+1)])
  
  seq_counts = {}
  
  
  for unique_seq in unique_seqs:
    matches = re.finditer(r'(?=(' + unique_seq + '))', seq)
    matches = [match.group(1) for match in matches]
    nmatches = len(matches)
    if nmatches > 1:
      if nmatches in seq_counts.keys(): 
        seq_counts[nmatches] = seq_counts[nmatches] + unique_seq
      else: 
        seq_counts[nmatches] = unique_seq
  
  max_count = max([x for x in seq_counts.keys()])
  print('Most repeating sequence found in', filename, 'is', 
       seq_counts[max_count], 'with', max_count, 'repeats')
       
  print(len(seq_counts[max_count]))




