# Affine Gap Alignment Implementation of CNV Detector

# Libraries
import os
import sys
import math
import copy
from affine_gap_alignment import BackTrack, Output, affine_gap_alignment

'''
Assumptions:
- a duplicated kmer must exist in the string (directly tandem) to the new insertion
- there are no SNPs in or outside of the CNV kmers
- there are no partial deletions/duplications
- s and t can be of different lengths (different windows), but are exactly the same apart from terminals and CNVs
'''


def get_gaps(gap_alignments):
    '''
    get the (start, end) indices of every gap in strings s and t
    '''
    gaps = []
    gap_start = None
    for gap_string in gap_alignments:
        gap_indices = []
        for idx, character in enumerate(gap_string):
            if character == "-":
                if gap_start == None:
                    gap_start = idx
            else:
                if gap_start != None:
                    gap_indices.append(gap_start, idx-1)
                    gap_start = None

        if gap_start != None:
            gap_indices.append((gap_start, len(gap_string) - 1))

        gaps.append(gap_indices)

    return gaps

def tandem_forward(gap_string: str, potential_kmer: str, kmer_sets: dict, k: int, start: int, end: int):
    '''
    given a string s, finds all tandems copies of potential_kmer, left_to_right
    '''

    kmer_num_instances = 0
    while (potential_kmer == gap_string[curr_idx:curr_idx+k]) & (curr_idx+k < end): # matches previous substring = hit
        kmer_num_instances +=1
        curr_idx +=k
    kmer_sets[(start, end)].append((potential_kmer, kmer_num_instances))

    return kmer_sets, curr_idx

def tandem_backward(gap_string: str, potential_kmer: str, kmer_sets: dict, k: int, start: int, end: int):
    '''
    given a string s, finds all tandems copies of potential_kmer, left_to_right
    '''

    kmer_num_instances = 0
    while (potential_kmer == gap_string[curr_idx-k:curr_idx]) & (curr_idx-k > start): 
        kmer_num_instances +=1
        curr_idx -=k
    kmer_sets[(start, end)].append((potential_kmer, kmer_num_instances))

    return kmer_sets, curr_idx

#def recursive_forward(k1, k2, curr_idx gap_string, potential_kmer, kmer_sets, start, end):
    '''
    starting at some determined kmer
    '''
    # base case: space left is less than k1 (not a possible kmer)
    #if (end - curr_idx) < k1: # not a good set

        





def CNV_detector(s: str, t: str, v: int, k1: int, k2: int, match_reward: int, mismatch_penality: int, gap_opening_penalty: int, gap_closing_penalty: int):
    '''
    s, t: input strings, not necessarily the same length
    v: v-proper global similarity threshold
    k1, k2: range of possible kmer lengths
    ur mom
    '''
    
    kmers, s_prime, t_prime = [], s, t
    
    # Affine Gap Alignment 
    gap_alignments = affine_gap_alignment(s, t, match_reward, mismatch_penalty, gap_opening_penalty, gap_extension_penalty)

    # Get positions in s_gap and t_gap where gaps occur --> these are where CNVs may exist
    gap_indices = get_gaps(gap_alignments)
    opps = [1,0]
    print("gap_indices", gap_indices)

    # gaps indicate a DELETION 
    # In each gap, for every kmer in length [k1, k2]
    kmer_deletions = {} 
    kmer_sets = {} #  gap start, end: (kmer, number tandem instances)
    for gap_set, idx in gap_indices: # consider every gap in s_gap
        gap_string = gap_alignments[opps[idx]] # we search opposite string for kmer
        for start, end in gap_set:
            gap = gap_string[start:end] 
            kmer_sets[(start, end)] = []
            # kmer in k1, k2 at every index
            # we assume kmer must start at first index of gap # FIX 
            # consider that the kmer can be a repeat of a substring before and after gap, or can be a complete deletion (exists only in gap portion)
            # start of string
            # we need to consider when multiple CNVs co-exist in the same gap
            # keep track of sets across all gaps 
            if start < k1: # if kmer can exist before 
                for k in range(k1, k2):
                    potential_kmer = gap[:k]
                    '''
                    checking gap modulo k == 0 will not work here because  kmers of different lengths can coexist in a single gap
                    '''
                    # matches previous substring --> keep going until repeats end (while in gap), or we hit end of gap
                    if potential_kmer == gap_string[start-k:k]:
                        temp_kmer_sets = copy.deepcopy(temp_kmer_sets)
                        curr_idx = start 
                        temp_potential_kmer = potential_kmer
                        temp_kmer_sets, curr_idx = tandem_forward(gap_string, temp_potential_kmer, temp_kmer_sets, k, start, end)
                        if curr_idx != end:
                            for k in range(k1, k2):
                                kmer_sets, curr_idx = tandem_backward(gap_string, potential_kmer, kmer_sets, k, start, end)

                                temp_kmer_sets, curr_idx = tandem_forward(gap_string, temp_potential_kmer, temp_kmer_sets, k, start, end)
                            #check if length left is less than k1
                            temp_potential_kmer = 
                            if (end - curr_idx) < k1: #not a possible kmer combination
                                break
                        if curr_idx == end:
                            kmer_sets = temp_kmer_sets
                            
                    else: #we assume the kmer must exist after the gap
                        # traverse backwards
                        curr_idx = end
                        while curr_idx != start:
                            kmer_sets, curr_idx = tandem_backward(gap_string, potential_kmer, kmer_sets, k, start, end)


                    
                    
                    



                    
                    #kmer_sets[(start, end)].append(potential_kmer, kmer_num_instances)  # alone instances


                    # Consider flanks


                    # TODO # LET ME COOK
                    



                #check if kmer exist in substring directly prior to gap
                
    return kmers, s_prime, t_prime



# Examples
match_reward, mismatch_penalty, gap_opening_penalty, gap_extension_penalty = 1, 5, 2, 1
s = "CCAT"
t = "GAT"


'''
while (potential_kmer == gap_string[curr_idx: curr_idx+k]) & (curr_idx+k < end):
                        kmer_num_instances +=1
                        curr_idx+=k
                        
'''

'''
def CNV_detector(s: str, t: str, v: int, k1: int, k2: int, match_reward: int, mismatch_penality: int, gap_opening_penalty: int, gap_closing_penalty: int):
    '''
    s, t: input strings, not necessarily the same length
    v: v-proper global similarity threshold
    k1, k2: range of possible kmer lengths
    ur mom
    '''
    
    kmers, s_prime, t_prime = [], s, t
    
    # Affine Gap Alignment 
    gap_alignments = affine_gap_alignment(s, t, match_reward, mismatch_penalty, gap_opening_penalty, gap_extension_penalty)

    # Get positions in s_gap and t_gap where gaps occur --> these are where CNVs may exist
    gap_indices = get_gaps(gap_alignments)
    opps = [1,0]
    print("gap_indices", gap_indices)

    # gaps indicate a DELETION 
    # In each gap, for every kmer in length [k1, k2]
    kmer_deletions = {} 
    kmer_sets = {} #  gap start, end: (kmer, number tandem instances)
    for gap_set, idx in gap_indices: # consider every gap in s_gap
        gap_string = gap_alignments[opps[idx]] # we search opposite string for kmer
        for start, end in gap_set:
            gap = gap_string[start:end] 
            kmer_sets[(start, end)] = []
            # kmer in k1, k2 at every index
            # we assume kmer must start at first index of gap # FIX 
            # consider that the kmer can be a repeat of a substring before and after gap, or can be a complete deletion (exists only in gap portion)
            # start of string
            # we need to consider when multiple CNVs co-exist in the same gap
            # keep track of sets across all gaps 
            if start < min(k1, k2): # if kmer can exist before 
                for k in range(k1, k2):
                    potential_kmer = gap[:k]
                    '''
                    checking gap modulo k == 0 will not work here because  kmers of different lengths can coexist in a single gap
                    '''
                    # matches previous substring --> keep going until repeats end (while in gap), or we hit end of gap
                    if potential_kmer == gap_string[start-k:k]:
                        curr_idx = start 
                        while curr_idx != end:
                            
                            kmer_sets, curr_idx = tandem_forward(gap_string, potential_kmer, kmer_sets, k, start, end)

                    else: #we assume the kmer must exist after the gap
                        # traverse backwards
                        curr_idx = end
                        while curr_idx != start:
                            kmer_sets, curr_idx = tandem_backward(gap_string, potential_kmer, kmer_sets, k, start, end)


                    
                    
                    



                    
                    #kmer_sets[(start, end)].append(potential_kmer, kmer_num_instances)  # alone instances


                    # Consider flanks


                    # TODO # LET ME COOK
                    



                #check if kmer exist in substring directly prior to gap
                
    return kmers, s_prime, t_prime

'''