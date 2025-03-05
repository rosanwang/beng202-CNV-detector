# Affine Gap Alignment Implementation of CNV Detector

# Libraries
import os
import sys
import math
from affine_gap_alignment import BackTrack, Output, affine_gap_alignment


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
    for gap_set, idx in gap_indices: # consider every gap in s_gap
        gap_string = gap_alignments[opps[idx]] # we search opposite string for kmer
        for start, end in gap_set:
            gap = gap_string[start:end] 
            # kmer in k1, k2 at every index
            # we assume kmer must start at first index of gap # FIX 
            # consider that the kmer can be a repeat of a substring before and after gap, or can be a complete deletion (exists only in gap portion)
            # start of string
            if start < min(k1, k2): # if kmer can exist before 
                for k in range(k1, k2):
                    potential_kmer = gap[:k]
                    # matches previous substring --> keep going until repeats end (while in gap), or we hit end of gap
                    curr_idx = start - k
                    while (potential_kmer == gap_string[curr_idx:curr_idx+k]) & (curr_idx+k < end): # matches previous substring = hit
                        kmer_deletions[potential_kmer] +=1
                        curr_idx +=k
                    # LET ME COOK 

                    # TODO 
                    



                #check if kmer exist in substring directly prior to gap
                
    return kmers, s_prime, t_prime



# Examples
match_reward, mismatch_penalty, gap_opening_penalty, gap_extension_penalty = 1, 5, 2, 1
s = "CCAT"
t = "GAT"