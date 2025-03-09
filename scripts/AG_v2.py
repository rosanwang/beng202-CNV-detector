# Affine Gap Alignment Implementation of CNV Detector

# Libraries
import os
import sys
import math
import copy
from itertools import product
from affine_gap_alignment import affine_gap_alignment
from global_alignment import global_alignment

'''
Assumptions:
- a duplicated kmer must exist in the string (directly tandem) to the new insertion
- there are no SNPs in or outside of the CNV kmers
- there are no partial deletions/duplications
- s and t can be of different lengths (different windows), but are exactly the same apart from terminals and CNVs
'''

def get_gaps(gap_alignments: tuple):
    '''
    get the (start, end) indices of every gap in strings s and t
    '''
    gaps = []
    gap_start = None
    for gap_string in gap_alignments:
        if "-" not in gap_string:
            gaps.append([])
            break
        gap_indices = []
        for idx, character in enumerate(gap_string):
            if character == "-":
                if gap_start == None:
                    gap_start = idx
            else:
                if gap_start != None:
                    #gap_indices.append((gap_start, idx-1))
                    gap_indices.append((gap_start, idx))
                    gap_start = None

        if gap_start != None:
            #gap_indices.append((gap_start, len(gap_string) - 1))
            gap_indices.append((gap_start, len(gap_string)))

        gaps.append(gap_indices)

    return gaps

def tandem_forward(gap_string: str, potential_kmer: str, kmer_sets: dict, k: int, curr_idx: int, start: int, end: int):
    '''
    given a string s, finds all tandems copies of potential_kmer, left_to_right
    '''
    kmer_num_instances = 0
    #print('potential kmer', potential_kmer)
    #print(gap_string[curr_idx:curr_idx+k])
   
    while (potential_kmer == gap_string[curr_idx:curr_idx+k]) & (curr_idx+k <= end): # matches previous substring = hit
        kmer_num_instances +=1
        curr_idx +=k
    #kmer_sets[(start, end)].append((potential_kmer, kmer_num_instances)) # rem

    print("kmer_num_instances", kmer_num_instances)
    print("")
    return kmer_num_instances, curr_idx

def tandem_backward(gap_string: str, potential_kmer: str, kmer_sets: dict, k: int, curr_idx: int, start: int, end: int):
    '''
    given a string s, finds all tandems copies of potential_kmer, left_to_right
    '''
    kmer_num_instances = 0
    while (potential_kmer == gap_string[curr_idx-k:curr_idx]) & (curr_idx-k >= start): 
        kmer_num_instances +=1
        curr_idx -=k
    #kmer_sets[(start, end)].append((potential_kmer, kmer_num_instances)) # rem

    return kmer_num_instances, curr_idx


def possible_kmer_sets(gap_alignments, gap_indices, opps, k1, k2, s, t):
    '''
    across all gaps in gap alignments, gets all possible kmer repeat sets
    '''
    # gaps indicate a DELETION 
    # In each gap, for every kmer in length [k1, k2]
    #  {s: {(start, end): [(kmer1, # tandem instances), (kmer2, # tandem instances)], ... }, t:}
    #pair_kmer_sets = []
    kmer_sets = {}
    for string_idx, gap_set in enumerate(gap_indices): # consider every gap in s_gap
        print("gap set", gap_set)
        print("idx", string_idx)
        original_string = gap_alignments[string_idx]
        #gap_string = gap_alignments[opps[string_idx]] # we search opposite string for kmer
        if string_idx == 0:
            gap_string = s
        else:
            gap_string = t
        print("gap string", gap_string)
        for start, end in gap_set:
            gap_length = end - start
            print("start", start, "end", end)
            gap = gap_string[start:end] 
            print("gap", gap)

            # get number of "-" prior to this in the original string (which has the gap)
            num_dash_org = original_string[:start].count("-")
            org_start, org_end = start - num_dash_org, end - num_dash_org
            kmer_sets[(org_start, org_end)] = []

            if start >= k1: # if kmer can exist before 
                kmer_range = [x for x in range(k1, k2+1) if x <= gap_length]
                print("kmer range", kmer_range)
                #for k in range(k1, k2):
                for k in kmer_range:
                    print("k", k)
                    potential_kmer1 = gap[:k]
                    print("potential_kmer1", potential_kmer1)
                    # matches previous substring --> keep going until repeats end (while in gap), or we hit end of gap
                    print("prev", gap_string[start-k:start])
                    if potential_kmer1 == gap_string[start-k:start]:
                        #temp_kmer_sets = copy.deepcopy(temp_kmer_sets)
                        curr_idx = start 
                        kmer1_num_instances, curr_idx = tandem_forward(gap_string, potential_kmer1, kmer_sets, k, curr_idx, start, end)

                        # kmer fills gap / is valid
                        if (curr_idx == end):
                            kmer_sets[(org_start, org_end)].append([(org_start, potential_kmer1, kmer1_num_instances, string_idx)]) #keep track of which string is modified

                        # there is a second kmer in gap, reverse direction
                        elif (curr_idx != end) & (end - curr_idx >= k1):
                            for k in range(k1, k2):
                                kmer_start = curr_idx # where new kmer starts
                                potential_kmer2 = gap[curr_idx:curr_idx+k]
                                kmer2_num_instances, curr_idx = tandem_backward(gap_string, potential_kmer2, kmer_sets, k, curr_idx, start, end)
                                if curr_idx == end: # valid, add
                                    kmer_sets[(org_start, org_end)].append([(kmer_start, potential_kmer1, kmer1_num_instances, string_idx), (potential_kmer2, kmer2_num_instances, string_idx)])
            
            # check substring after gap
            if end <= len(gap_string) - k1:
                kmer_range = [x for x in range(k1, k2+1) if x <= gap_length]
                print("kmer range", kmer_range)
                for k in kmer_range:
                    potential_kmer1 = gap[end-k:]
                    # matches substring after --> keep going until repeats end (while in gap), or we hit end of gap
                    if potential_kmer1 == gap_string[end:end+k]:
                        #temp_kmer_sets = copy.deepcopy(temp_kmer_sets)
                        curr_idx = end 
                        kmer1_num_instances, curr_idx = tandem_backward(gap_string, potential_kmer1, kmer_sets, k, curr_idx, start, end)

                        # kmer fills gap / is valid
                        if (curr_idx == start):
                            kmer_sets[(org_start, org_end)].append([(org_start, potential_kmer1, kmer1_num_instances, string_idx)])

                        # there is a second kmer in gap, reverse direction
                        elif (curr_idx != start) & (curr_idx - start >= k1):
                            for k in range(k1, k2):
                                kmer_start = curr_idx # where new kmer starts
                                potential_kmer2 = gap[curr_idx-k:curr_idx]
                                kmer2_num_instances, curr_idx = tandem_forward(gap_string, potential_kmer2, kmer_sets, k, curr_idx, start, end)
                                if curr_idx == start: # valid, add
                                    kmer_sets[(org_start, org_end)].append([(kmer_start, potential_kmer1, kmer1_num_instances, string_idx), (potential_kmer2, kmer2_num_instances, string_idx)])
   
    return kmer_sets  

def all_kmer_combinations(kmer_sets):
    '''
    across all gaps in gap alignments, gets all combinations of cnv possibilities
    '''
    kmer_combinations = [[]]
    
    for key, lists_for_key in kmer_sets.items():
        new_combinations = []
        
        for combo in kmer_combinations:
            # don't pick any list from this key (keep combination as is)
            new_combinations.append(combo)
            
            # or, pick one of the lists from this key and add to the combination
            for lst in lists_for_key:
                #key_lst = [key + lst[0]] # keep track of gap
                new_combinations.append(combo + lst)
        
        kmer_combinations = new_combinations
    
    # Remove the empty combination
    if kmer_combinations and len(kmer_combinations[0]) == 0:
        kmer_combinations = kmer_combinations[1:]
    
    return kmer_combinations

def v_threshold_mod_combinations(s: str, t: str, v: int, all_mod_combinations):
    '''
    given all possible kmer modifications across all gaps, if meets v threshold for global alignment, gets the modified strings
    '''
    v_thresh_mods = [] #[(s_mod, num_mods)]

    #for combintation in all_mod_combinations:
    # we ONLY modify s

    min_num_mods = 1000000
    best_mod_s = s
    best_kmers = []

    for combination in all_mod_combinations:
        mod_s = (s, t)
        num_mods = 0
        kmers = []
        for modification in combination:
            kmer_start, kmer, num_repeats, gap_string = modification[0], modification[1], modification[2], modification[3]
            rep_kmer = kmer*num_repeats
            num_mods += num_repeats
            if gap_string == 0: # s is the string where the gap is (insert into s)
                mod_s = s[:kmer_start] + rep_kmer + s[kmer_start:] #FIX
                kmers.append((kmer, "ins"))
                print("mod_s insertion", mod_s)
            if gap_string == 1: # the gap is in t --> delete from s
                mod_s = s[:kmer_start] + s[kmer_start+ len(rep_kmer):]
                print("mod_s deletion", mod_s)
                kmers.append((kmer, "del"))
        
        # calculate global alignment for this combination
        print("mod_s", mod_s)
        global_alignment_score, aligned_s, aligned_t = global_alignment(match_reward, mismatch_penalty, 1, mod_s, t)
        if global_alignment_score <= v:
            if num_mods < min_num_mods:
                best_mods = (aligned_s, aligned_t)
                best_kmers = kmers

    return best_mods, best_kmers

def CNV_detector(s: str, t: str, v: int, k1: int, k2: int, match_reward: int, mismatch_penality: int, gap_opening_penalty: int, gap_closing_penalty: int):
    '''
    s, t: input strings, not necessarily the same length
    v: v-proper global similarity threshold
    k1, k2: range of possible kmer lengths
    '''

    # Affine Gap Alignment 
    gap_alignments = affine_gap_alignment(match_reward, mismatch_penalty, gap_opening_penalty, gap_extension_penalty, s, t)
    print("gap alignments", gap_alignments)
    print("")
    # ('-ACCATCTT-', 'TACCATCATC')

    # get the number of "-" in each alignment

    # Get positions in s_gap and t_gap where gaps occur --> these are where CNVs may exist
    gap_indices = get_gaps(gap_alignments)
    opps = [1,0]
    print("gap_indices", gap_indices)
    print("")
    # [[(0, 1), (9, 10)], []]

    possible_mods = possible_kmer_sets(gap_alignments, gap_indices, opps, k1, k2, s, t)
    print("possible_mods", possible_mods)
    print("")
                    
    all_mod_combinations = all_kmer_combinations(possible_mods)
    print("all_mod_combinations", all_mod_combinations)

    # check which combinations meet threshold global alignment score >= v
    # get alterated s that requires minimum number of modifications
    #transform t accordingly
    best_mod_s, best_kmers = v_threshold_mod_combinations(s, t, v, all_mod_combinations)

    return  best_kmers, best_mod_s

# Examples
match_reward, mismatch_penalty, gap_opening_penalty, gap_extension_penalty = 0, 10, 5, 1
s = "ACCATCTT"
t = "TACCATCATC"
k1, k2 = 2, 3
v = 6

x = CNV_detector(s, t, v, k1, k2, match_reward, mismatch_penalty, gap_opening_penalty, gap_extension_penalty)
print(x)

