from global_alignment import global_alignment
from utils import get_search_space_left, get_search_space_right, split_string
import re


def lcp_array(comb, suffix_a):
    """
    Purpose:
        generates the lcp array

    Args:
        comb (str): combination of both strings with "$" appended to the  (s + "$" + t + "#")
        suffix_a (array): suffix array 

        suffix array 
        - all the suffixes in sorted order with the index in which they occur 

    Returns:
        array: the lcp array 
    """

    rank = [0] * len(comb)
    lcp = [0] * len(comb)
    
    for i, suffix in enumerate(suffix_a):
        rank[suffix] = i
    
    idx = 0
    for i in range(len(comb)):
        if rank[i] > 0:
            j = suffix_a[rank[i] - 1]
            while (i + idx < len(comb)) and (j + idx < len(comb)) and (comb[i + idx] == comb[j + idx]):
                idx += 1
            lcp[rank[i]] = idx
            if idx > 0:
                idx -= 1
    return lcp


def lus(s, t, k1):
    """
    Purpose:
        find optimal unique common shared subsequence 
        only returns one substring 
    Args:
        s: string S 
        t: string T
        k1: int (such that k1 < k2 -- see tandem_transform_lcs)

    Returns: 
        the longest common unique subsequence between strings S and T. 
        if there are two longest common subsequences of the same length, 
        return the one which maximizes the global alignment score.
        
        our code breaks if there are no matches between S and T 
        
    """

    # find longest unique substrings between s and t
    comb = s + "$" + t + "#"

    # craete suffix array
    suffixes = [(comb[i:], i) for i in range(len(comb))]
    suffixes.sort()
    suffix_array = [suffix[1] for suffix in suffixes]
    
    # create lcp array
    lcp = lcp_array(comb, suffix_array)
    
    # find all unique substrings
    all_shared_substrings = []

    for i in range(1, len(comb)):
        suff_one = suffix_array[i]
        suff_two = suffix_array[i - 1]

        if (suff_one < len(s) and suff_two > len(s)) or (suff_one > len(s) and suff_two < len(s)):
            if lcp[i] > 0:
                if (i-1) >= 0:
                     prev = lcp[i - 1] 
                else:
                     prev = 0
                if (i+1) < len(lcp):
                     next = lcp[i+1]
                else:
                     next = 0

                if lcp[i] > max(prev, next):
                    shared = comb[suff_one : suff_one + lcp[i]]
                    all_shared_substrings.append(shared)
    
    # find longest string with most optimal score
    max_len = 0
    for string in all_shared_substrings:
        max_len = max(len(string), max_len)

    longest = None
    cur_score = float("-inf")

    for string in all_shared_substrings:
        if len(string) == max_len and len(string) > k1:
            sl, sr, tl, tr = split_string(s, t, string)
            score = global_alignment(sl, tl)[0] + global_alignment(sr, tr)[0]
            if score > cur_score:
                longest = string
    return longest


def tandem_transform_lcs(s: str, t: str, k1: int, k2: int, v: int):
        """
        Purpose:
                Finds longest unique substring between s and t and greedily applies most optimal tandem transformation
                on left and right flank of LUS.
        Args:
                s: string S 
                t: string T
                k1: int (such that k1 < k2 -- see tandem_transform_lcs)
                k2: int (longest possible kmer length)
                v: similarity threshold
        
        Returns: 
                score[0]: alignment score
                score[1]: globally aligned s string
                score[2]: globally aligned t string
                kmers: kmers used in tandem transformaation
                common_alignment: lus used
        """
    
        # keep track of kmers
        kmers = []
        max_iters = 100
        iters = 0

        # This while loop ensures that all potential tandem transformations around alignment are explored
        while iters < max_iters:

                score = global_alignment(s, t)
                if score[0] >= v:
                        if iters == 0:
                                print("Strings similar enough. No transformations needed")
                        return (score[0], score[1], score[2], kmers, "")

                # find best possible alignment between s and t longer than k1, uses greedy approach
                common_alignment = lus(s, t, k1)
                
                # return if no common alignments exist
                if not common_alignment:
                        return (score[0], score[1], score[2], kmers, common_alignment)
                
                indices = [(x.start(), x.end() - 1) for x in re.finditer(re.escape(common_alignment), s)][0]

        
                # align s and t to common prefix, get the left flanking and right flanking sequence
                sl, sr, tl, tr = split_string(s, t, common_alignment)

                # left side applied transformations

                # left kmer search space
                search_space_l = get_search_space_left(sl, tl, k1, k2, common_alignment)

                # test left side search space values and find the one that leads to highest global alignment score when applied
                max_score = float("-inf")
                target_kmer = None
                transformation = ""
                for transforms in search_space_l:
                        if transforms[-1] == "del":
                                s_trans = sl[:len(sl) - len(transforms[0])]
                        else:
                                s_trans = sl + transforms[0]
                        reconstructed = s_trans + common_alignment  + s[indices[-1] + 1:]

                        # apply transformation and check alignment score
                        score = global_alignment(reconstructed, t)

                        if score[0] > v:
                                # if threshold is reached, a tandem transformation has been found
                                kmers.append(transforms)
                                return (score[0], score[1], score[2], kmers, common_alignment)
                        else:
                                # assign the highest scoring transformation 
                                if score[0] > max_score:
                                        transformation = reconstructed
                                        target_kmer = transforms

                # right kmer search space
                if search_space_l:
                        s = transformation
                        kmers.append(target_kmer)
                indices = [(x.start(), x.end() - 1) for x in re.finditer(re.escape(common_alignment), s)][0]
                search_space_r = get_search_space_right(sr, tr, k1, k2, common_alignment)


                # test right side search space values and select the one that leads to highest global alignment score when applied
                transformation_r = ""
                for transforms in search_space_r:

                        if transforms[-1] == "del":
                                s_trans = sr[len(transforms[0]):]
                        else:
                                s_trans = transforms[0] + sr
        
                        reconstructed = s[:indices[-1] + 1] + s_trans

                        # apply transformation and check alignment score
                        score = global_alignment(reconstructed, t)

                        if score[0] > v:
                                # if threshold is reached, a tandem transformation has been found
                                kmers.append(transforms)
                                return (score[0], score[1], score[2], kmers, common_alignment)
                        else:
                                # assign the highest scoring transformation 
                                if score[0] > max_score:
                                        transformation_r = reconstructed
                                        target_kmer = transforms

                # transformation on right side complete, final s transformation for this iteration
                if search_space_r:
                        s = transformation_r
                        kmers.append(target_kmer)

                # no more possible tandem duplications or deletions
                if not search_space_l and not search_space_r:
                        #print("All tandem transformations exhausted for this alignment.")
                        return (score[0], score[1], score[2], kmers, common_alignment)
                
                iters += 1
                
        return 
