from tandem_transform import tandem_transform_lcs
from utils import split_string
from global_alignment import global_alignment
from copy import deepcopy
import sys

def run_lcp(s: str, t: str, k1: int, k2: int, v: int, depth: int = 0, max_depth: int = 5, all_kmers=None):
    """
    Purpose:
        Uses divide and conquer to recursively run tandem_transform_lcs
    
    Args:
        s: String s
        t: String t
        k1: int (such that k1 < k2 -- see tandem_transform_lcs)
        k2: int (longest possible kmer length)
        v: similarity threshold 
        depth: Current recursion depth
        max_depth: maximum recursion depth to prevent infinite recursion
        all_kmers: List of kmers used 
        
    Returns:
        final_score: final alignment score
        final_s: final s aligned string
        final_t: final t aligned string
        all_kmers: all kmers used in transformation
    """

    if all_kmers is None:
        all_kmers = []
    
    # Base case: max recursion depth reached
    if depth >= max_depth:
        print("No tandem transformation possible.")
        return

    # run initial tandem transform 
    result = tandem_transform_lcs(s, t, k1, k2, v)

    # Update kmers_used
    if result[3]:
        all_kmers.extend(result[3])

    # base case: s and t are past similarity threshold
    if result[0] > v:
        return result[0], result[1], result[2], all_kmers
    
    # second base case: no possible tandem transformations, as s and t below smallest possible kmer transformation
    if len(s) < k1 and len(t) < k1:
        return result[0], result[1], result[2], all_kmers
    
    # third base case: no common shared subsequence
    lus = result[-1]
    if not lus:
        return result[0], result[1], result[2], all_kmers
    
    # begin divide and conquer approach
    s = result[1].replace("-","")
    t = result[2].replace("-","")
    s_left, s_right, t_left, t_right = split_string(s,t,lus)

    # keep track of scores
    final_score = result[0] 
    final_s = result[1]
    final_t = result[2]
    final_kmers = deepcopy(all_kmers)

    # recursively call for left side
    left_side = None
    if s_left and t_left and (len(s_left) >= k1 or len(t_left) >= k1):

        left_side = run_lcp(s_left, t_left, k1, k2, v, depth + 1, max_depth, deepcopy(all_kmers))
        s_trans = left_side[1].replace("-", "") + lus + s_right
        t_trans = left_side[2].replace("-", "") + lus + t_right

        # did left transform result in an optimal global alignment score?
        left_score, ls, lt = global_alignment(s_trans, t_trans)
            
        if left_score > final_score:
            final_score = left_score
            final_s = ls
            final_t = lt
            final_kmers = left_side[3]

    # recursively call for right side
    right_side = None
    if s_right and t_right and (len(s_right) >= k1 or len(t_right) >= k1):
        right_side = run_lcp(s_right, t_right, k1, k2, v, depth + 1, max_depth, deepcopy(all_kmers))

        s_trans = s_left + lus + right_side[1].replace("-", "")
        t_trans = t_left + lus + right_side[2].replace("-", "")
            
        # did right transform result in a optimal global alignment score?
        right_score, rs, rt = global_alignment(s_trans, t_trans)
            
        if right_score > final_score:
            final_score = right_score
            final_s = rs
            final_t = rt
            final_kmers = right_side[3]

    return final_score, final_s, final_t, final_kmers


if __name__ == "__main__":
    # which test case do you want to use?
    i = 15
    # homedir = "/Users/karenpu/Desktop/beng202-CNV-detector/"
    homedir = "/Users/rosanwang/Documents/school/ucsd/year 1/winter/BENG 202/beng202-CNV-detector"

    # with open(f"{homedir}/test_cases/test_{i}.txt", "r") as f:
    #     inputs = f.readlines()
    with open(f"{homedir}/test_cases/inputs/input_{i}.txt", "r") as f:
        inputs = f.readlines()

    print(f"Test case description:\n {inputs[0]}")

    s = inputs[1].replace("\n","")
    t = inputs[2].replace("\n","")
    k1 = inputs[3].replace("\n","")
    k2 = inputs[4].replace("\n","")
    v  = inputs[5].replace("\n","")
    final_score, s_final, t_final, kmers_final = run_lcp(s, t, int(k1), int(k2), int(v))

    # no possible tandem transformation if similarity threshold isn't reached at this point
    if final_score < int(v):
        print("Couldn't find a good tandem transformation, but here's the best attempt:\n")
        print(f"Final Score: {final_score}")
        print(f"Original Strings: \n{s}\n{t}")
        print(f"Aligned Transformation:\n{s_final}\n{t_final}")
        print(f"Kmers Used: {kmers_final}")
        print("=====================================\n")
    else:
        print("\n===== LCP Transformation Results =====")
        print(f"Final Score: {final_score}")
        print(f"Original Strings: \n{s}\n{t}")
        print(f"Aligned Transformation:\n{s_final}\n{t_final}")
        print(f"Kmers Used: {kmers_final}")
        print("=====================================\n")
