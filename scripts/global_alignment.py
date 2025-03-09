# Global Alignment Functions
import sys
import os
import math

def BackTrack(s, t, match_reward, mismatch_penalty, indel_penalty):
    backtrack, scores = dict(), dict()
    len_s, len_t = len(s) + 1, len(t) + 1

    scores[(0,0)] = 0
    backtrack[(0, 0)] = ""
    for i in range(1, len_s):
        scores[(i, 0)] = -indel_penalty + scores[(i-1,0)]
        backtrack[(i, 0)] = "del" 
    for j in range(1, len_t):
        scores[(0, j)] = -indel_penalty + scores[(0, j-1)]
        backtrack[(0, j)] = "ins" 
    for i in range(1, len_s): #going across
        for j in range(1, len_t):
            if s[i-1] == t[j-1]:
                match = match_reward
            else:
                match = -mismatch_penalty

            scores[(i, j)] = max( (scores[(i-1, j)] - indel_penalty), (scores[(i, j-1)]- indel_penalty), (scores[(i-1, j-1)] + match) )

            if scores[(i, j)] == (scores[(i-1, j)] - indel_penalty):
                backtrack[(i, j)] = "del" #down
            elif scores[(i, j)] == (scores[(i, j-1)] - indel_penalty):
                backtrack[(i, j)] = "ins"  #left
            elif scores[(i, j)] == (scores[(i-1, j-1)] + match):
                backtrack[(i, j)] = "m" #diagonal down
    return backtrack, scores

def Output(backtrack, s, t, i, j):
    if (i == 0) & (j == 0):
        return "", ""
    if backtrack[(i,j)] == "del": #deleted in reference to the top sequence s (missing from t, present in s)
        out = Output(backtrack, s, t, i-1, j)
        return out[0]+s[i-1], out[1]+"-"
    
    elif backtrack[(i,j)] == "ins":
        out = Output(backtrack, s, t, i, j-1)
        return out[0] + "-", out[1] + t[j-1]
    else: 
        out = Output(backtrack, s, t, i-1, j-1) 
        return out[0] + s[i-1], out[1] + t[j-1]

def global_alignment(match_reward, mismatch_penalty, indel_penalty, s, t): 
    backtrack, scores = BackTrack(s, t, match_reward, mismatch_penalty, indel_penalty)
    s_a, t_a = Output(backtrack, s, t, len(s), len(t))
    alignment_score = scores[len(s), len(t)]
    #return (alignment_score, s_a, t_a)
    return alignment_score


