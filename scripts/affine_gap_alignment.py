# Affine Gap Alignment Functions
import sys
import os
import math

def BackTrack(s, t, match_reward, mismatch_penalty, gap_opening_penalty, gap_extension_penalty):
    lower_backtrack, middle_backtrack, upper_backtrack = {(0,0): ""}, {(0,0): ""}, {(0,0): ""}
    lower_scores, middle_scores, upper_scores = {(0,0): 0}, {(0,0): 0}, {(0,0): 0} # set middle (0,0) to 0
    len_s = len(s) + 1
    len_t = len(t) + 1
    #first row of lower is all zero score
    for j in range(0, len_t):
        #lower_scores[(0, j)] = -gap_opening_penalty # set 0
        lower_scores[(0, j)] = 0 # not even needed
        lower_backtrack[(0, j)] = "" #these dont need backtrack, would just be a start
    #first col of upper is all zero score
    for i in range(0, len_s):
        #upper_scores[(i, 0)] = -gap_opening_penalty #0
        upper_scores[(i, 0)] = 0 # not even needed
        upper_backtrack[(i, 0)] = "" #these dont need backtrack, would just be a start
    ########## now we can update the rest of the (i,j)s for i in range(1, len_s), j in range(1, len_t)
    for i in range(0, len_s): #going across first
        for j in range(0, len_t):

            if s[i-1] == t[j-1]:
                match = match_reward
            else:
                match = -mismatch_penalty
            if i != 0:
                if i == 1: 
                    lower_scores[(i,j)] = middle_scores[(i-1,j)] - gap_opening_penalty
                    lower_backtrack[(i, j)] = "m" 
                else:
                    lower_scores[(i,j)] = max( (lower_scores[(i-1,j)] - gap_extension_penalty), (middle_scores[(i-1,j)] - gap_opening_penalty) )

                    if lower_scores[(i,j)] == (lower_scores[(i-1,j)] - gap_extension_penalty):
                        lower_backtrack[(i, j)] = "del" #down
                    elif lower_scores[(i,j)] == (middle_scores[(i-1,j)] - gap_opening_penalty):
                        lower_backtrack[(i, j)] = "m" 
            if j != 0:
                if j == 1:
                    upper_scores[(i,j)] = middle_scores[(i,j-1)] - gap_opening_penalty
                    upper_backtrack[(i, j)] = "m"

                else:
                    upper_scores[(i,j)] = max( (upper_scores[(i,j-1)] - gap_extension_penalty), (middle_scores[(i,j-1)] - gap_opening_penalty) )

                    if upper_scores[(i,j)] == (upper_scores[(i,j-1)] - gap_extension_penalty):
                        upper_backtrack[(i, j)] = "ins" #left
                    elif upper_scores[(i,j)] == (middle_scores[(i,j-1)] - gap_opening_penalty):
                        upper_backtrack[(i, j)] = "m"
            #update middle  middle(i,j) must be calculated AFTER lower(i,j) and upper(i,j)
            if (i, j) != (0,0):
                if i == 0: #only consider upper
                    middle_scores[(i,j)] = upper_scores[(i,j)] # 0 weight edge
                    middle_backtrack[(i, j)] = "ins"

                elif j == 0: #only consider lower
                    middle_scores[(i,j)] = lower_scores[(i,j)]
                    middle_backtrack[(i, j)] = "del"

                else:
                    middle_scores[(i,j)] = max( lower_scores[(i,j)] , (middle_scores[(i-1,j-1)] + match), upper_scores[(i,j)] )

                    if middle_scores[(i,j)] == lower_scores[(i,j)]:
                        middle_backtrack[(i, j)] = "del" #down
                    elif middle_scores[(i,j)] == upper_scores[(i,j)]:
                        middle_backtrack[(i, j)] = "ins" #left    #switched order
                    elif middle_scores[(i,j)] == (middle_scores[(i-1,j-1)] + match):
                        middle_backtrack[(i, j)] = "m"
                
    return lower_backtrack, middle_backtrack, upper_backtrack, lower_scores, middle_scores, upper_scores

def Output(lower_backtrack, middle_backtrack, upper_backtrack, curr_back, s, t, i, j):

    if (i == 0) & (j == 0): # keep  #FIX
        return "", ""
    
    if curr_back[(i,j)] == "del": #down, go to lower_backtrack
        if curr_back == middle_backtrack:
            out = Output(lower_backtrack, middle_backtrack, upper_backtrack, lower_backtrack, s, t, i, j)
            return out
        else: 
            out = Output(lower_backtrack, middle_backtrack, upper_backtrack, lower_backtrack, s, t, i-1, j)
            return out[0]+ s[i-1], out[1]+"-"
    
    elif curr_back[(i,j)] == "ins": #left, go to upper_backtrack
        if curr_back == middle_backtrack: #we are in middle, and backtracking to upper (gap closing)
            out = Output(lower_backtrack, middle_backtrack, upper_backtrack, upper_backtrack, s, t, i, j)
            return out
        else: # we are in upper and backtracking middle 
            out = Output(lower_backtrack, middle_backtrack, upper_backtrack, upper_backtrack, s, t, i, j-1)
            return out[0]+ "-", out[1] + t[j-1]

    else: 
        if curr_back == middle_backtrack: 
            out = Output(lower_backtrack, middle_backtrack, upper_backtrack, middle_backtrack, s, t, i-1, j-1)
            return out[0] + s[i-1], out[1] + t[j-1]
        elif curr_back == lower_backtrack: 
            out = Output(lower_backtrack, middle_backtrack, upper_backtrack, middle_backtrack, s, t, i-1, j)
            return out[0]+s[i-1], out[1]+"-"
        else:
            out = Output(lower_backtrack, middle_backtrack, upper_backtrack, middle_backtrack, s, t, i, j-1)
            return out[0]+ "-", out[1] + t[j-1]


def affine_gap_alignment(match_reward, mismatch_penalty, gap_opening_penalty, gap_extension_penalty, s, t):
    lower_backtrack, middle_backtrack, upper_backtrack, lower_scores, middle_scores, upper_scores = BackTrack(s, t, match_reward, mismatch_penalty, gap_opening_penalty, gap_extension_penalty)
    alignment_score = middle_scores[(len(s), len(t))]  
    curr_back = middle_backtrack
    s_a, t_a = Output(lower_backtrack, middle_backtrack, upper_backtrack, curr_back, s, t, len(s), len(t))
    #return (alignment_score, s_a, t_a)
    return (s_a, t_a)

