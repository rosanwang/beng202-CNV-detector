# global alignment functions, utilizing dynamic programming
# Ch 5 implementation

def GlobalAlignmentBacktrack(v, w, indel, match_reward, mismatch_penalty):
    vLen = len(v)
    wLen = len(w)
    backtrack = [[0] * (vLen + 1) for _ in range(wLen + 1)]
    s = [[0] * (vLen + 1) for _ in range(wLen + 1)]

    for i in range(1, wLen + 1):
        s[i][0] = s[i-1][0] - indel
        backtrack[i][0] = "down"
    for j in range(1, vLen + 1):
        s[0][j] = s[0][j-1] - indel
        backtrack[0][j] = "right"

    for i in range(1, wLen + 1):
        for j in range(1, vLen + 1):
            w_letter = w[i-1]
            v_letter = v[j-1]

            # determine if match or mismatch
            if w_letter == v_letter:
                diag = s[i-1][j-1] + match_reward
            else:
                diag = s[i-1][j-1] - mismatch_penalty

            right = s[i][j-1] - indel
            down = s[i-1][j] - indel

            s[i][j] = max(right, down, diag)

            if s[i][j] == down:
                backtrack[i][j] = "down"
            elif s[i][j] == right:
                backtrack[i][j] = "right"
            elif s[i][j] ==  diag:
                backtrack[i][j] = "diagonal"

    return backtrack, s

def OutputGlobalAlignment(Backtrack, str2, str1):
    align1 = ""
    align2 = ""
    i = len(str2)
    j = len(str1)

    while not(i == 0 and j == 0):
        if Backtrack[i][j] == "down": # insertion in string 1 
            i = i-1
            align2 = str2[i] + align2
            align1 = "-" + align1

        elif Backtrack[i][j] == "right": # deletion in string 1 
            j = j-1
            align2 = "-" + align2
            align1 = str1[j] + align1
        else:
            i = i-1
            j = j-1
            align2 = str2[i] + align2
            align1 = str1[j] + align1

    return align1, align2

# match_reward=match, mismatch_penalty=mismatch, indel_penalty=gap
def global_alignment( s: str, t: str, match_reward=2, mismatch_penalty=1, indel_penalty=2):
        backtrack, score = GlobalAlignmentBacktrack(s, t, indel_penalty, match_reward, mismatch_penalty)
        alignments = OutputGlobalAlignment(backtrack, t, s)
        return score[len(t)][len(s)], alignments[0], alignments[1]