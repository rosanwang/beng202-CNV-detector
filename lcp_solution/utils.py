
def split_string(s, t, lus):
    """
    Purpose:
        Helper function: splits s and t after common subsequence alignment
    Args:
        s: string S 
        t: string T
        lus: identified lus
        
    Returns: 
        s_left: subsequence to left of lus in s
        s_right: subsequence to right of lus in s
        t_left: subsequence to left of lus in t
        t_right: subsequence to right of lus in t
    """

    s_idx = s.find(lus)
    t_idx = t.find(lus)

    # Extract left and right parts
    s_left, s_right = s[:s_idx], s[s_idx + len(lus):]
    t_left, t_right = t[:t_idx], t[t_idx + len(lus):]

    return s_left, s_right, t_left, t_right


def get_search_space_right(s_sub, t_sub, k1, k2, common_alignment):
    """
    Purpose:
        Helper function: identifies all possible tandem duplicatations and deletions at right of lus
    Args:
        s_sub: subsequence of s
        t_sub: subsequence of t
        k1: int (such that k1 < k2 -- see tandem_transform_lcs)
        k2: int (longest possible kmer length)
        common_alignment: the lus
        
    Returns: 
        search_space: list of tuples informing specifc kmer transformation
            (kmer, "dup" or "del")
    """
    temp = []
    for i in range(k1, k2+1):
        temp.append(common_alignment[len(common_alignment) - i:])

    search_space = []
    for kmer in temp:
        # tandem deletion
        kmer = kmer.replace(" ", "")

        if s_sub[:len(kmer)] == kmer and t_sub[:len(kmer)] != kmer:
            search_space.append((kmer, 'del'))

        # tandem duplication
        if s_sub[:len(kmer)] != kmer and t_sub[:len(kmer)] == kmer:
            search_space.append((kmer, 'dup'))

    return search_space

def get_search_space_left(s_sub, t_sub, k1, k2, common_alignment):
    """
    Purpose:
        Helper function: identifies all possible tandem duplicatations and deletions at left of lus
    Args:
        s_sub: subsequence of s
        t_sub: subsequence of t
        k1: int (such that k1 < k2 -- see tandem_transform_lcs)
        k2: int (longest possible kmer length)
        common_alignment: the lus
        
    Returns: 
        search_space: list of tuples informing specifc kmer transformation
            (kmer, "dup" or "del")
    """

    temp = []
    for i in range(k1, k2+1):
        temp.append(common_alignment[:i])
        
    search_space = []
    for kmer in temp:
        # tandem deletion
        kmer = kmer.replace(" ", "")
        if s_sub[len(s_sub) - len(kmer):] == kmer and t_sub[len(t_sub) - len(kmer):] != kmer:
            search_space.append((kmer, 'del'))

        # tandem duplication
        if s_sub[len(s_sub) - len(kmer):] != kmer and t_sub[len(t_sub) - len(kmer):] == kmer:
            search_space.append((kmer, 'dup'))

    return search_space
        