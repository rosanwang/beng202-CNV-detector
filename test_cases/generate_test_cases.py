import random

def totalkmerlen(kmers):
    total_kmer_len = 0
    for mer in kmers:
        total_kmer_len = total_kmer_len + len(mer)
    return total_kmer_len

def generateSeq(length):
    seq = ""
    bps = ["A", "T", "C", "G"]
    for i in range(length):
        base = random.choice(bps)
        seq = seq + base
    return seq 

def makeTestCase(kmers, s_len, t_len):
    kmer_len = totalkmerlen(kmers)
    if s_len < kmer_len or t_len < kmer_len:
        raise Exception("total kmer length must be less than length of s or t!")
    
    adjacent = [True, False]
    s = ""
    t = ""
    counter = 0
    tandom_operations = ["delete", "insert"]
    while len(s) < s_len and len(t) < t_len:
        operation = random.choice(tandom_operations)
        kmer = random.choice(kmers)
        adj = random.choice(adjacent)

        s_ind = s.find(kmer)
        t_ind = t.find(kmer)
        print(kmer)
        print(operation)
        print(adj)
        
        if operation == "delete":
            if t_ind == -1 or not adj:
                # insert one copies of kmer in t 
                t = t + kmer 
                # insert two copy of kmer in s 
                s = s + kmer + kmer
            elif s_ind != -1 and t_ind != -1 : 
                # insert one copy of kmer in s 
                s = s[:s_ind + len(kmer) ]  + kmer + kmer + s[s_ind + len(kmer):]
                t = t[:t_ind + len(kmer) ]  + kmer + t[t_ind + len(kmer) :]
            else: 
                raise Exception("funny seeing you here")
            
        else: # operation = "insert"

            if t_ind == -1 or not adj:
                # insert two copies of kmer in t 
                t = t + kmer + kmer
                # insert one copy of kmer in s 
                s = s + kmer 
            elif s_ind != -1 and t_ind != -1: 
                # insert one copy of kmer in s 
                s = s[:s_ind + len(kmer) ]  + kmer + s[s_ind + len(kmer):]
                t = t[:t_ind + len(kmer) ]  + kmer + kmer + t[t_ind + len(kmer) :]
            else: 
                raise Exception("funny seeing you here")
       
        if counter % 10 != 0:
            # insert common sequence between s and t 
            seq_len = random.randint(5, 8)
            seq = generateSeq(seq_len)
            print("common seq:" + seq)
            s = s + seq
            t = t + seq

        counter = counter + 1


    return s, t

if __name__ == "__main__":
    # kmers = ["ATCAGAGTTA", "GCAGCTGCAGA", "AGATATTTACA"]
    kmers = ["GAT", "GGAGTCCGTTG", "TACGTGCGT"]
    # k1 = 3, k2 = 11 
    out = makeTestCase(kmers, 100, 100)
    print(out)

# TA
# delete
# GGA
# delete
# common seq:GGTCGGAT
# TA
# insert
# common seq:GCCAAAG

# 'TATAAGGAGGAGGTCGGATGCCAAAG'
# 'TATATAGAGGTCGGATGCCAAAG'