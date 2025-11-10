#!/Users/wuchtys/opt/miniconda3/bin/python3
import sys
from bm_preproc import *


def boyer_moore_with_counts(p, p_bm, t):
    m = len(p)
    n = len(t)
    i = 0
    occurrences = []
    n_align = 0
    n_comp = 0

    while i <= n - m:
        n_align += 1
        shift = 1
        mismatched = False

        # scan right-to-left
        for j in range(m - 1, -1, -1):
            n_comp += 1
            if p[j] != t[i + j]:
                # bad-character and good-suffix rules
                bc = p_bm.bad_character_rule(j, t[i + j])
                gs = p_bm.good_suffix_rule(j)
                if bc > shift:
                    shift = bc
                if gs > shift:
                    shift = gs
                mismatched = True
                break

        if not mismatched:
            occurrences.append(i)
            ms = p_bm.match_skip()  # standard BM shift after a full match
            if ms > shift:
                shift = ms

        i += shift

    return occurrences, n_align, n_comp


# def boyer_moore (p, bm_p, t):   
#     i = 0                                    #running index
#     occurrences = []
    
#     while i < len(t) - len(p) + 1:
#         shift = 1
#         mismatch = False
#         for j in range(len(p) - 1, -1, -1):  #running through p from
#             if not p[j] == t[i+j]:           #the end of p 
#                 skip_bc = bm_p.bad_character_rule(j, t[i+j]) 
#                 skip_gs = bm_p.good_suffix_rule(j)
#                 shift = max(shift, skip_bc, skip_gs) #which gives
#                 mismatch = True                      #biggest shift
#                 break
#         if not mismatch:                   #match! 
#             occurrences.append(i)
#             skip_ms = bm_p.match_skip()
#             shift = max(shift, skip_ms)     
#         i += shift
#     return occurrences


############################ Tests performed ##############################

# Test 1
def read_fasta(path):
    seq = []
    with open(path) as fh:
        for line in fh:
            if not line or line[0] == '>':
                continue
            seq.append(line.strip())
    return "".join(seq)


t = read_fasta("chr1.GRCh38.excerpt.fasta")
p = "GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG"
p_bm = BoyerMoore(p)
occ, n_algn, n_char = boyer_moore_with_counts(p, p_bm, t)
print(occ, n_algn, n_char)
## Test1 output -> [56922] 127974 165191

# Test 2
p = "word"
t = "there would have been a time for such a word"
lowercase_alphabet = "abcdefghijklmnopqrstuvwxyz "
p_bm = BoyerMoore(p, lowercase_alphabet)
occurrences, n_algn, n_char = boyer_moore_with_counts(p, p_bm, t)
print(occurrences, n_algn, n_char)  
## Test2 output -> [40] 12 15

# Test 3
t = 'GCTACGATCTAGAATCTA'
p = 'TCTA'
p_bm = BoyerMoore(p)
print (boyer_moore_with_counts(p, p_bm, t))
## Test2 output -> ([7, 14], 7, 16)

