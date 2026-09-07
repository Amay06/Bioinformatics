#Code for Q2
seq1 = "AATCTATA"
seq2 = "AAG--ATA"

match_score = 1
mismatch_score = 0
origination_penalty = -2
length_penalty = -1

score = 0
in_gap = False

for a, b in zip(seq1, seq2):
    if a == '-' or b == '-':
        if not in_gap:
            score += origination_penalty   
            in_gap = True
        score += length_penalty            
    else:
        in_gap = False
        if a == b:
            score += match_score
        else:
            score += mismatch_score

print("Total alignment score:", score)