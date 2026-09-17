seq1 = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"
seq2 = "MVHWTAEEKQLITGLWGKVNVAECGAEALARLLIVYPWTQRFFASFGNLSSPTAILGNPMVRAHGKKVLTSFGDAVKNLDNIKNTFSQLSELHCDKLHVDPENFRLLGDILIIVLAAHFSKDFTPECQAAWQKLVRVVAHALARKYH"

similar_groups = [
    ["K", "R", "H"],
    ["D", "E"],
    ["I", "L", "V", "M"],
    ["F", "Y", "W"],
    ["S", "T"],
    ["A", "G"] ]

matches = 0
similars = 0
gaps = 0
total_len = len(seq1)

for i in range(total_len):
    ch1 = seq1[i]
    ch2 = seq2[i]
    
    if ch1 == "-" or ch2 == "-":
        gaps = gaps + 1
    elif ch1 == ch2:
        matches = matches + 1
        similars = similars + 1
    else:
        for g in similar_groups:
            if ch1 in g and ch2 in g:
                similars = similars + 1
                break

original_query_len = 147
query_aligned_bases = 0
for ch in seq1:
    if ch != "-":
        query_aligned_bases = query_aligned_bases + 1

ident_pct = (matches / total_len) * 100
sim_pct = (similars / total_len) * 100
gap_pct = (gaps / total_len) * 100
cov_pct = (query_aligned_bases / original_query_len) * 100

print("Length:", total_len)
print("Matches:", matches)
print("Similars:", similars)
print("Gaps:", gaps)
print("Identity %:", ident_pct)
print("Similarity %:", sim_pct)
print("Gap %:", gap_pct)
print("Query Coverage %:", cov_pct)