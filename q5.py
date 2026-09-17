seq1 = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"
seq2 = "MVHWTAEEKQLITGLWGKVNVAECGAEALARLLIVYPWTQRFFASFGNLSSPTAILGNPMVRAHGKKVLTSFGDAVKNLDNIKNTFAQLSELHCDKLHVDPENFRLLGDILIIVLAAHFSKDFTPECQAAWQKLVRVVAHALAHKYH"

dict1 = {}
dict2 = {}

for i in range(len(seq1) - 4):
    pp = seq1[i:i+5]
    if pp in dict1:
        dict1[pp] = dict1[pp] + 1
    else:
        dict1[pp] = 1

for i in range(len(seq2) - 4):
    pp = seq2[i:i+5]
    if pp in dict2:
        dict2[pp] = dict2[pp] + 1
    else:
        dict2[pp] = 1

ans = []
for k in dict1:
    if k in dict2:
        ans.append(k)

print("Matching pentapeptides:", len(ans))
print("Pentapeptides seq1_count seq2_count")

for x in ans:
    print(x, "\t\t", dict1[x], "\t\t", dict2[x])