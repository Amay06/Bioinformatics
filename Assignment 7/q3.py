import matplotlib.pyplot as plt

hydrophobicity_scale = {
    "A": 13.85,
    "D": 11.61,
    "C": 15.37,
    "E": 11.38,
    "F": 13.93,
    "G": 13.34,
    "H": 13.82,
    "I": 15.28,
    "K": 11.58,
    "L": 14.13,
    "M": 13.86,
    "N": 13.02,
    "P": 12.35,
    "Q": 12.61,
    "R": 13.10,
    "S": 13.39,
    "T": 12.70,
    "V": 14.56,
    "W": 15.48,
    "Y": 13.88,
}

seq = "ALLSFERKYRVRGGTLIGGDLFDFWVGPYFVGFFGVSAIFFIFLGVSLIGYAASQGPTWDPFAISINPPDLKYGLGAAPLLEGGFWQAITVCALGAFISWMLREVEISRKLGIGWHVPLAFCVPIFMFCVLQVFRPLLLGSWGHAFPYGILSHLDWVNNFGYQYLNWHYNPGHMSSVSFLFVNAMALGLHGGLILSVANPGDGDKVKTAEHENQYFRDVVGYSIGALSIHRLGLFLASNIFLTGAFGTIASGPFWTRGWPEWWGWWLDIPFWS"

def get_profile(s, win):
    scores = []
    pos = []
    half = win // 2
    for i in range(len(s) - win + 1):
        window_seq = s[i : i + win]
        tot = sum(hydrophobicity_scale[aa] for aa in window_seq)
        scores.append(tot / win)
        pos.append(i + half + 1)
    return pos, scores


pos9, prof9 = get_profile(seq, 9)
pos19, prof19 = get_profile(seq, 19)

# plot profiles
plt.figure(figsize=(12, 5))
plt.plot(pos9, prof9, label="Window = 9", color="blue", alpha=0.7)
plt.plot(pos19, prof19, label="Window = 19 (TM)", color="red", linewidth=2)
plt.axhline(
    y=13.8, color="gray", linestyle="--", label="Hydrophobic threshold (~13.8)"
)
plt.title("Hydrophobicity Profile - 1PRC:L")
plt.xlabel("Residue Position")
plt.ylabel("Mean Hydrophobicity")
plt.legend()
plt.tight_layout()
plt.show()

# scan continuous TM regions where window 19 stays above ~13.8
threshold = 13.8
print("Residues where Window 19 > 13.8:")
high_tm = [p for p, sc in zip(pos19, prof19) if sc >= threshold]

# group consecutive points into contiguous blocks
blocks = []
current = []
for p in high_tm:
    if not current or p == current[-1] + 1:
        current.append(p)
    else:
        blocks.append((current[0] - 9, current[-1] + 9))
        current = [p]
if current:
    blocks.append((current[0] - 9, current[-1] + 9))

print("\nCandidate Transmembrane Segments (approximate residue bounds):")
for i, (start, end) in enumerate(blocks, start=1):
    print(f"TM {i}: residues {start} to {end} -> {seq[start-1:end]}")