import matplotlib.pyplot as plt

hydro_scale = {
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

# Sequences given in the fasta file
s1 = "FDCAEYRSTNIYGYGLYEVSMKPAKNTGIVSSFFTYTGPAHGTQWEIDIEFLGKDTTKVQFNYYTNGVGGHEKVISLGFDASKGFHTYAFDWQPGYIKWYVDGVLK"
s2 = "KASEDLVKKHAGVLGAILKKKGHHEAELKPLAQSHATKAHKNIFISEAIIHVLHSRHPGDFGADAQGAMNKALELFRKDIAAKYKELGY"
s3 = "TVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA"

sequences = [("seq1", s1), ("seq2", s2), ("seq3", s3)]


# Sliding window
def profile(seq, win):
    scores = []
    pos = []
    half = win // 2
    for i in range(len(seq) - win + 1):
        window_seq = seq[i : i + win]
        total = 0
        for aa in window_seq:
            total = total + hydro_scale[aa]
        avg = total / win
        scores.append(avg)
        pos.append(i + half + 1)  
    return pos, scores


# Secondary structure scanning:
# Window 5-7 for beta sheets, 9-11 or 19 for alpha helices
# Using window 7 for beta and 11 for alpha
for name, s in sequences:
    print("=" * 40)
    print("Sequence:", name)
    print("Length:", len(s))

    pos_beta, prof_beta = profile(s, 7)
    pos_alpha, prof_alpha = profile(s, 11)

    # find regions above average hydrophobicity threshold (~13.6)
    threshold = 13.6

    print("Predicted Beta-strand candidate centers (window=7, cutoff > 13.6):")
    betapeaks = []
    for p, score in zip(pos_beta, prof_beta):
        if score > threshold:
            betapeaks.append((p, round(score, 2)))
    print(betapeaks)

    print("Predicted Alpha-helix candidate centers (window=11, cutoff > 13.6):")
    alpha_peaks = []
    for p, score in zip(pos_alpha, prof_alpha):
        if score > threshold:
            alpha_peaks.append((p, round(score, 2)))
    print(alpha_peaks)

    # plot
    plt.figure(figsize=(10, 4))
    plt.plot(pos_beta, prof_beta, label="Beta window (w=7)", color="blue")
    plt.plot(pos_alpha, prof_alpha, label="Alpha window (w=11)", color="red")
    plt.axhline(
        y=threshold, color="gray", linestyle="--", label="Hydrophobic cutoff"
    )
    plt.title(f"Hydrophobicity Profile - {name}")
    plt.xlabel("Residue Position")
    plt.ylabel("Mean Hydrophobicity")
    plt.legend()
    plt.tight_layout()
    plt.show()