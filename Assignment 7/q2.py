import math

hydrophobocity_values = {"A": 13.85,
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
    "Y": 13.88,}

def cal_mu(seg, angle_deg):
    delta = math.radians(angle_deg)
    sin_sum = 0.0
    cos_sum = 0.0
    for n, aa in enumerate(seg, start=1):
        h = hydrophobocity_values[aa]
        sin_sum += h * math.sin(n * delta)
        cos_sum += h * math.cos(n * delta)
    return (math.sqrt(sin_sum**2 + cos_sum**2)) / len(seg)

# β-strands (length 6, delta = 160 deg)
strands = [
    ("seq1 strand 1 (9-14)", "TNIYGY"),
    ("seq1 strand 2 (29-34)", "NTGIVS"),
    ("seq1 strand 3 (73-78)", "VISLGF"),
    ("seq1 strand 4 (97-102)", "IKWYVD"),
    ("seq2 strand (11-16)", "GVLGAI"),
    ("seq3 strand (7-12)", "SIAAAT"),
]

# α-helices (length 8, delta = 100 deg)
helices = [
    ("seq2 helix 1 (10-17)", "AGVLGAIL"),
    ("seq2 helix 2 (44-51)", "IFISEAII"),
    ("seq2 helix 2b (48-55)", "SEAIIHVL"),
    ("seq2 helix 3 (68-75)", "NKALELFR"),
]

print("--- β-strands (N=6, angle=160 deg) ---")
for label, seq in strands:
    val = cal_mu(seq, 160)
    print(f"{label:26} : {seq} -> {val:.3f}")

print("\n--- α-helices (N=8, angle=100 deg) ---")
for label, seq in helices:
    val = cal_mu(seq, 100)
    print(f"{label:26} : {seq} -> {val:.3f}")