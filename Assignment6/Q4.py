import math

AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")

def read_blosum62(filename="BLOSUM62"):
    mtx = {}
    headers = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if not headers and any(aa in parts for aa in AMINO_ACIDS):
                headers = parts
                continue
            if headers and parts:
                row_aa = parts[0]
                for col_aa, val in zip(headers, parts[1:]):
                    try:
                        v = float(val)
                        mtx[(row_aa, col_aa)] = v
                        mtx[(col_aa, row_aa)] = v
                    except ValueError:
                        continue
    return mtx

def read_clustal_aln(filepath):
    seqs = {}
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('*') or line.startswith(':') or line.startswith('.'):
                continue
            tokens = line.split()
            if len(tokens) >= 2:
                seq_id = tokens[0]
                seq_chunk = tokens[1]
                if all(c in "ACDEFGHIKLMNPQRSTVWY-XBZ" for c in seq_chunk.upper()):
                    seqs[seq_id] = seqs.get(seq_id, "") + seq_chunk
    return list(seqs.values())

def calculate_conservation_measures(aln_sequences, blosum_matrix):
    num_positions = len(aln_sequences[0])
    results = []

    for pos in range(num_positions):
        col = [seq[pos].upper() for seq in aln_sequences if seq[pos] != '-']
        N = len(col)
        
        if N == 0:
            results.append((pos + 1, 0.0, 0.0, 0.0))
            continue

        counts = {aa: col.count(aa) for aa in AMINO_ACIDS}
        freqs = {aa: counts[aa] / N for aa in AMINO_ACIDS}

        entropy = sum(f * math.log(f) for f in freqs.values() if f > 0)

        variance = math.sqrt(sum((f - 0.05) ** 2 for f in freqs.values()))

        sop_score = 0.0
        if N > 1:
            total_pair_val = 0.0
            total_pairs = 0
            for i in range(N):
                for j in range(i + 1, N):
                    a1, a2 = col[i], col[j]
                    score = blosum_matrix.get((a1, a2), blosum_matrix.get((a2, a1), 0.0))
                    total_pair_val += score
                    total_pairs += 1
            sop_score = total_pair_val / total_pairs if total_pairs > 0 else 0.0

        results.append((pos + 1, entropy, variance, sop_score))

    return results

if __name__ == "__main__":
    mtx = read_blosum62("BLOSUM62")
    alignment = read_clustal_aln("set1.aln")
    scores = calculate_conservation_measures(alignment, mtx)

    print(f"{'Position':<10} {'Entropy':<12} {'Variance':<12} {'Sum-of-Pairs':<12}")
    print("-" * 48)
    for pos, ent, var, sop in scores[:25]:
        print(f"{pos:<10} {ent:<12.4f} {var:<12.4f} {sop:<12.4f}")