def read_entropy(file_path):

    scores = []
    with open(file_path, 'r') as f:
        lines = f.readlines()
        
    for line in lines[1:]:
        parts = line.strip().split(',')
        if len(parts) >= 3:
            scores.append(float(parts[2].strip()))
    return scores

def compare_tools(clustal_csv, mafft_csv, muscle_csv):
   
    clustal_scores = read_entropy(clustal_csv)
    mafft_scores = read_entropy(mafft_csv)
    muscle_scores = read_entropy(muscle_csv)
  

    # Use the minimum length to prevent IndexError across different algorithms
    min_length = min(len(clustal_scores), len(mafft_scores), len(muscle_scores))
    
    # Print the table header
    print(f"{'Position':<10} {'Clustal':<10} {'MAFFT':<10} {'MUSCLE':<10}")
    
    # Initialize the tracking lists
    similar = []
    different = []

    for i in range(min_length):
        c = clustal_scores[i]
        m = mafft_scores[i]
        mu = muscle_scores[i]
        
        # Print the formatted row values
        print(f"{i+1:<10} {c:<10.4f} {m:<10.4f} {mu:<10.4f}")
        
        # Strict restriction replaced with < 0.0001 tolerance check
        if abs(c - m) < 0.0001 and abs(c - mu) < 0.0001:
            similar.append(i + 1)
        else:
            different.append(i + 1)


    print("\nSimilar conservation scores:")
    print(similar)
    
    print("\nDifferent conservation scores:")
    print(different)


compare_tools(
    "set2_scores_clustal.csv", 
    "set2_scores_mafft.csv", 
    "set2_scores_muscle.csv"
)