DNA = "ATGGATTCCTTATCTTTACCCACCGTGAGTGACGTTACTGTTCCGGCTTTGAACGTTGACAGTCTCGTGCGAGACTATGTCAGCAATGTGAGAGCCGATGATAGCAATAACGTCAGTAGATTCCTCGGTGAAGTGGCCTTGAGAGAAATTAAATCTCAAGTTGACACCAGCAATGGTGATTTCCAGAAGTTAAACGTCGGTTTTCGTTTGACTCCTAATGAGAAGAACGCTTTGAAGGCGAATTTTCCTGGTCTTGAAATTGCGTTTAAAGACTCGTGTCATTCGTCTCATAGTTTTGCCGCACCACATAGAGTCTGCGAGACCCTAGAGATATACAATCGGTTTAAAACAAAAACCGATCGCATTATCGACCTTGGTGGCAATTATGTCACCCATGCTAAACAAGGTCGATCCAATGTGCATTCCTGCTGTCCCATCCTTGATGTGCGAGACGATGCTAGGCATACCGATCGTTATATATCTCTAGCTGCCTCTGTTGAGAATCGTCACAGAGAGTTAACAGTAGATTTCTGTTGCCATAAGTTCGAAGAATGTGATGTCAAAGCACCGTTTGCTATGGCAGTGCATTCCATCAGTGACATTCCTATTTCAACGGTTGCAACACATTGCGTTAGGCGCGGTGTGCGGAAACTTATTGCTTCCGTCATGATGGATCCGCTTATGATGCCTTATGACAAAGGTCATATACTCTTACTCAACGTTGAATGGGAAAAGGACGACATTGAAGAAGGAAAAACCCTGATTCATTTCCATTTTGTTGACGCCCCCGGATTAAGCTATCCGCACGATTTTAATGTACTGTCCCAATACATGATCACGAATCAGGTTATTGTTAATAATACCTACTCTTGTAGGGTAGAGAGGACAGCCTGTTTATCAGGTGTATATATCGTGGAGATGACTCTGTCCATGACAGATGGTTGTTCTTTAGCTTATCTG"

# Define the length of the segment to analyze
segment_length = 20

# Initialize variables to keep track of the top 5 segments, their counts, and their melting temperatures
top_segments = []
top_counts = []
top_tm = []

# Iterate through the DNA sequence with a sliding window of 20 bp
for i in range(len(DNA) - segment_length + 1):  # Ensure the last segment is fully covered
    segment = DNA[i:i + segment_length]
    count = 0
    
    # Count the occurrences of "GT" and "TG" in the current segment
    for j in range(segment_length - 1):
        if segment[j:j + 2] == "GT" or segment[j:j + 2] == "TG":
            count += 1
    
    # Update the top segments and counts
    if len(top_segments) < 5 or count > min(top_counts):
        if len(top_segments) == 5:
            min_count_index = top_counts.index(min(top_counts))
            del top_segments[min_count_index]
            del top_counts[min_count_index]
            del top_tm[min_count_index]
        top_segments.append(segment)
        top_counts.append(count)
        # Calculate the melting temperature (Tm) for the current segment
        A_count = segment.count("A")
        T_count = segment.count("T")
        C_count = segment.count("C")
        G_count = segment.count("G")
        tm = 2 * (A_count + T_count) + 4 * (C_count + G_count)
        top_tm.append(tm)

# Print the top 5 segments, their counts, and their melting temperatures
print("Top 5 segments with highest GT/TG content and their melting temperatures:")
for segment, count, tm in zip(top_segments, top_counts, top_tm):
    print(f"Segment: {segment}, Count: {count}, Melting Temperature: {tm} °C")
