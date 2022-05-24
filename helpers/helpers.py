# In this file you can find sorted helpers for the BLAST implementation


# 1. Reading the blossum file

def convert_blosum_txt_to_dict_correct(filepath):
    """
    Convert a BLOSUM file to a dictionary.
    """
    with open(filepath, 'r') as f:
        lines = f.readlines()
    list_amino_acids = lines[6].strip().split()
    dict_blosum = {}
    for line in lines[7:]:
        line = line.strip().split()
        dict_blosum[line[0]] = {amino_acid: int(score) for amino_acid, score in zip(list_amino_acids, line[1:])}
    return dict_blosum


# 2. Indexing a sequence

def index_sequence_by_kmers_correct(sequence, k):
    """
    Index a sequence by k-mers.
    """
    dict_kmer_positions = {}
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        if kmer not in dict_kmer_positions:
            dict_kmer_positions[kmer] = []
        dict_kmer_positions[kmer].append(i+1)

    return dict_kmer_positions


# 3. Searching for all similar kmers with a given threshold


def generate_all_possible_kmers_of_size_k_correct(k, dict_blosum):
    """
    Generate all possible kmers of size k.
    """
    all_possible_characters = list(dict_blosum.keys())
    all_possible_kmers = [""]
    for _ in range(k):
        all_possible_kmers = [word_ + character for word_ in all_possible_kmers
                              for character in all_possible_characters]
    return all_possible_kmers


def find_similar_kmers_for_kmer_correct(kmer, kmer_similarity_threshold, all_possible_kmers, dict_blosum):
    """
    Find similar words.
    """
    kmers_within_threshold = []
    for word_ in all_possible_kmers:
        difference = sum(dict_blosum[word_[i]][kmer[i]] for i in range(len(kmer)))
        if difference >= kmer_similarity_threshold:
            kmers_within_threshold.append(word_)
    return kmers_within_threshold


def find_similar_kmers_for_sequence_correct(sequence, kmer_similarity_threshold, all_possible_kmers, dict_blosum):
    """
    Find similar words.
    """
    dict_similar_kmers = {}
    kmer_size = len(all_possible_kmers[0])
    for i in range(len(sequence) - kmer_size + 1):
        kmer = sequence[i:i+kmer_size]
        if kmer not in dict_similar_kmers:
            kmers_within_threshold = find_similar_kmers_for_kmer_correct(kmer, kmer_similarity_threshold,
                                                                 all_possible_kmers, dict_blosum)
            dict_similar_kmers[kmer] = kmers_within_threshold
    return dict_similar_kmers


def create_index_pairs(query_sequence, database_sequence, kmer_size, kmer_similarity_threshold, dict_blosum):
    """
    Create an index of the query sequence and the database sequence.
    """
    dict_both_indexes = {}
    indexed_database_sequence = index_sequence_by_kmers_correct(database_sequence, kmer_size)
    indexed_query_sequence = index_sequence_by_kmers_correct(query_sequence, kmer_size)
    all_possible_kmers = generate_all_possible_kmers_of_size_k_correct(kmer_size, dict_blosum)
    dict_similar_kmers_query = find_similar_kmers_for_sequence_correct(query_sequence, kmer_similarity_threshold,
                                                               all_possible_kmers, dict_blosum)
    for kmer, kmer_indexes_query in indexed_query_sequence.items():
        indexes_database_all_similar = []
        for similar_kmer in dict_similar_kmers_query[kmer]:
            indexes_similar_kmer_in_database = indexed_database_sequence.get(similar_kmer, [])
            indexes_database_all_similar.extend(indexes_similar_kmer_in_database)
        dict_both_indexes[kmer] = kmer_indexes_query, indexes_database_all_similar

    return dict_both_indexes




def main():
    dict_blosum = convert_blosum_txt_to_dict_correct("BLOSUM62.txt")
    print(dict_blosum)

    database = "DPPEGVVDPP"
    query = "RPPQGLF"

    dict_both_indexes = create_index_pairs(query, database, 3, 12, dict_blosum)
    print(dict_both_indexes)



if __name__ == "__main__":
    main()