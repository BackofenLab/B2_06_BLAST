# In this file you can find sorted helpers for the BLAST implementation

import random
random.seed(42)


# 1. Reading the blosum file

def random_protein_seq(length):
    """
    Generate a random protein sequence.
    """
    return ''.join(random.choice('ACDEFGHIKLMNPQRSTVWY') for _ in range(length))


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


def create_index_pairs_correct(query_sequence, database_sequence, kmer_size, kmer_similarity_threshold, dict_blosum):
    """
    Create an index of the query sequence and the database sequence.
    """
    dict_both_indexes = {}
    indexed_database_sequence = index_sequence_by_kmers_correct(database_sequence, kmer_size)
    indexed_query_sequence = index_sequence_by_kmers_correct(query_sequence, kmer_size)
    all_possible_kmers = list(indexed_database_sequence.keys()) + list(indexed_query_sequence.keys())
    dict_similar_kmers_query = find_similar_kmers_for_sequence_correct(query_sequence, kmer_similarity_threshold,
                                                                       all_possible_kmers, dict_blosum)
    for kmer, kmer_indexes_query in indexed_query_sequence.items():
        indexes_database_all_similar = []
        for similar_kmer in dict_similar_kmers_query[kmer]:
            indexes_similar_kmer_in_database = indexed_database_sequence.get(similar_kmer, [])
            indexes_database_all_similar.extend(indexes_similar_kmer_in_database)
        dict_both_indexes[kmer] = kmer_indexes_query, indexes_database_all_similar

    return dict_both_indexes


# 4. Merging hits

def merge_single_hit_with_single_hit_correct(single_hit1, singe_hit2, max_distance):
    """
    Merge two extended hits.
    """
    dist_first = abs(single_hit1[0] - singe_hit2[0]) - 1
    dist_second = abs(single_hit1[1] - singe_hit2[1]) - 1
    if dist_first + dist_second <= max_distance:
        return True, sorted([single_hit1, singe_hit2])
    return False, None


def merge_multiple_hits_with_single_hit_correct(multiple_hits, single_hit, max_distance):
    """
    Merge multiple hits with a single hit.
    """
    for hit in multiple_hits:
        flag_can_be_merged, new_multiple_hits = merge_single_hit_with_single_hit_correct(hit, single_hit, max_distance)
        if flag_can_be_merged:
            multiple_hits.append(single_hit)
            return True, sorted(multiple_hits)
    return False, None


def merge_two_extended_hits(extended_hit1, extended_hit2, max_distance):
    """
    Merge two extended hits.
    """
    if (type(extended_hit1) == tuple) and (type(extended_hit2) == tuple):
        can_be_merged, merged = merge_single_hit_with_single_hit_correct(extended_hit1, extended_hit2, max_distance)
        return can_be_merged, merged
    elif type(extended_hit1) == list and type(extended_hit2) == tuple:
        can_be_merged, merged = merge_multiple_hits_with_single_hit_correct(extended_hit1, extended_hit2, max_distance)
        return can_be_merged, merged
    elif type(extended_hit1) == tuple and type(extended_hit2) == list:
        can_be_merged, merged = merge_multiple_hits_with_single_hit_correct(extended_hit2, extended_hit1, max_distance)
        return can_be_merged, merged
    else:
        for single_hit in extended_hit1:
            can_be_merged, merged = merge_single_hit_with_single_hit_correct(extended_hit2, single_hit, max_distance)
            if can_be_merged:
                merged = extended_hit1 + extended_hit2
                return True, sorted(merged)
    return False, None


def merging_one_iteration(list_extended_hits, max_distance):
    """
    Merge hits in one iteration.
    """
    for i in range(len(list_extended_hits) - 1):
        for j in range(i + 1, len(list_extended_hits)):
            can_be_merged, merged = merge_two_extended_hits(list_extended_hits[i], list_extended_hits[j], max_distance)
            if can_be_merged:
                first_element_to_remove = list_extended_hits[i]
                second_element_to_remove = list_extended_hits[j]
                list_extended_hits.remove(first_element_to_remove)
                list_extended_hits.remove(second_element_to_remove)
                list_extended_hits.append(merged)
                return False, list_extended_hits
    return True, list_extended_hits


def create_extended_hits_correct(query_sequence, database_sequence, kmer_size,
                                 kmer_similarity_threshold, dict_blosum, max_distance):
    dict_index_pairs = create_index_pairs_correct(query_sequence, database_sequence,
                                                  kmer_size, kmer_similarity_threshold, dict_blosum)

    flat_list_hits = []
    for kmer, (kmer_indexes_query, indexes_database_all_similar) in dict_index_pairs.items():
        all_pairs = [(i, j) for i in kmer_indexes_query for j in indexes_database_all_similar]
        flat_list_hits.extend(all_pairs)

    print(flat_list_hits)

    list_extended_hits = flat_list_hits
    flag_all_potential_merges_are_done = False
    while not flag_all_potential_merges_are_done:
        flag_all_potential_merges_are_done, list_extended_hits = merging_one_iteration(list_extended_hits, max_distance)

    return list_extended_hits


def main():

    ###############  1 Start ################
    dict_blosum = convert_blosum_txt_to_dict_correct("BLOSUM62.txt")
    print(dict_blosum)

    ###############  1 End ##################

    ###############  2 Start ################

    database = "DPPEGVVDPP"
    query = "RPPQGLF"

    indexes_db = index_sequence_by_kmers_correct(database, 3)
    indexes_query = index_sequence_by_kmers_correct(query, 3)

    print(indexes_db)
    print(indexes_query)

    ###############  2 End ##################


    ###############  3 Start ################

    all_existing_kmers = list(indexes_db.keys()) + list(indexes_query.keys())
    all_similar_kmers = find_similar_kmers_for_kmer_correct("PQG", 13, all_existing_kmers, dict_blosum)
    print(all_similar_kmers)

    all_similar_kmers_for_sequence = find_similar_kmers_for_sequence_correct("DPPEGVVDPP", 13, all_existing_kmers, dict_blosum)
    print(all_similar_kmers_for_sequence)

    dict_both_indexes = create_index_pairs_correct(query, database, 3, 5, dict_blosum)
    print(dict_both_indexes)

    ###############  3 End ##################


    ###############  4 Start ################


    extension, extedned_hit = merge_single_hit_with_single_hit_correct((1, 1), (2, 8), 2)
    print(extension, extedned_hit)



    ###############  4 End ################





    extended_hits = create_extended_hits_correct(query, database, 3, 12, dict_blosum, 2)
    print(extended_hits)


if __name__ == "__main__":
    main()