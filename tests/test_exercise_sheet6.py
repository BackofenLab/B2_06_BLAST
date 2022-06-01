import pytest
import os
from helpers.helpers import (
    convert_blosum_txt_to_dict_correct,
    index_sequence_by_kmers_correct,
    find_similar_kmers_for_kmer_correct,
    find_similar_kmers_for_sequence_correct,
    create_index_pairs_correct,
    random_protein_seq,
    merge_single_hit_with_single_hit_correct,
    merge_extended_hit_with_single_hit_correct,
    merge_two_extended_hits_correct,
    merge_two_hits_correct,
    create_extended_hits_correct
)

from exercise_sheet6 import (
    read_blosum62,
    index_sequence_by_kmers,
    find_similar_kmers_for_kmer,
    find_similar_kmers_for_sequence,
    create_index_pairs,
    merge_single_hit_with_single_hit,
    merge_extended_hit_with_single_hit,
    merge_two_extended_hits,
    merge_two_hits,
    create_extended_hits
)
import random

random.seed(41)


TESTFOLDER = os.path.dirname(__file__)
PARENTDIR = os.path.dirname(TESTFOLDER)
HELPERDIR = os.path.join(PARENTDIR, "helpers")


def blosum62():
    file = os.path.join(HELPERDIR, "BLOSUM62.txt")
    assert os.path.exists(file)
    return file


def test_exercise_3a():
    blosum_dict = blosum62()
    expected = convert_blosum_txt_to_dict_correct(blosum_dict)
    actual = read_blosum62(blosum_dict)
    assert actual == expected


@pytest.mark.parametrize(
    "sequence,k_mer_length",
    [
        (random_protein_seq(10), 3),
        (random_protein_seq(100), 5),
        (random_protein_seq(100), 4),

    ]
)
def test_exercise_3b(sequence, k_mer_length):
    expected = index_sequence_by_kmers_correct(sequence, k_mer_length)
    actual = index_sequence_by_kmers(sequence, k_mer_length)
    for key,value in expected.items():
        expected[key] = sorted(value)
    if actual is not None:
        for key, value in actual.items():
            actual[key] = sorted(value)
    assert actual == expected


def random_mutated_subseq(sequence, length, mutations):
    alphabet = list(convert_blosum_txt_to_dict_correct(blosum62()))
    index = random.randint(0, len(sequence) - length)
    subseq = list(sequence[index:index+length])
    for _ in range(mutations):
        idx = random.randrange(0, len(subseq))
        a_idx = random.randrange(0, len(alphabet))
        letter = alphabet[a_idx]
        subseq[idx] = letter
    return "".join(subseq)







@pytest.mark.parametrize(
    "kmer, kmer_similarity_threshold",
    [
        (random_protein_seq(7), 5),
        (random_protein_seq(5), 5),
        (random_protein_seq(10), 5),
        (random_protein_seq(3), 3),

    ]
)
def test_exercise_3c1(kmer, kmer_similarity_threshold):
    all_possible_kmers = set(random_mutated_subseq(kmer, len(kmer), 3) for _ in range(10))
    dict_blosum = convert_blosum_txt_to_dict_correct(blosum62())
    expected = find_similar_kmers_for_kmer_correct(
        kmer,
        kmer_similarity_threshold,
        all_possible_kmers,
        dict_blosum
    )
    actual = find_similar_kmers_for_kmer(
        kmer,
        kmer_similarity_threshold,
        all_possible_kmers,
        dict_blosum
    )
    assert sorted(actual) == sorted(expected)


@pytest.mark.parametrize(
    "query, kmer_similarity_threshold, kmer_length",
    [
        (random_protein_seq(20), 5, 7),
        (random_protein_seq(30), 5, 5),
        (random_protein_seq(10), 5, 5),
        (random_protein_seq(30), 10, 3),

    ]
)
def test_exercise_3c2(query, kmer_similarity_threshold, kmer_length):
    dict_blosum = convert_blosum_txt_to_dict_correct(blosum62())
    all_possible_kmers = list(set(random_mutated_subseq(query, kmer_length, 3) for _ in range(10)))
    expected = find_similar_kmers_for_sequence_correct(
        query,
        kmer_similarity_threshold,
        all_possible_kmers,
        dict_blosum
    )
    actual = find_similar_kmers_for_sequence(
        query,
        kmer_similarity_threshold,
        all_possible_kmers,
        dict_blosum
    )
    for key, value in expected.items():
        expected[key] = sorted(value)
    if actual is not None:
        for key, value in actual.items():
            actual[key] = sorted(value)
    assert actual == expected


@pytest.mark.parametrize(
    "database, query_size, kmer_size, kmer_similarity_threshold, mutations",
    [
        (random_protein_seq(100), 20, 7, 10, 20),
        (random_protein_seq(100), 10, 5, 20, 5),
        (random_protein_seq(200), 50, 10, 15, 10),
        (random_protein_seq(30), 10, 3, 20, 10),
        (random_protein_seq(100), 20, 7, 10, 100),
    ]
)
def test_exercise_3c3(database, query_size, kmer_size, kmer_similarity_threshold, mutations):
    dict_blosum = convert_blosum_txt_to_dict_correct(blosum62())
    query = random_mutated_subseq(database, query_size, mutations)
    expected = create_index_pairs_correct(
        query, database, kmer_size, kmer_similarity_threshold, dict_blosum
    )
    actual = create_index_pairs(
        query, database, kmer_size, kmer_similarity_threshold, dict_blosum
    )
    for key,value in expected.items():
        expected[key] = sorted(value[0]), sorted(value[1])
    if actual is not None:
        for key, value in actual.items():
            actual[key] = sorted(value[0]), sorted(value[1])
    assert actual == expected


@pytest.mark.parametrize(
    "hit1, hit2, max_distance",
    [
        ((1, 1), (2, 2), 2),
        ((1, 1), (2, 8), 2),
        ((1, 4), (3, 8), 3)
    ]
)
def test_exercise_3d1(hit1, hit2, max_distance):
    expected = merge_single_hit_with_single_hit_correct(
        hit1,
        hit2,
        max_distance
    )
    actual = merge_single_hit_with_single_hit(
        hit1,
        hit2,
        max_distance
    )
    assert actual == expected


@pytest.mark.parametrize(
    "single_hit, extended_hit, max_distance",
    [
        ((1, 1), [(2, 2), (3, 3)], 2),
        ((7, 8), [(2, 2), (3, 3)], 2),
    ]
)
def test_exercise_3d2(single_hit, extended_hit, max_distance):
    expected = merge_extended_hit_with_single_hit_correct(extended_hit[:], single_hit, max_distance)
    actual = merge_extended_hit_with_single_hit(extended_hit[:], single_hit, max_distance)
    assert actual == expected


@pytest.mark.parametrize(
    "extended_hit1, extended_hit2, max_distance",
    [
        ([(4, 4), (5, 5)], [(2, 2), (3, 3)], 2),
        ([(7, 8), (9, 9)], [(2, 2), (3, 3)], 2),
    ]
)
def test_exercise_3d3(extended_hit1, extended_hit2, max_distance):
    expected = merge_two_extended_hits_correct(extended_hit1[:], extended_hit2[:], max_distance)
    actual = merge_two_extended_hits(extended_hit1[:], extended_hit2[:], max_distance)
    assert actual == expected


@pytest.mark.parametrize(
    "hit1, hit2, max_distance",
    [
        ([(4, 4), (5, 5)], [(2, 2), (3, 3)], 2),
        ([(7, 8), (9, 9)], [(2, 2), (3, 3)], 2),
        ((1, 1), [(2, 2), (3, 3)], 2),
        ((7, 8), [(2, 2), (3, 3)], 2),
        ((1, 1), (2, 2), 2),
        ((1, 1), (2, 8), 2),
        ((1, 4), (3, 8), 3)
    ]
)
def test_exercise_3d4(hit1, hit2, max_distance):
    expected = merge_two_hits_correct(hit1, hit2, max_distance)
    actual = merge_two_hits(hit1, hit2, max_distance)
    assert actual == expected



@pytest.mark.parametrize(
    "query_size, database, kmer_size, kmer_similarity_threshold, max_distance, mutations",
    [
        (20, random_protein_seq(100), 3, 10, 2, 20),
        (10, random_protein_seq(100), 5, 15, 4, 5),
        (15, random_protein_seq(200), 5, 50, 3, 10),
        (30, random_protein_seq(300), 7, 10, 2, 10),
        (10, random_protein_seq(100), 5, 20, 2, 100),
    ]
)
def test_exercise_3e(query_size, database, kmer_size, kmer_similarity_threshold, max_distance, mutations):
    blosum_dict = convert_blosum_txt_to_dict_correct(blosum62())
    query = random_mutated_subseq(database, query_size, mutations)
    expected = create_extended_hits_correct(query, database, kmer_size, kmer_similarity_threshold, blosum_dict, max_distance)
    actual = create_extended_hits(query, database, kmer_size, kmer_similarity_threshold, blosum_dict, max_distance)
    if actual is None:
        assert False
    for entry in expected:
        assert entry in actual











