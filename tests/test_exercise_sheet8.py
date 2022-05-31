import pytest
import os
from helpers.helpers import (
    convert_blosum_txt_to_dict_correct,
    index_sequence_by_kmers_correct,
    find_similar_kmers_for_kmer_correct,
    find_similar_kmers_for_sequence_correct,
    create_index_pairs_correct
)
from exercise_sheet8 import (
    read_blosum62,
    sequence_to_kmers,
    find_similar_kmers_for_kmer,
    find_similar_kmers_for_sequence,
    create_index_pairs
)
import random


TESTFOLDER = os.path.dirname(__file__)
PARENTDIR = os.path.dirname(TESTFOLDER)
HELPERDIR = os.path.join(PARENTDIR, "helpers")


def blosum62():
    file = os.path.join(HELPERDIR, "BLOSUM62.txt")
    assert os.path.exists(file)
    return file


def test_exercise_3a(blosum62):
    expected = convert_blosum_txt_to_dict_correct(blosum62)
    actual = read_blosum62(blosum62)
    assert actual == expected


@pytest.mark.parametrize(
    "sequence,k_mer_length",
    [
        (random_protein_seq(), 5),

    ]
)
def test_exercise_3b(sequence, k_mer_length):
    expected = index_sequence_by_kmers_correct(sequence, k_mer_length)
    actual = sequence_to_kmers(sequence, k_mer_length)
    assert actual == expected


def test_exercise_3c1(kmer, kmer_similarity_threshold, all_possible_kmers):
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
    assert actual == expected


def test_exercise_3c2(query, kmer_similarity_threshold, all_possible_kmers):
    dict_blosum = convert_blosum_txt_to_dict_correct(blosum62())
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
    assert actual == expected


def test_exercise_3c3(query, database, kmer_size, kmer_similarity_threshold):
    dict_blosum = convert_blosum_txt_to_dict_correct(blosum62())
    expected = create_index_pairs_correct(
        query, database, kmer_size, kmer_similarity_threshold, dict_blosum
    )
    actual = create_index_pairs(
        query, database, kmer_size, kmer_similarity_threshold, dict_blosum
    )
    assert actual == expected








