from typing import Iterable, Dict

def read_blosum62(path: str):
    return None


def sequence_to_kmers(sequence: str, k_mer_lenght: int):
    return None


def find_similar_kmers_for_kmer(
        kmer: str,
        kmer_similarity_threshold: int,
        all_possible_kmers: Iterable[str],
        blosum_dict: Dict[str, Dict[str, int]]
):
    return None


def find_similar_kmers_for_sequence(
        query: str,
        kmer_similarity_threshold: int,
        all_possible_kmers: Iterable[str],
        blosum_dict: Dict[str, Dict[str, int]]
):
    return None


def create_index_pairs(
        query: str,
        database: str,
        kmer_size: int,
        kmer_similarity_threshold: int,
        blosum_dict: Dict[str, Dict[str, int]]
):
    return None


def merge_single_hit_with_single_hit(
        hit1,
        hit2,
        max_distance
):
    return None


def merge_extended_hit_with_single_hit(
    extended_hit,
    single_hit,
    max_distance
):
    return None


def merge_two_extended_hits(
        extended_hit1,
        extended_hit2,
        max_distance
):
    return None


def merge_two_hits(
        extended_hit1,
        extended_hit2,
        max_distance
):
    return None


def create_extended_hits(
        query: str,
        database: str,
        kmer_size: int,
        kmer_similarity_threshold: int,
        blosum_dict: Dict[str, Dict[str, int]],
        max_distance: int
):
    return None