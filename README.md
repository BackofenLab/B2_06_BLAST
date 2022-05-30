<img src="./figures/banner.png" alt="UniFreiburg Banner"/>

Lehrstuhl für Bioinformatik - Institut für Informatik - *http://www.bioinf.uni-freiburg.de*

---
## Bioinformatics 2
###### SS 2022
##### Exercise sheet 6: BLAST
---

### _Exercise 3 - Programming Assignment_
In this exercise you will implement a simplified version of the BLAST algorithm.
Please note that the implementation is created for the educational purposes and only mimics the BLAST algorithm.
Also note that in this sheet we will be working with the protein sequences.

**a)** In order to start searching for the hits we need to get the similarity between amino acids.
The first helper function will read the BLOSUM62 matrix from the file `BLOSUM62.txt` (https://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt) and return a dictionary with the similarity between amino acids.
The keys of the dictionary are the amino acids and the values are sub-dictionaries. The sub-dictionaries have the amino acids as keys and the similarity as values.


<details>
  <summary>Example: (Spoiler)</summary>

  ```
   >>> blosum62 = read_blosum62("BLOSUM62.txt")
   >>> print(blosum62)
   {'A': {'A': 4, 'R': -1, 'N': -2, 'D': -2, 'C': 0, 'Q': -1, 'E': -1, 'G': 0, 'H': -2, 'I': -1, 'L': -1, 'K': -1, 'M': -1, 'F': -2, 'P': -1, 'S': 1, 'T': 0, 'W': -3, 'Y': -2, 'V': 0, 'B': -2, 'Z': -1, 'X': 0, '*': -4}, 'R': {'A': -1, 'R': 5, 'N': 0, 'D': -2, 'C': -3, ...
  ```

</details>

if you are puzzled with this step you can simply import the `convert_blosum_txt_to_dict_correct` from the helpers file
```angular2html
from helpers import convert_blosum_txt_to_dict_correct
```

**b)** In order to search for the similar sub-sequences between the database and the query we need to be able to index a sequence.
Here you are asked to implement the `index_sequence_by_kmers` function. The function takes the sequence and the kmer length and returns a dictionary with the kmers as keys and the positions (list of indexes of kmers starts) as values. In biology the numeration starts with 1.

<details>
  <summary>Example: (Spoiler)</summary>

  ```
    >>>database = "DPPEGVVDPP"
    >>>query = "RPPQGLF"

    >>>indexes_db = index_sequence_by_kmers(database, 3)
    >>>indexes_query = index_sequence_by_kmers(query, 3)
    
    >>>print(indexes_db)
    >>>print(indexes_query)
    
    {'DPP': [1, 8], 'PPE': [2], 'PEG': [3], 'EGV': [4], 'GVV': [5], 'VVD': [6], 'VDP': [7]}
    {'RPP': [1], 'PPQ': [2], 'PQG': [3], 'QGL': [4], 'GLF': [5]}
  
  ```

</details>


**c)** Now we need to search for similar kmers in the database and the query. In order to do so lets implement a couple of helpers.

**c0)** List all existing kmers(use both database and query). Note that you can use the index function from the previous part and just the keys.

**c1)** First lets try to list all kmers which are similar to the given one within the given threshold.
In order to do that implement the `find_similar_kmers_for_kmer` which takes four arguments: kmer, kmer similarity threshold, all existing kmers and blosum dictionary.

<details>
  <summary>Example: (Spoiler)</summary>

  ```
    >>>all_similar_kmers = find_similar_kmers_for_kmer("PQG", 13, all_existing_kmers, dict_blosum)
    >>>print(all_similar_kmers)
 
    ['PEG', 'PQG']
  
  ```

</details>


**c2)** Now we can implement the function which will of all similar kmers in the sequence for each kmer in the sequence.
In order to do that implement the `find_similar_kmers_for_sequence` which takes three arguments: sequence, kmer similarity threshold and blosum dictionary.
This function should return a dictionary with the kmer as key and the list of similar kmers as value.


<details>
  <summary>Example: (Spoiler)</summary>

  ```
    >>>all_similar_kmers_for_sequence = find_similar_kmers_for_sequence_correct("DPPEGVVDPP", 13, all_existing_kmers, dict_blosum)
    >>>print(all_similar_kmers_for_sequence)
 
    {'DPP': ['DPP'], 'PPE': ['PPE', 'PPQ'], 'PEG': ['PEG', 'PQG'], 'EGV': ['EGV'], 'GVV': ['GVV'], 'VVD': ['VVD'], 'VDP': ['VDP']}
  
  ```

</details>

**c3)** Finally we are ready to find the hits in both database and query sequences.
In order to do that implement the `create_index_pairs_correct` function. The function has 5 arguments as the input: query sequence, database sequence, kmer_size, kmer similarity threshold and the blosum dictionary.
The function returns a dictionary with kmers as keys and tuples of indexes as values. Each tuple has the list of indexes of the kmer in the query sequence as the first value and a list of indexes of the similar tuples is the database as the second value.

<details>
  <summary>Example: (Spoiler)</summary>

  ```
    >>>dict_both_indexes = create_index_pairs_correct(query, database, 3, 5, dict_blosum)
    >>>print(dict_both_indexes)print(dict_both_indexes)
    
    {'RPP': ([1], [1, 8]), 'PPQ': ([2], [1, 8, 2]), 'PQG': ([3], [3]), 'QGL': ([4], [4]), 'GLF': ([5], [5])}
  ```

</details>
