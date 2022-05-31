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
    >>>all_similar_kmers_for_sequence = find_similar_kmers_for_sequence("DPPEGVVDPP", 13, all_existing_kmers, dict_blosum)
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
    >>>dict_both_indexes = create_index_pairs(query, database, 3, 5, dict_blosum)
    >>>print(dict_both_indexes)
    
    {'RPP': ([1], [1, 8]), 'PPQ': ([2], [1, 8, 2]), 'PQG': ([3], [3]), 'QGL': ([4], [4]), 'GLF': ([5], [5])}
  ```

</details>

**d)** In this part we will be extending the kmer hits. In other words we will be building extended hits which represent longer interval of similar sequences in both query and database sequences.
We agreed to have a single hit a tuple of two indexes (query index and database index). The extended hit will be nothing but a list of such tuples.


**d1)** In order to start building the extended hits we first need a way to merge two single hits.
In order to do that implement the `merge_single_hit_with_single_hit` function which takes two arguments: hit1 and hit2 and the max_distance threshold.
The function returns two values. The first value is the flag which is True if the hits were merged and False if they were not. The second value is the extended hit as a list if the hits can be merged and None otherwise.
The sum of differences between the indexes of the query sequence and the database sequence should be less than the max distance threshold.
In order to compute the distance use the following formula (we want the distance of adjacent hits to be zero):

```
  >>>distance = abs(hit1[0] - hit2[0]) - 1  + abs(hit1[1] - hit2[1]) - 1
```

<details>
  <summary>Example: (Spoiler)</summary>

  ```
    >>>extension, extedned_hit = merge_single_hit_with_single_hit_correct((1, 1), (2, 2), 2)
    >>>print(extension, extedned_hit)
    
    True [(1, 1), (2, 2)]
    
    >>>extension, extedned_hit = merge_single_hit_with_single_hit_correct((1, 1), (2, 8), 2)
    >>>print(extension, extedned_hit)
    
    False None
    
    
  ```

</details>

**d2)** Now when we know how to merge two single hits we can try to merge an extended hit with a single hit.
In order to do that implement the `merge_extended_hit_with_single_hit` function which takes three arguments: extended hit, single hit and max distance threshold.
The function returns two values. The first value is the flag which is True if the hits were merged and False if they were not. The second value is the extended hit as a list if the hits can be merged and None otherwise.


**d3)** Now when we know how to merge two single hits and an extended hit with a single hit we can try to merge an extended hit with an extended hit.
In order to do that implement the `merge_two_extended_hits` function which takes three arguments: extended hit1, extended hit2 and max distance threshold.
The function returns two values. The first value is the flag which is True if the hits were merged and False if they were not. The second value is the extended hit as a list if the hits can be merged and None otherwise.

**d4)** Now when we know how to merge hits of all types implement a function which tries to merge two hits(they can be both single or extended).
In order to do that implement the `merge_two_hits` function which takes three arguments: hit1, hit2 and max distance threshold.
The function returns two values. The first value is the flag which is True if the hits were merged and False if they were not. The second value is the extended hit as a list if the hits can be merged and None otherwise.


**e)** Now when all the helper functions are implemented we can finally assemble all the extended hits.
In order to do that implement the `create_extended_hits_correct` function which takes six arguments: query sequence, database sequence, kmer size, kmer similarity threshold, blosum dictionary, max distance threshold.
Use the helper functions from the previous steps.




We will be taking care of the statistical significance of the results in the next sheet.


