<img src="./figures/banner.png" alt="UniFreiburg Banner"/>

Lehrstuhl für Bioinformatik - Institut für Informatik - *http://www.bioinf.uni-freiburg.de*

---
## Bioinformatics 2
###### SS 2022
##### Exercise sheet 6: BLAST
---

### _Exercise 3 - Programming Assignment_
In this exercise you will implement a simplified version of the BLAST algorithm.
In these sheet we will be working with the protein sequences.

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
