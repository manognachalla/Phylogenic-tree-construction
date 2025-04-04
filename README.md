# Phylogenetic Tree Construction Tool

A C++ program for constructing phylogenetic trees from sequence data using Neighbor-Joining, Fitch-Margoliash, UPGMA, or Minimum Evolution algorithms. The program supports FASTA  format sequence files and can also generate random trees for testing purposes.

## Features

- Four tree construction algorithms:
  - Neighbor-Joining (NJ): Fast and widely used method
  - Fitch-Margoliash (FM): More accurate branch lengths through weighted least squares optimization
  - UPGMA: Produces ultrametric trees with molecular clock assumption
  - Minimum Evolution (ME): Optimizes tree topology based on minimum total branch length
- Multiple distance calculation methods:
  - Fractional k-mer count (default)
  - Mahalanobis distance
  - Cosine distance
- Input formats:
  - FASTA format
  - Random distance matrix generation

## Compilation

To compile the program, use g++ with C++11 support:

```bash
g++ -o phylo_tree main.cpp tree.cpp neighbor_joining.cpp fitch_margoliash.cpp upgma.cpp minimum_evolution.cpp tree_io.cpp operations.cpp eval.cpp -std=c++11
```

## Usage

### Basic Usage

```bash
./phylo_tree <input_file> [options]
```

### Command Line Options

- Algorithm Selection:
  - `-nj` : Use Neighbor-Joining algorithm (default)
  - `-fm` : Use Fitch-Margoliash algorithm
  - `-upgma` : Use UPGMA algorithm
  - `-me` : Use Minimum Evolution algorithm

- Distance Calculation Methods:
  - `-m` : Use Mahalanobis distance
  - `-c` : Use Cosine distance
  - (default: fractional k-mer count)

- K-mer Length:
  - `-k <INT>` : Set k-mer length (default: 8)

- Verbose Output:
  - `-v` : Enable verbose output

### Examples

1. Basic usage with FASTA file (using default Neighbor-Joining):
```bash
./phylo_tree sequences.fasta
```

2. Using Fitch-Margoliash algorithm:
```bash
./phylo_tree sequences.fasta -fm
```

3. Using UPGMA algorithm:
```bash
./phylo_tree sequences.fasta -upgma
```

4. Using Minimum Evolution algorithm:
```bash
./phylo_tree sequences.fasta -me
```

5. Using Mahalanobis distance with k-mer length 6:
```bash
./phylo_tree sequences.fasta -m -k 6
```

6. Generate a random tree with 10 leaves using UPGMA:
```bash
./phylo_tree -random 10 -upgma
```

### Output

The program generates a Newick format tree file named `output.txt` in the current directory. For PAML files with multiple replicates, the file will contain one tree per line.

## Algorithm Details

### Neighbor-Joining (NJ)
- Fast and widely used method
- Constructs trees based on the principle of minimum evolution
- Good for large datasets
- Computationally efficient: O(n³) time complexity

### Fitch-Margoliash (FM)
- Uses weighted least squares optimization
- More accurate branch lengths
- Better for non-ultrametric distances
- Computationally more intensive than NJ
- Features:
  - Branch length optimization after each node addition
  - Iterative tree improvement

### UPGMA (Unweighted Pair Group Method with Arithmetic Mean)
- Produces ultrametric trees (molecular clock assumption)
- Simple and intuitive hierarchical clustering method
- Good for datasets where evolutionary rates are roughly constant
- Features:
  - Guaranteed ultrametric property
  - Height-based tree construction
  - Computationally efficient: O(n²) time complexity
- Best suited for:
  - Closely related sequences
  - Data that follows a molecular clock
  - Situations where ultrametric trees are desired

### Minimum Evolution (ME)
- Optimizes tree topology to minimize total branch length
- Based on the principle that shorter trees are more likely to be correct
- Features:
  - Progressive clustering approach
  - Efficient distance matrix updates
  - Balanced consideration of all possible tree topologies
- Advantages:
  - Good for finding parsimonious trees
  - Effective for datasets with varying evolutionary rates
  - Less sensitive to long branch attraction
- Computationally efficient: O(n³) time complexity

## Distance Calculation Methods

1. **Fractional k-mer Count (Default)**
   - Counts k-mer frequencies and normalizes them
   - Good balance between speed and accuracy

2. **Mahalanobis Distance (-m)**
   - Takes into account the correlation between k-mers
   - Better for sequences with varying composition

3. **Cosine Distance (-c)**
   - Treats k-mer profiles as vectors
   - Good for comparing sequence composition patterns

## Input File Format

### FASTA Format
```
>Sequence1
ATGCTAGCTAGCT
>Sequence2
ATGCTAGCTAGCT
```

## Notes

- For large sequences, increasing k-mer length might improve accuracy but will increase memory usage
- The Fitch-Margoliash algorithm is slower but generally produces more accurate branch lengths
- Use verbose mode (-v) to monitor progress for large datasets
- Output trees are in Newick format and can be visualized using tools like FigTree or iTOL

## Contributors
