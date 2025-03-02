#include "tree.h"
#include <iostream>
#include <random>
#include <ctime>

void help() {
    std::cout << "\nArgument help:\n"
              << "1st argument:\n"
              << "            filename of the sequences ['.fasta'] format.\n"
              << "            or\n"
              << "            [-random INT] : generate a random distance matrix of size INT x INT and create a Newick format tree with INT leaf nodes\n\n"
              << "Additional arguments: \n"
              << "Algorithm selection:\n"
              << "            [-nj] : Neighbor-Joining algorithm (default)\n"
              << "            [-fm] : Fitch-Margoliash algorithm\n"
              << "            [-upgma] : UPGMA algorithm\n"
              << "            [-me] : Minimum Evolution algorithm\n\n"
              << "Methods for calculating the distance matrix based on kmer profiles of sequences:  \n\n"
              << "            [-m] : mahalanobis; \n"
              << "            [-c] : cosine. \n"
              << "            (default: fractional k-mer count)\n\n"
              << "kmer-length (default 8): \n"
              << "            [-k INT]:\n\n"
              << "Number of replicates to parse in .paml files of synthetic sequences (default 1): \n"
              << "            [-replicates INT]\n"
              << "            Outputs INT Newick trees each based on a different set of replicate sequences.\n\n"
              << "Verbose:    [-v]\n";
}

std::vector<dmatrix_row> random_distance_matrix(int size) {
    std::vector<dmatrix_row> D(size);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);

    for (int i = 0; i < size; i++) {
        D[i].distances.resize(size);
        D[i].id = i;
        D[i].sum = 0;
        
        for (int j = 0; j < size; j++) {
            if (i == j) {
                D[i].distances[j] = 0;
            }
            else if (j > i) {
                D[i].distances[j] = dis(gen);
                D[i].sum += D[i].distances[j];
            }
            else {
                D[i].distances[j] = D[j].distances[i];
                D[i].sum += D[i].distances[j];
            }
        }
    }
    return D;
}

void random_newick_tree(int size, std::string algorithm, std::string output, bool verbose) {
    std::vector<dmatrix_row> D = random_distance_matrix(size);
    if (algorithm == "fm") {
        fitch_margoliash_tree(D, output, verbose);
    } else if (algorithm == "upgma") {
        upgma_tree(D, output, verbose);
    } else if (algorithm == "me") {
        minimum_evolution_tree(D, output, verbose);
    } else {
        neighbor_joining_tree(D, output, verbose);
    }
}

void fasta_to_newick(std::string& filename, int kmer_length, std::string method, std::string algorithm, std::string output, bool verbose) {
    sequence sequences = read_fasta(filename);
    std::vector<std::vector<float>> frequencies = count_kmer_frequencies(sequences, kmer_length);
    std::vector<dmatrix_row> D = distance_matrix(frequencies, sequences, kmer_length, method);
    Tree tree(sequences);
    
    if (algorithm == "fm") {
        fitch_margoliash(D, tree, verbose);
    } else if (algorithm == "upgma") {
        upgma(D, tree, verbose);
    } else if (algorithm == "me") {
        minimum_evolution(D, tree, verbose);
    } else {
        neighbor_joining(D, tree, verbose);
    }
    
    std::vector<std::string> to_write = {tree.newick};
    write_to_file(output, to_write);
}


int main(int argc, char** argv) {
    if (argc < 2) {
        help();
        return 1;
    }

    std::string input = argv[1];
    std::string method = "fractional";
    std::string algorithm = "nj";  // default to neighbor-joining
    std::string output = "output.txt";
    int kmer_length = 8;
    int n_replicates = 1;
    bool verbose = false;

    for (int i = 2; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "-m") method = "mahalanobis";
        else if (arg == "-c") method = "cosine";
        else if (arg == "-nj") algorithm = "nj";
        else if (arg == "-fm") algorithm = "fm";
        else if (arg == "-upgma") algorithm = "upgma";
        else if (arg == "-me") algorithm = "me";
        else if (arg == "-k" && i + 1 < argc) kmer_length = std::stoi(argv[++i]);
        else if (arg == "-replicates" && i + 1 < argc) n_replicates = std::stoi(argv[++i]);
        else if (arg == "-v") verbose = true;
    }

    if (input == "-random" && argc > 2) {
        int size = std::stoi(argv[2]);
        random_newick_tree(size, algorithm, output, verbose);
    }
    else {
        fasta_to_newick(input, kmer_length, method, algorithm, output, verbose);
    }

    return 0;
}
