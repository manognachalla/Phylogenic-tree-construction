#include "tree.h"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <iostream>

std::vector<std::vector<float>> count_kmer_frequencies(sequence& sequences, int& kmer_length) {
    std::cout << "Reading sequences, counting K-mers of length: " << kmer_length << "..." << std::endl;
    std::map<std::string, int> unique_kmers;
    std::vector<std::vector<float>> kmer_frequencies(sequences.seq.size());
    int index, flag;
    int unique_counter = 0;

    for (int i = 0; i < sequences.seq.size(); i++) {
        for (int j = 0; j <= sequences.seq[i].length() - kmer_length; j++) {
            flag = 0;
            std::string kmer = sequences.seq[i].substr(j, kmer_length);
            for (int s = 0; s < kmer_length; s++) {
                if (kmer[s] != 'A' && kmer[s] != 'C' && kmer[s] != 'G' && kmer[s] != 'T') {
                    flag = 1;
                    break;
                }
            }

            if (!unique_kmers[kmer] && flag != 1) {
                unique_kmers[kmer] = unique_counter;
                for (int k = 0; k < kmer_frequencies.size(); k++) {
                    kmer_frequencies[k].push_back(0);
                }
                kmer_frequencies[i][unique_counter] += 1;
                unique_counter += 1;
            }
            else if (flag != 1) {
                index = unique_kmers[kmer];
                kmer_frequencies[i][index] += 1;
            }
        }
    }
    return kmer_frequencies;
}

std::vector<dmatrix_row> distance_matrix(std::vector<std::vector<float>>& frequencies, sequence& sequences, int kmer_length, std::string method) {
    std::vector<dmatrix_row> D(frequencies.size());
    float distance;

    for (int i = 0; i < frequencies.size(); i++) {
        D[i].distances.resize(frequencies.size());
        D[i].id = i;
        D[i].sum = 0;

        for (int j = 0; j < frequencies.size(); j++) {
            if (i == j) {
                D[i].distances[j] = 0;
                continue;
            }

            if (method == "cosine") {
                float dot = 0, norm1 = 0, norm2 = 0;
                for (int k = 0; k < frequencies[i].size(); k++) {
                    dot += frequencies[i][k] * frequencies[j][k];
                    norm1 += frequencies[i][k] * frequencies[i][k];
                    norm2 += frequencies[j][k] * frequencies[j][k];
                }
                distance = 1 - (dot / (std::sqrt(norm1) * std::sqrt(norm2)));
            }
            else {
                float total1 = 0, total2 = 0;
                for (float f : frequencies[i]) total1 += f;
                for (float f : frequencies[j]) total2 += f;

                std::vector<float> freq1 = frequencies[i];
                std::vector<float> freq2 = frequencies[j];
                
                for (int k = 0; k < freq1.size(); k++) {
                    freq1[k] /= total1;
                    freq2[k] /= total2;
                }

                if (method == "mahalanobis") {
                    distance = 0;
                    for (int k = 0; k < freq1.size(); k++) {
                        if (freq1[k] + freq2[k] > 0) {
                            distance += std::pow(freq1[k] - freq2[k], 2) / (freq1[k] + freq2[k]);
                        }
                    }
                    distance = std::sqrt(distance);
                }
                else {
                    distance = 0;
                    for (int k = 0; k < freq1.size(); k++) {
                        distance += std::abs(freq1[k] - freq2[k]);
                    }
                }
            }
            D[i].distances[j] = distance;
            D[i].sum += distance;
        }
    }
    return D;
}

void neighbor_joining(std::vector<dmatrix_row>& D, Tree& tree, bool verbose) {
    int n = D.size();
    std::vector<bool> active(n, true);
    int iterations = 0;

    while (n > 2) {
        float min_q = std::numeric_limits<float>::max();
        int min_i = -1, min_j = -1;

        for (int i = 0; i < D.size(); i++) {
            if (!active[i]) continue;
            for (int j = i + 1; j < D.size(); j++) {
                if (!active[j]) continue;

                float q = (n - 2) * D[i].distances[j] - D[i].sum - D[j].sum;
                if (q < min_q) {
                    min_q = q;
                    min_i = i;
                    min_j = j;
                }
            }
        }

        float d_ij = D[min_i].distances[min_j];
        float d_i = (d_ij + (D[min_i].sum - D[min_j].sum) / (n - 2)) / 2;
        float d_j = d_ij - d_i;

        tree.joinNodes(D[min_i].id, D[min_j].id, d_i, d_j);

        active[min_j] = false;
        D[min_i].id = tree.tree.size() - 1;

        for (int k = 0; k < D.size(); k++) {
            if (!active[k] || k == min_i) continue;
            float d_ik = D[min_i].distances[k];
            float d_jk = D[min_j].distances[k];
            float d_new = (d_ik + d_jk - d_ij) / 2;
            
            D[min_i].distances[k] = d_new;
            D[k].distances[min_i] = d_new;
        }

        D[min_i].sum = 0;
        for (int k = 0; k < D.size(); k++) {
            if (active[k] && k != min_i) {
                D[min_i].sum += D[min_i].distances[k];
            }
        }

        n--;
        iterations++;
        
        if (verbose && iterations % 100 == 0) {
            std::cout << "Iteration: " << iterations << std::endl;
        }
    }
}

void neighbor_joining_tree(std::vector<dmatrix_row>& D, std::string output, bool verbose) {
    Tree tree(D);
    neighbor_joining(D, tree, verbose);
    std::vector<std::string> to_write = {tree.newick};
    write_to_file(output, to_write);
}
