#include "tree.hpp"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <iostream>
#include <vector>
#include <queue>
#include <tuple>
#include <limits>


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
                if (norm1 == 0 || norm2 == 0) {
                    distance = 1.0;  // Maximum distance for sequences with no k-mers
                } else {
                    distance = 1 - (dot / (std::sqrt(norm1) * std::sqrt(norm2)));
                }
            }
            else {
                float total1 = 0, total2 = 0;
                for (float f : frequencies[i]) total1 += f;
                for (float f : frequencies[j]) total2 += f;

                if (total1 == 0 || total2 == 0) {
                    distance = 1.0;  // Maximum distance for sequences with no k-mers
                    D[i].distances[j] = distance;
                    D[i].sum += distance;
                    continue;
                }

                std::vector<float> freq1 = frequencies[i];
                std::vector<float> freq2 = frequencies[j];
                
                for (int k = 0; k < freq1.size(); k++) {
                    freq1[k] = total1 > 0 ? freq1[k] / total1 : 0;
                    freq2[k] = total2 > 0 ? freq2[k] / total2 : 0;
                }

                if (method == "mahalanobis") {
                    distance = 0;
                    for (int k = 0; k < freq1.size(); k++) {
                        float sum = freq1[k] + freq2[k];
                        if (sum > 0) {
                            distance += std::pow(freq1[k] - freq2[k], 2) / sum;
                        }
                    }
                    distance = std::sqrt(distance);
                }
                else {  // fractional
                    distance = 0;
                    for (int k = 0; k < freq1.size(); k++) {
                        distance += std::abs(freq1[k] - freq2[k]);
                    }
                    distance /= 2.0;  // Normalize to [0,1] range
                }
            }
            D[i].distances[j] = distance;
            D[i].sum += distance;
        }
    }
    return D;
}