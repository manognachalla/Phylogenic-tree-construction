#include "tree.h"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <iostream>
#include <vector>
#include <queue>
#include <tuple>
#include <limits>



using NodePair = std::tuple<float, int, int>;  // (Q-value, index1, index2)

// Priority queue: Min-heap to store Q-matrix values
std::priority_queue<NodePair, std::vector<NodePair>, std::greater<NodePair>> pq;

void compute_initial_q_matrix(std::vector<dmatrix_row>& D) {
    int n = D.size();
    
    // Populate the priority queue with initial Q-values
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            float q = (n - 2) * D[i].distances[j] - D[i].sum - D[j].sum;
            pq.push(std::make_tuple(q, i, j));
        }
    }
}

void neighbor_joining(std::vector<dmatrix_row>& D, Tree& tree, bool verbose) {
    int n = D.size();
    std::vector<bool> active(n, true);  // Tracks active nodes

    compute_initial_q_matrix(D);

    int iterations = 0;
    while (n > 2) {
        float min_q;
        int min_i, min_j;

        // Extract the best pair (lazy deletion: ignore invalid pairs)
        do {
            std::tie(min_q, min_i, min_j) = pq.top();
            pq.pop();
        } while (!active[min_i] || !active[min_j]);  // Ensure pair is valid

        float d_ij = D[min_i].distances[min_j];
        float d_i = (d_ij + (D[min_i].sum - D[min_j].sum) / (n - 2)) / 2;
        float d_j = d_ij - d_i;

        // Merge nodes in the tree
        tree.joinNodes(D[min_i].id, D[min_j].id, d_i, d_j);

        active[min_j] = false;  // Mark node j as inactive
        D[min_i].id = tree.tree.size() - 1;  // Assign new merged node ID

        // Update distances for the new merged node
        for (int k = 0; k < D.size(); k++) {
            if (!active[k] || k == min_i) continue;

            float d_ik = D[min_i].distances[k];
            float d_jk = D[min_j].distances[k];
            float d_new = (d_ik + d_jk - d_ij) / 2;

            D[min_i].distances[k] = d_new;
            D[k].distances[min_i] = d_new;
        }

        // Update sum of distances for the new node
        D[min_i].sum = 0;
        for (int k = 0; k < D.size(); k++) {
            if (active[k] && k != min_i) {
                D[min_i].sum += D[min_i].distances[k];
            }
        }

        // Push updated Q-values into the priority queue
        for (int k = 0; k < D.size(); k++) {
            if (!active[k] || k == min_i) continue;

            float q = (n - 2) * D[min_i].distances[k] - D[min_i].sum - D[k].sum;
            pq.push(std::make_tuple(q, min_i, k));
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
