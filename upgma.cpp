#include "tree.hpp"
#include <algorithm>
#include <limits>
#include <iostream>

void upgma(std::vector<dmatrix_row>& D, Tree& tree, bool verbose) {
    int n = D.size();
    std::vector<float> heights(n, 0.0f);  // Height of each node from leaves
    std::vector<int> cluster_size(n, 1);  // Size of each cluster

    for (int step = 0; step < n - 1; step++) {
        // Find minimum distance
        float min_dist = std::numeric_limits<float>::max();
        int min_i = -1, min_j = -1;
        
        for (int i = 0; i < D.size(); i++) {
            for (int j = 0; j < i; j++) {
                if (D[i].distances[j] < min_dist) {
                    min_dist = D[i].distances[j];
                    min_i = i;
                    min_j = j;
                }
            }
        }

        if (verbose) {
            std::cout << "Step " << step + 1 << ": Joining clusters " << min_i << " and " << min_j 
                      << " (distance = " << min_dist << ")\n";
        }

        // Calculate height for the new node
        float new_height = min_dist / 2.0f;

        // Join the nodes in the tree
        float dist_i = new_height - heights[min_i];
        float dist_j = new_height - heights[min_j];
        tree.joinNodes(min_i, min_j, dist_i, dist_j);

        // Update distance matrix
        std::vector<float> new_distances;
        for (int k = 0; k < D.size(); k++) {
            if (k != min_i && k != min_j) {
                // Calculate new distance using weighted average
                float new_dist = (D[min_i].distances[k] * cluster_size[min_i] + 
                                D[min_j].distances[k] * cluster_size[min_j]) / 
                               (cluster_size[min_i] + cluster_size[min_j]);
                new_distances.push_back(new_dist);
            }
        }

        // Update cluster sizes
        cluster_size.push_back(cluster_size[min_i] + cluster_size[min_j]);

        // Remove larger index first to avoid shifting issues
        if (min_i > min_j) {
            D.erase(D.begin() + min_i);
            D.erase(D.begin() + min_j);
        } else {
            D.erase(D.begin() + min_j);
            D.erase(D.begin() + min_i);
        }

        // Add new row
        dmatrix_row new_row;
        new_row.distances = new_distances;
        new_row.id = D.size();
        D.push_back(new_row);

        // Update heights
        heights.push_back(new_height);
    }
}

void upgma_tree(std::vector<dmatrix_row>& D, std::string output, bool verbose) {
    std::vector<std::string> names;
    for (int i = 0; i < D.size(); i++) {
        names.push_back(std::to_string(i));
    }
    Tree tree(D, names);
    upgma(D, tree, verbose);
    std::vector<std::string> to_write = {tree.newick};
    write_to_file(output, to_write);
}
