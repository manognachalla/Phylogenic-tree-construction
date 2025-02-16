#include "tree.h"
#include <algorithm>
#include <limits>
#include <iostream>
#include <map>

void minimum_evolution(std::vector<dmatrix_row>& D, Tree& tree, bool verbose) {
    int n = D.size();
    std::map<int, std::string> labels;
    
    // Initialize labels with node names from the tree
    for (int i = 0; i < n; i++) {
        labels[i] = tree.tree[i].name;
    }
    
    while (D.size() > 1) {
        // Find minimum distance pair
        double min_dist = std::numeric_limits<double>::max();
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
            std::cout << "Merging nodes " << labels[D[min_i].id] << " and " 
                      << labels[D[min_j].id] << " (distance = " << min_dist << ")\n";
        }

        // Calculate branch lengths
        float dist_i = min_dist / 2.0f;
        float dist_j = min_dist / 2.0f;

        // Join the nodes in the tree
        tree.joinNodes(D[min_i].id, D[min_j].id, dist_i, dist_j);
        
        // Update distance matrix
        std::vector<float> new_distances;
        for (int k = 0; k < D.size(); k++) {
            if (k != min_i && k != min_j) {
                // Use simple average for new distances
                float new_dist = (D[min_i].distances[k] + D[min_j].distances[k]) / 2.0f;
                new_distances.push_back(new_dist);
            }
        }
        
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
    }
}

void minimum_evolution_tree(std::vector<dmatrix_row>& D, std::string output, bool verbose) {
    Tree tree(D);
    minimum_evolution(D, tree, verbose);
    std::vector<std::string> to_write = {tree.newick};
    write_to_file(output, to_write);
}
