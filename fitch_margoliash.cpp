#include "tree.hpp"
#include <algorithm>
#include <limits>
#include <iostream>
#include <map>
#include <cmath>

float calculate_tree_fit(const Tree& tree, const std::vector<dmatrix_row>& D) {
    float fit = 0.0f;
    for (int i = 0; i < D.size(); i++) {
        for (int j = 0; j < i; j++) {  // Only use lower triangle
            float diff = D[i].distances[j] - D[j].distances[i];  // Should be close to 0
            fit += diff * diff;
        }
    }
    return fit;
}

void optimize_branch_lengths(Tree& tree, const std::vector<dmatrix_row>& D) {
    // Simple gradient descent for branch length optimization
    float learning_rate = 0.01f;
    int max_iterations = 100;
    
    for (int iter = 0; iter < max_iterations; iter++) {
        float old_fit = calculate_tree_fit(tree, D);
        
        // Update each branch length
        for (auto& node : tree.tree) {
            if (!node.isleaf) {
                // Try adjusting child1 distance
                float old_dist = node.child1_distance;
                node.child1_distance *= (1.0f + learning_rate);
                float new_fit = calculate_tree_fit(tree, D);
                
                if (new_fit > old_fit) {
                    // Try decreasing instead
                    node.child1_distance = old_dist * (1.0f - learning_rate);
                    new_fit = calculate_tree_fit(tree, D);
                    if (new_fit > old_fit) {
                        // Revert if no improvement
                        node.child1_distance = old_dist;
                    }
                }
                
                // Same for child2 distance
                old_dist = node.child2_distance;
                node.child2_distance *= (1.0f + learning_rate);
                new_fit = calculate_tree_fit(tree, D);
                
                if (new_fit > old_fit) {
                    node.child2_distance = old_dist * (1.0f - learning_rate);
                    new_fit = calculate_tree_fit(tree, D);
                    if (new_fit > old_fit) {
                        node.child2_distance = old_dist;
                    }
                }
            }
        }
    }
}

void fitch_margoliash(std::vector<dmatrix_row>& D, Tree& tree, bool verbose) {
    int n = D.size();
    std::map<int, int> current_to_original;  // Maps current indices to original node IDs
    std::vector<int> active_indices(n);
    
    // Initialize the mapping and active indices
    for (int i = 0; i < n; i++) {
        current_to_original[i] = i;
        active_indices[i] = i;
    }
    
    // Initialize the distance matrix to ensure symmetry
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
            D[j].distances[i] = D[i].distances[j];
        }
    }
    
    while (D.size() > 1) {
        // Find minimum weighted distance pair
        double min_dist = std::numeric_limits<double>::max();
        int min_i = -1, min_j = -1;
        
        for (int i = 0; i < D.size(); i++) {
            for (int j = 0; j < i; j++) {
                double weighted_dist = D[i].distances[j] / std::sqrt(D[i].distances[j]);
                if (weighted_dist < min_dist) {
                    min_dist = weighted_dist;
                    min_i = i;
                    min_j = j;
                }
            }
        }
        
        if (verbose) {
            std::cout << "Merging nodes " << tree.tree[current_to_original[min_i] + 1].name 
                     << " and " << tree.tree[current_to_original[min_j] + 1].name 
                     << " (distance = " << min_dist << ")\n";
            std::cout << "Current matrix size: " << D.size() << std::endl;
        }
        
        // Calculate branch lengths
        float dist_i = min_dist / 2.0f;
        float dist_j = min_dist / 2.0f;
        
        // Join the nodes in the tree
        int orig_i = current_to_original[min_i];
        int orig_j = current_to_original[min_j];
        tree.joinNodes(orig_i, orig_j, dist_i, dist_j);
        int new_node_id = tree.tree.size() - 1;
        
        // Create new distance matrix
        std::vector<dmatrix_row> new_D;
        new_D.reserve(D.size() - 1);
        
        // Copy unmerged rows and update their distances
        for (int i = 0; i < D.size(); i++) {
            if (i != min_i && i != min_j) {
                dmatrix_row new_row;
                new_row.id = D[i].id;
                new_row.distances.reserve(D.size() - 2);  // Reserve space to avoid reallocations
                
                // Copy distances for unmerged nodes
                for (int j = 0; j < D.size(); j++) {
                    if (j != min_i && j != min_j) {
                        new_row.distances.push_back(D[i].distances[j]);
                    }
                }
                
                // Add distance to new merged node using weighted average
                float new_dist = (D[i].distances[min_i] + D[i].distances[min_j]) / 2.0f;
                new_row.distances.push_back(new_dist);
                new_D.push_back(new_row);
            }
        }
        
        // Add the new merged node
        dmatrix_row new_row;
        new_row.id = new_node_id;
        new_row.distances.reserve(new_D.size());  // Reserve space to avoid reallocations
        
        for (const auto& row : new_D) {
            new_row.distances.push_back(row.distances.back());
        }
        new_row.distances.push_back(0.0f);  // Distance to itself
        new_D.push_back(new_row);
        
        // Update active indices
        std::vector<int> new_active_indices;
        new_active_indices.reserve(active_indices.size() - 1);  // Reserve space
        for (int i = 0; i < active_indices.size(); i++) {
            if (i != min_i && i != min_j) {
                new_active_indices.push_back(active_indices[i]);
            }
        }
        new_active_indices.push_back(new_node_id);
        active_indices = new_active_indices;
        
        // Update the mapping
        current_to_original.clear();
        for (int i = 0; i < active_indices.size(); i++) {
            current_to_original[i] = active_indices[i];
        }
        
        // Replace old distance matrix
        D = std::move(new_D);  // Use move semantics to avoid copying
    }
    
    // Final optimization of branch lengths
    if (verbose) {
        std::cout << "Optimizing branch lengths...\n";
    }
    optimize_branch_lengths(tree, D);
}

void fitch_margoliash_tree(std::vector<dmatrix_row>& D, std::string output, bool verbose) {
    std::vector<std::string> names;
    names.reserve(D.size());  // Reserve space
    for (int i = 0; i < D.size(); i++) {
        names.push_back(std::to_string(i));
    }
    Tree tree(D, names);
    fitch_margoliash(D, tree, verbose);
    std::vector<std::string> to_write = {tree.newick};
    write_to_file(output, to_write);
}
