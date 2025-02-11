#include "tree.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <limits>

// Calculate the tree fit using least squares criterion
float calculate_tree_fit(const Tree& tree, const std::vector<dmatrix_row>& D) {
    float fit = 0.0;
    int n = D.size();
    
    // Calculate pairwise distances in the current tree
    std::vector<std::vector<float>> tree_distances(n, std::vector<float>(n, 0.0));
    
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            // Calculate path length between leaves i and j in the tree
            float distance = 0.0;
            int current_i = i;
            int current_j = j;
            std::vector<int> path_i, path_j;
            
            // Trace path from i to root
            while (current_i != -1) {
                path_i.push_back(current_i);
                current_i = tree.tree[current_i].parent;
            }
            
            // Trace path from j to root
            while (current_j != -1) {
                path_j.push_back(current_j);
                current_j = tree.tree[current_j].parent;
            }
            
            // Find lowest common ancestor
            int lca = -1;
            for (int pi : path_i) {
                for (int pj : path_j) {
                    if (pi == pj) {
                        lca = pi;
                        break;
                    }
                }
                if (lca != -1) break;
            }
            
            // Sum up distances
            current_i = i;
            while (current_i != lca) {
                if (tree.tree[current_i].parent != -1) {
                    const auto& parent_node = tree.tree[tree.tree[current_i].parent];
                    if (parent_node.child1 == current_i)
                        distance += parent_node.child1_distance;
                    else
                        distance += parent_node.child2_distance;
                }
                current_i = tree.tree[current_i].parent;
            }
            
            current_j = j;
            while (current_j != lca) {
                if (tree.tree[current_j].parent != -1) {
                    const auto& parent_node = tree.tree[tree.tree[current_j].parent];
                    if (parent_node.child1 == current_j)
                        distance += parent_node.child1_distance;
                    else
                        distance += parent_node.child2_distance;
                }
                current_j = tree.tree[current_j].parent;
            }
            
            tree_distances[i][j] = tree_distances[j][i] = distance;
            
            // Add to fit criterion (weighted least squares)
            float observed = D[i].distances[j];
            float weight = 1.0 / (observed * observed); // Fitch-Margoliash weighting
            fit += weight * std::pow(observed - distance, 2);
        }
    }
    
    return fit;
}

// Optimize branch lengths using least squares
void optimize_branch_lengths(Tree& tree, const std::vector<dmatrix_row>& D) {
    const float epsilon = 1e-6;
    const int max_iterations = 100;
    int iterations = 0;
    float prev_fit = std::numeric_limits<float>::max();
    
    while (iterations < max_iterations) {
        for (size_t node_idx = 0; node_idx < tree.tree.size(); ++node_idx) {
            if (!tree.tree[node_idx].isleaf) {
                // Store original values
                float original_child1 = tree.tree[node_idx].child1_distance;
                float original_child2 = tree.tree[node_idx].child2_distance;
                
                // Simple hill climbing
                float best_fit = calculate_tree_fit(tree, D);
                float best_child1 = original_child1;
                float best_child2 = original_child2;
                
                for (float delta : {-0.1f, 0.1f}) {
                    // Try adjusting child1 distance
                    tree.tree[node_idx].child1_distance = std::max(0.0f, original_child1 + delta);
                    float fit = calculate_tree_fit(tree, D);
                    if (fit < best_fit) {
                        best_fit = fit;
                        best_child1 = tree.tree[node_idx].child1_distance;
                    }
                    
                    // Try adjusting child2 distance
                    tree.tree[node_idx].child1_distance = original_child1;
                    tree.tree[node_idx].child2_distance = std::max(0.0f, original_child2 + delta);
                    fit = calculate_tree_fit(tree, D);
                    if (fit < best_fit) {
                        best_fit = fit;
                        best_child2 = tree.tree[node_idx].child2_distance;
                    }
                }
                
                tree.tree[node_idx].child1_distance = best_child1;
                tree.tree[node_idx].child2_distance = best_child2;
            }
        }
        
        float current_fit = calculate_tree_fit(tree, D);
        if (std::abs(current_fit - prev_fit) < epsilon) {
            break;
        }
        prev_fit = current_fit;
        iterations++;
    }
}

void fitch_margoliash(std::vector<dmatrix_row>& D, Tree& tree, bool verbose) {
    // Start with a star tree
    int n = D.size();
    
    // Initialize with two nodes joined
    float min_distance = std::numeric_limits<float>::max();
    int min_i = 0, min_j = 1;
    
    // Find closest pair
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (D[i].distances[j] < min_distance) {
                min_distance = D[i].distances[j];
                min_i = i;
                min_j = j;
            }
        }
    }
    
    // Join first two nodes
    tree.joinNodes(min_i, min_j, min_distance/2, min_distance/2);
    
    // Add remaining nodes one by one
    std::vector<bool> added(n, false);
    added[min_i] = added[min_j] = true;
    
    for (int iter = 2; iter < n; iter++) {
        if (verbose && iter % 10 == 0) {
            std::cout << "Adding node " << iter << " of " << n << std::endl;
        }
        
        // Find next node to add
        float min_total_error = std::numeric_limits<float>::max();
        int best_node = -1;
        float best_distance = 0;
        int best_attachment = -1;
        
        for (int i = 0; i < n; i++) {
            if (added[i]) continue;
            
            // Try attaching to each existing branch
            for (size_t node_idx = 0; node_idx < tree.tree.size(); ++node_idx) {
                if (tree.tree[node_idx].isleaf) continue;
                
                // Try different attachment points
                float original_child1 = tree.tree[node_idx].child1_distance;
                float original_child2 = tree.tree[node_idx].child2_distance;
                
                // Temporarily attach the new node
                tree.joinNodes(i, node_idx, D[i].distances[node_idx]/2, D[i].distances[node_idx]/2);
                float error = calculate_tree_fit(tree, D);
                
                if (error < min_total_error) {
                    min_total_error = error;
                    best_node = i;
                    best_distance = D[i].distances[node_idx]/2;
                    best_attachment = node_idx;
                }
                
                // Restore original tree
                tree.tree.pop_back();
                tree.tree[node_idx].child1_distance = original_child1;
                tree.tree[node_idx].child2_distance = original_child2;
            }
        }
        
        // Add the best node
        tree.joinNodes(best_node, best_attachment, best_distance, best_distance);
        added[best_node] = true;
        
        // Optimize branch lengths
        optimize_branch_lengths(tree, D);
    }
}

void fitch_margoliash_tree(std::vector<dmatrix_row>& D, std::string output, bool verbose) {
    Tree tree(D);
    fitch_margoliash(D, tree, verbose);
    std::vector<std::string> to_write = {tree.newick};
    write_to_file(output, to_write);
}
