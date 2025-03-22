#include <iostream>
#include <vector>
#include <string>
#include <limits>
#include <cmath>
#include <iomanip>
#include <map>
#include <sstream>
#include <algorithm>
#include <numeric>

using namespace std;

struct TreeNode {
    string name;
    TreeNode* left = nullptr;
    TreeNode* right = nullptr;
    double branch_length_left = 0.0;
    double branch_length_right = 0.0;

    bool is_leaf() const {
        return left == nullptr && right == nullptr;
    }
};

// Recursively build pairwise distances in tree
void traverse_distances(TreeNode* node, vector<pair<string, double>> path,
                        map<string, int>& label_index, vector<vector<double>>& dist_matrix, double acc_dist) {
    if (!node) return;

    if (node->is_leaf()) {
        for (auto& [other, d] : path) {
            int i = label_index[node->name];
            int j = label_index[other];
            dist_matrix[i][j] = dist_matrix[j][i] = acc_dist + d;
        }
        path.emplace_back(node->name, 0.0);
    } else {
        auto path_left = path;
        path_left.emplace_back(node->left->name, node->branch_length_left);
        traverse_distances(node->left, path_left, label_index, dist_matrix, acc_dist + node->branch_length_left);

        auto path_right = path;
        path_right.emplace_back(node->right->name, node->branch_length_right);
        traverse_distances(node->right, path_right, label_index, dist_matrix, acc_dist + node->branch_length_right);
    }
}

// Compute least squares error
double compute_least_squares_error(const vector<vector<double>>& D, TreeNode* tree, const vector<string>& labels) {
    int n = labels.size();
    vector<vector<double>> dist_matrix(n, vector<double>(n, 0.0));
    map<string, int> label_index;
    for (int i = 0; i < n; ++i)
        label_index[labels[i]] = i;

    traverse_distances(tree, {}, label_index, dist_matrix, 0.0);

    double error = 0.0;
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            error += pow(D[i][j] - dist_matrix[i][j], 2);

    return error;
}

// Neighbor joining tree construction
TreeNode* neighbor_joining(vector<vector<double>> D, vector<string> labels) {
    int n = D.size();
    vector<TreeNode*> nodes(n);
    for (int i = 0; i < n; ++i)
        nodes[i] = new TreeNode{labels[i]};

    vector<int> active(n);
    iota(active.begin(), active.end(), 0);

    while (active.size() > 2) {
        int m = active.size();
        vector<double> total_d(n, 0.0);
        for (int i : active)
            for (int j : active)
                total_d[i] += D[i][j];

        vector<vector<double>> Q(n, vector<double>(n, numeric_limits<double>::infinity()));
        for (int i : active) {
            for (int j : active) {
                if (i < j)
                    Q[i][j] = (m - 2) * D[i][j] - total_d[i] - total_d[j];
            }
        }

        int min_i = -1, min_j = -1;
        double min_val = numeric_limits<double>::infinity();
        for (int i : active) {
            for (int j : active) {
                if (i < j && Q[i][j] < min_val) {
                    min_val = Q[i][j];
                    min_i = i;
                    min_j = j;
                }
            }
        }

        TreeNode* new_node = new TreeNode;
        double delta = (total_d[min_i] - total_d[min_j]) / (m - 2);
        double limb_i = 0.5 * D[min_i][min_j] + 0.5 * delta;
        double limb_j = 0.5 * D[min_i][min_j] - 0.5 * delta;

        new_node->left = nodes[min_i];
        new_node->right = nodes[min_j];
        new_node->branch_length_left = limb_i;
        new_node->branch_length_right = limb_j;

        vector<double> new_row(n, 0.0);
        for (int k : active) {
            if (k != min_i && k != min_j)
                new_row[k] = 0.5 * (D[min_i][k] + D[min_j][k] - D[min_i][min_j]);
        }

        for (int k : active) {
            if (k != min_i && k != min_j) {
                D[k][min_i] = D[min_i][k] = new_row[k];
                D[k][min_j] = D[min_j][k] = new_row[k];
            }
        }

        nodes[min_i] = new_node;
        active.erase(remove(active.begin(), active.end(), min_j), active.end());
    }

    TreeNode* root = new TreeNode;
    root->left = nodes[active[0]];
    root->right = nodes[active[1]];
    root->branch_length_left = D[active[0]][active[1]] / 2;
    root->branch_length_right = D[active[0]][active[1]] / 2;

    return root;
}

// Gradient-free branch length optimization
void optimize_branch_lengths(TreeNode* tree, const vector<vector<double>>& D, const vector<string>& labels, double lr = 0.01, int iterations = 100) {
    for (int iter = 0; iter < iterations; ++iter) {
        function<void(TreeNode*)> update_branch = [&](TreeNode* node) {
            if (!node || node->is_leaf()) return;

            double orig_left = node->branch_length_left;
            double orig_right = node->branch_length_right;
            double base_error = compute_least_squares_error(D, tree, labels);

            node->branch_length_left += lr;
            if (compute_least_squares_error(D, tree, labels) > base_error)
                node->branch_length_left = orig_left;

            node->branch_length_right += lr;
            if (compute_least_squares_error(D, tree, labels) > base_error)
                node->branch_length_right = orig_right;

            update_branch(node->left);
            update_branch(node->right);
        };
        update_branch(tree);
    }
}

TreeNode* run_fitch_margoliash(vector<vector<double>> D, vector<string> labels, double& final_error) {
    TreeNode* tree = neighbor_joining(D, labels);
    optimize_branch_lengths(tree, D, labels);
    final_error = compute_least_squares_error(D, tree, labels);
    return tree;
}

// Convert tree to Newick string
string to_newick(TreeNode* node) {
    if (node->is_leaf()) return node->name;

    string left = to_newick(node->left);
    string right = to_newick(node->right);
    ostringstream ss;
    ss << "(" << left << ":" << fixed << setprecision(6) << node->branch_length_left
       << "," << right << ":" << node->branch_length_right << ")";
    return ss.str();
}

int main() {
    vector<string> labels = {"A", "B", "C", "D"};
    vector<vector<double>> D = {
        {0.0, 5.0, 9.0, 9.0},
        {5.0, 0.0, 10.0, 10.0},
        {9.0, 10.0, 0.0, 8.0},
        {9.0, 10.0, 8.0, 0.0}
    };

    double error;
    TreeNode* tree = run_fitch_margoliash(D, labels, error);

    cout << "Fitch-Margoliash least squares error: " << error << endl;
    string newick = to_newick(tree) + ";";
    cout << "Newick format:\n" << newick << endl;

    return 0;
}