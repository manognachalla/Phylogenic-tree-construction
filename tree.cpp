#include "tree.hpp"

Tree::Tree(const sequence& sequences) {
    node root;
    root.id = sequences.seq.size();
    root.parent = -1;
    root.level = 0;
    root.isleaf = false;
    root.name = "root";
    root.subtree = "";
    tree.push_back(root);

    for (int i = 0; i < sequences.seq.size(); i++) {
        node leaf;
        leaf.id = i;
        leaf.parent = sequences.seq.size();
        leaf.level = 1;
        leaf.isleaf = true;
        leaf.name = sequences.name[i];
        leaf.subtree = sequences.name[i];
        tree.push_back(leaf);
    }
}

Tree::Tree(const std::vector<dmatrix_row>& D, const std::vector<std::string>& names) {
    node root;
    root.id = D.size();
    root.parent = -1;
    root.level = 0;
    root.isleaf = false;
    root.name = "root";
    root.subtree = "";
    tree.push_back(root);

    for (int i = 0; i < D.size(); i++) {
        node leaf;
        leaf.id = i;
        leaf.parent = D.size();
        leaf.level = 1;
        leaf.isleaf = true;
        leaf.name = names[i];
        leaf.subtree = names[i];
        tree.push_back(leaf);
    }
}

void Tree::joinNodes(int child1, int child2, float child1_distance, float child2_distance) {
    node new_node;
    new_node.id = tree.size();
    new_node.parent = tree[child1].parent;
    new_node.level = tree[child1].level;
    new_node.isleaf = false;
    new_node.child1 = child1;
    new_node.child2 = child2;
    new_node.child1_distance = child1_distance;
    new_node.child2_distance = child2_distance;
    new_node.name = "(" + tree[child1].name + "," + tree[child2].name + ")";

    tree[child1].parent = new_node.id;
    tree[child1].level++;
    tree[child2].parent = new_node.id;
    tree[child2].level++;

    new_node.subtree = "(" + tree[child1].subtree + ":" + std::to_string(child1_distance) + 
                      "," + tree[child2].subtree + ":" + std::to_string(child2_distance) + ")";
    
    tree.push_back(new_node);
    newick = new_node.subtree + ";";
}
