#include <iostream>
#include <vector>
#include <cmath>
#include <set>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <cstdlib>
#include <ctime>

using namespace std;

class Node {
public:
    Node* left;
    Node* right;
    int value; // Only used for leaf taxa
    bool isLeaf;
    Node* nearestNeighbor;
    double distToNN;

    Node(int val) : left(nullptr), right(nullptr), value(val), isLeaf(true), nearestNeighbor(nullptr), distToNN(0.0) {}
    Node(Node* l, Node* r) : left(l), right(r), isLeaf(false), nearestNeighbor(nullptr), distToNN(0.0) {}

    vector<Node*> getLeaves() const {
        vector<Node*> leaves;
        if (isLeaf) {
            leaves.push_back(const_cast<Node*>(this));
        } else {
            if (left) {
                vector<Node*> l_leaves = left->getLeaves();
                leaves.insert(leaves.end(), l_leaves.begin(), l_leaves.end());
            }
            if (right) {
                vector<Node*> r_leaves = right->getLeaves();
                leaves.insert(leaves.end(), r_leaves.begin(), r_leaves.end());
            }
        }
        return leaves;
    }

    int countLeaves() const {
        return getLeaves().size();
    }

    double getDistance(Node* other) {
        vector<Node*> a = this->getLeaves();
        vector<Node*> b = other->getLeaves();

        double total = 0.0;
        for (Node* x : a) {
            for (Node* y : b) {
                total += abs(x->value - y->value);
            }
        }

        return total / (a.size() * b.size());
    }

    string toNewick() const {
        if (isLeaf) {
            return to_string(value);
        } else {
            return "(" + left->toNewick() + "," + right->toNewick() + ")";
        }
    }
};

class UPGMA {
public:
    Node* tree;

    UPGMA(const vector<int>& taxa) {
        set<Node*> clusters;
        for (int val : taxa) {
            clusters.insert(new Node(val));
        }

        buildTree(clusters);
    }

    ~UPGMA() {
        deleteTree(tree);
    }

    void buildTree(set<Node*>& clusters) {
        for (Node* node : clusters) {
            updateDistances(node, clusters);
        }

        while (clusters.size() > 1) {
            Node* c1 = nullptr;
            Node* c2 = nullptr;
            double minDist = 1e9;

            for (Node* node : clusters) {
                if (node->distToNN < minDist) {
                    minDist = node->distToNN;
                    c1 = node;
                    c2 = node->nearestNeighbor;
                }
            }

            clusters.erase(c1);
            clusters.erase(c2);

            Node* merged = new Node(c1, c2);
            clusters.insert(merged);

            updateDistances(merged, clusters);

            for (Node* node : clusters) {
                if (node->nearestNeighbor == c1 || node->nearestNeighbor == c2) {
                    updateDistances(node, clusters);
                }
            }
        }

        tree = *clusters.begin();
    }

    void updateDistances(Node* target, const set<Node*>& clusters) {
        target->nearestNeighbor = nullptr;
        target->distToNN = 1e9;

        for (Node* node : clusters) {
            if (node == target) continue;

            double d = target->getDistance(node);
            if (d < target->distToNN) {
                target->distToNN = d;
                target->nearestNeighbor = node;
            }
        }
    }

    string getNewick() const {
        return tree->toNewick() + ";";
    }

private:
    void deleteTree(Node* node) {
        if (!node) return;
        deleteTree(node->left);
        deleteTree(node->right);
        delete node;
    }
};

int main() {
    srand(time(0));
    vector<int> taxa;
    for (int i = 0; i < 6; ++i) {
        taxa.push_back(rand() % 100);
    }

    cout << "Initial taxa: ";
    for (int val : taxa) cout << val << " ";
    cout << endl;

    UPGMA upgma(taxa);

    cout << "UPGMA Tree (Newick format): " << endl;
    cout << upgma.getNewick() << endl;

    return 0;
}