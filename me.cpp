#include <iostream>
#include <vector>
#include <limits>
#include <cmath>
#include <map>

using namespace std;

pair<int, int> findMinPair(const vector<vector<double>>& distMatrix) {//minimum element in the distance matrix
    int n = distMatrix.size();
    double minDist = numeric_limits<double>::max();
    pair<int, int> minPair;
    
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (distMatrix[i][j] < minDist) {
                minDist = distMatrix[i][j];
                minPair = {i, j};
            }
        }
    }
    return minPair;
}

typedef vector<vector<double>> Matrix;
Matrix updateDistanceMatrix(const Matrix& distMatrix, int u, int v) {//new distance matrix after merging nodes
    int n = distMatrix.size();
    Matrix newMatrix(n - 1, vector<double>(n - 1, 0.0));
    
    int newIdx = 0;
    for (int i = 0; i < n; i++) {
        if (i == u || i == v) continue;
        int newJ = 0;
        for (int j = 0; j < n; j++) {
            if (j == u || j == v) continue;
            newMatrix[newIdx][newJ] = distMatrix[i][j];
            newJ++;
        }
        newIdx++;
    }
    return newMatrix;
}

void constructMETree(Matrix distMatrix) {//construct the Minimum Evolution tree using NJ method
    int n = distMatrix.size();
    map<int, string> labels;
    for (int i = 0; i < n; i++) labels[i] = "Node" + to_string(i);
    
    while (distMatrix.size() > 1) {
        pair<int, int> minPair = findMinPair(distMatrix);
        int u = minPair.first, v = minPair.second;
        
        cout << "Merging " << labels[u] << " and " << labels[v] << endl;
        labels[u] = "(" + labels[u] + ", " + labels[v] + ")";
        
        distMatrix = updateDistanceMatrix(distMatrix, u, v);
    }
    
    cout << "Final Tree: " << labels[0] << endl;
}

int main() {
    //example distance matrix (symmetric)
    Matrix distMatrix = {
        {0, 5, 9, 9, 8},
        {5, 0, 10, 10, 9},
        {9, 10, 0, 8, 7},
        {9, 10, 8, 0, 3},
        {8, 9, 7, 3, 0}
    };
    
    constructMETree(distMatrix);
    return 0;
}
