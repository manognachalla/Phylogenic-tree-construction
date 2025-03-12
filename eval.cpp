#include "tree.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <unordered_map>
#include <sstream>
#include <stack>  // ✅ FIX: Include stack for Newick parsing
#include <algorithm>  // ✅ FIX: Include algorithm for std::remove
#include <cmath>  // ✅ FIX: Include cmath for isnan()

using namespace std;

// Function to parse a Newick tree and store its structure
void parseNewick(const std::string &tree, std::unordered_map<std::string, int> &cladeCounts) {
    std::stack<std::string> cladeStack;
    std::string clade;
    
    for (char c : tree) {
        if (c == '(') {
            cladeStack.push("");
        } else if (c == ',') {
            if (!cladeStack.empty()) {
                cladeStack.top() += clade + ",";
                clade.clear();
            }
        } else if (c == ')') {
            if (!cladeStack.empty()) {
                cladeStack.top() += clade;
                std::string completeClade = cladeStack.top();
                cladeStack.pop();
                if (!cladeStack.empty()) {
                    cladeStack.top() += completeClade;
                }
                cladeCounts[completeClade]++;
                clade.clear();
            }
        } else {
            clade += c;
        }
    }
}



// Function to compute bootstrap confidence scores
void computeBootstrapSupport(const vector<string> &bootstrapTrees, int numBootstraps) {
    unordered_map<string, int> cladeCounts;

    for (const string &tree : bootstrapTrees) {
        unordered_map<string, int> singleTreeCounts;
        parseNewick(tree, singleTreeCounts);
        for (const auto &entry : singleTreeCounts) {
            cladeCounts[entry.first]++;
        }
    }

    cout << "\nBootstrap Confidence Scores:\n";
    for (const auto &entry : cladeCounts) {
        double confidence = (entry.second / (double)numBootstraps) * 100;
        cout << "Clade: " << entry.first << " - Support: " << confidence << "%\n";
    }
}



// Function to check if a nucleotide is a purine (A or G)
bool isPurine(char base) {
    return base == 'A' || base == 'G';
}

// Function to count transitions and transversions between two sequences
pair<int, int> countTransitionsTransversions(const string &seq1, const string &seq2) {
    int transitions = 0, transversions = 0;
    for (size_t i = 0; i < seq1.length(); i++) {
        char b1 = seq1[i], b2 = seq2[i];
        if (b1 != b2) {
            if ((isPurine(b1) && isPurine(b2)) || (!isPurine(b1) && !isPurine(b2))) {
                transitions++;
            } else {
                transversions++;
            }
        }
    }
    return make_pair(transitions, transversions);
}

// Function to compute transition/transversion ratio for all sequence pairs
void computeTransitionTransversionRatio(const vector<string> &names, const vector<string> &sequences) {
    cout << "Transition/Transversion Ratios:\n";
    for (size_t i = 0; i < sequences.size(); i++) {
        for (size_t j = i + 1; j < sequences.size(); j++) {
            pair<int, int> result = countTransitionsTransversions(sequences[i], sequences[j]);
            int transitions = result.first;
            int transversions = result.second;
            double ratio = (transversions == 0) ? 0.0 : (double)transitions / transversions;
            cout << names[i] << " vs " << names[j] << ": " << ratio << endl;
        }
    }
}

vector<vector<string>> bootstrapSequences(const vector<string> &sequences, int numBootstrap) {
    vector<vector<string>> resampled;  // ✅ Now it matches expected return type
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(0, sequences[0].size() - 1);

    for (int i = 0; i < numBootstrap; i++) {
        vector<string> resampledSet;  // ✅ Each bootstrap set should be a vector<string>
        for (const auto &seq : sequences) {
            string resampledSeq;
            for (size_t j = 0; j < seq.size(); j++) {
                int randIndex = dis(gen);
                resampledSeq += seq[randIndex];
            }
            resampledSet.push_back(resampledSeq);  // ✅ Store resampled sequence
        }
        resampled.push_back(resampledSet);  // ✅ Store the entire resampled set
    }
    return resampled;
}





// Function to process bootstrap analysis and save to a file
void performBootstrapAnalysis(const vector<string> &sequences, int numBootstrap, const string &outputFile) {
    vector<vector<string>> bootstrappedSets = bootstrapSequences(sequences, numBootstrap);

    ofstream bootstrapFile(outputFile);
    for (size_t i = 0; i < bootstrappedSets.size(); i++) {
        for (size_t j = 0; j < bootstrappedSets[i].size(); j++) {
            bootstrapFile << ">Bootstrap_" << i + 1 << "_Seq" << j + 1 << "\n" << bootstrappedSets[i][j] << "\n";
        }
    }
    bootstrapFile.close();

    cout << "Bootstrap resampling completed. Output saved to " << outputFile << endl;
}

