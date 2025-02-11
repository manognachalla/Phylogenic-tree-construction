#include "tree.h"
#include <iostream>
#include <fstream>
#include <algorithm>

sequence read_fasta(std::string filename) {
    sequence sequence_list;
    std::ifstream input(filename.c_str());
    
    if (!input.good()) {
        std::cerr << "Error opening '" << filename << std::endl;
    }
    else if (input.good()) {
        std::cout << "Reading '" << filename << "'...";
    }
    
    std::string line, content, name;
    while (std::getline(input, line).good()) {
        if (line.length() <= 1) {
            if (!content.empty()) {
                sequence_list.seq.push_back(content);
                sequence_list.name.push_back(name);
            }
            content.clear();
            name.clear();
        }
        if (line[0] == '>') {
            line = line.substr(1, line.length() -1);
            name = line;
        }
        else {
            content += line;
        }
    }
    input.close();
    std::cout << "  Number of sequences: " << sequence_list.seq.size() << std::endl;
    return sequence_list;
}

std::vector<sequence> read_paml(std::string filename, int n_replicates) {
    std::vector<sequence> sequence_list(n_replicates);
    std::ifstream input;
    input.open(filename);
    int sequence_number = 0;
    int batch = 0;
    int n_seq;
    std::string line, name;

    if (!input.good()) {
        std::cerr << "Error opening '" << filename << std::endl;
    }
    else if (input.good()) {
        std::cout << "\nReading '" << filename << "'...";
    }

    while (std::getline(input, line).good() && (batch < n_replicates)) {
        if (!line.empty() && line.size() > 2 && line.size() < 100) {
            n_seq = std::stoi(line.substr(0,7));
        }

        if (!line.empty() && line[0] != ' ' && line.length() > 100) {
            name = line.substr(0, 10);
            line.erase(line.begin(), line.begin() + 30);
            name.erase(std::remove_if(name.begin(), name.end(), static_cast<int(*)(int)>(&std::isspace)), name.end());
            line.erase(std::remove_if(line.begin(), line.end(), static_cast<int(*)(int)>(&std::isspace)), line.end());
            sequence_list[batch].seq.push_back(line);
            sequence_list[batch].name.push_back(name);
            sequence_number += 1;
            
            if (sequence_number >= n_seq) {
                batch += 1;
                sequence_number = 0;
            }
        }
    }

    std::cout << "  Number of unique sequences: " << n_seq << std::endl;
    input.close();
    return sequence_list;
}

void write_to_file(std::string filename, std::vector<std::string> to_write) {
    std::ofstream output(filename);
    for (const auto& line : to_write) {
        output << line << std::endl;
    }
    output.close();
}
