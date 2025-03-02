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

void write_to_file(std::string filename, std::vector<std::string> to_write) {
    std::ofstream output(filename);
    for (const auto& line : to_write) {
        output << line << std::endl;
    }
    output.close();
}
