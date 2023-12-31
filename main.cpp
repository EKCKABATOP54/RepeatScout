#include <iostream>
#include "repeat_scout.h"
#include <fstream>

int main() {
    config c;
    c.max_offset = 10;
    c.gap_penalty = -5;
    c.cap_penalty = -20;
    c.match = 1;
    c.mismatch = -1;
    c.max_extend_len = 100;


    std::string genome_file = "/home/androposh/CLionProjects/RepeatScout/humanFormattedUppercase.fa";
    std::ifstream fs(genome_file);

    if (!fs.is_open()) {
        std::cerr << "Unable to open file: " << genome_file << std::endl;
        return 1; // return an error code
    }

    // Read the contents of the file into a string
    std::vector<genome_token> genome((std::istreambuf_iterator<char>(fs)),
                                     (std::istreambuf_iterator<char>()));

    // Close the file
    fs.close();

    std::vector<size_t> poses;

    std::string pos_file = "/home/androposh/CLionProjects/RepeatScout/poses.txt";


    std::ifstream fs2(pos_file);
    if (!fs2.is_open()) {
        std::cerr << "Unable to open file: " << pos_file << std::endl;
        return 1; // return an error code
    }

    size_t x;
    while (fs2 >> x) {
        std::string s;
        poses.push_back(x);
        fs2 >> s;
    }
/*
    for(auto pos: poses) {
        std::vector<genome_token> rp(genome.begin()+pos, genome.begin()+ pos + 10);
        for(auto xx:rp){
            std::cout << xx.c;
        }
        std::cout << std::endl;
    }
    */

//ACCAGCCTGGCC
    repeat_scout rs(9, c, genome, poses);
    rs.extend_right();
}
