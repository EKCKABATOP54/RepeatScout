#include <iostream>
#include "repeat_scout.h"
#include <fstream>

bool read_genome_from_file(const std::string& genome_file, std::vector<genome_token>& genome) {
    std::ifstream fs(genome_file);

    if (!fs.is_open()) {
        std::cerr << "Unable to open file: " << genome_file << std::endl;
        return false;
    }

    genome = std::vector<genome_token> ((std::istreambuf_iterator<char>(fs)),
                                     (std::istreambuf_iterator<char>()));

    return true;
}

bool read_repeat_pos_from_file(const std::string& pos_file, std::vector<size_t>& poses) {
    std::ifstream fs2(pos_file);
    if (!fs2.is_open()) {
        std::cerr << "Unable to open file: " << pos_file << std::endl;
        return false; // return an error code
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

    return true;

}

int main() {
    std::ofstream outFile("aboba.txt");
    //std::cout.rdbuf(outFile.rdbuf());
    std::iostream::sync_with_stdio(false);
    config c;
    c.max_offset = 5;
    c.gap_penalty = -5;
    c.cap_penalty = -10;
    c.match = 1;
    c.mismatch = -2;
    c.max_extend_len = 5000;

    std::vector<genome_token> genome;
    bool genome_read_res = read_genome_from_file("/home/androposh/CLionProjects/RepeatScout/humanFormattedUppercase.fa", genome);
    if(!genome_read_res) {
        std::cout << "Eror";
        return 0;
    }

    std::vector<size_t> poses;

    bool read_pos_from_file_res = read_repeat_pos_from_file("/home/androposh/CLionProjects/RepeatScout/poses.txt", poses);
    if(!read_pos_from_file_res) {
        std::cout << "Eror";
        return 0;
    }

//ACCAGCCTGGCC
    repeat_scout rs(12, c, genome, poses);
    rs.extend_right();
    rs.extend_left();


    std::size_t max_left_extention_len = rs.left_q.size()-rs.exact_repeat_size;
    for(int n=0;n<rs.exact_repeat_pos.size();n++) {
        std:size_t n_left_len = rs.exact_repeat_pos[n]-rs.best_left_border[n];
        max_left_extention_len = std::max(max_left_extention_len, n_left_len);
    }

    auto total_q = rs.left_q;
    total_q.insert(total_q.end(), rs.right_q.begin()+rs.exact_repeat_size, rs.right_q.end());

    std::cout <<  "Size of Q: " << total_q.size() << '\n';
    std::cout << "Printing Q and seqs\n";
    std::cout << std::string(max_left_extention_len-(rs.left_q.size()-rs.exact_repeat_size), ' ');
    for(auto e:total_q) {
        std::cout << e.c;
    }
    std::cout << '\n';

    /*
    for(int n=0;n<rs.exact_repeat_pos.size();n++) {
        for(auto p=rs.exact_repeat_pos[n];p<=rs.best_right_border[n];p++) {
            std::cout << rs.genome[p].c;
        }
        std::cout << '\n';
    }
*/

    for(int n=0;n<rs.exact_repeat_pos.size();n++) {
        std::size_t space_cnt = max_left_extention_len - (rs.exact_repeat_pos[n]-rs.best_left_border[n]);
        std::cout << std::string(space_cnt, ' ');

        for(auto p=rs.best_left_border[n];p<=rs.best_right_border[n];p++) {
            std::cout << rs.genome[p].c;
        }
        std::cout << '\n';
    }

}
