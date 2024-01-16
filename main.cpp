#include <iostream>
#include "repeat_scout.h"
#include <fstream>

bool read_genome_from_file(const std::string& genome_file, std::vector<genome_token>& genome) { //file format - one token in line
    std::ifstream fs(genome_file);

    if (!fs.is_open()) {
        std::cerr << "Unable to open file: " << genome_file << std::endl;
        return false;
    }

    std::string line;
    while (std::getline(fs, line)) {
        genome.emplace_back(line);
    }


    //genome = std::vector<genome_token> ((std::istreambuf_iterator<char>(fs)),
    //                                 (std::istreambuf_iterator<char>()));

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


void build_repeat_pos(const std::string& repeat, std::vector<genome_token>& genome_tokens, std::vector<size_t> &repeat_poses_tok) { //can split some tokens, to align repeat with beggining of token
    //WARNING: not efficient version
    std::string whole_genome;
    for(const auto&tok :genome_tokens) {
        whole_genome += tok.to_string();
    }
    std::vector<size_t> positions;

    size_t pos = whole_genome.find(repeat, 0);
    while (pos != std::string::npos) {
        positions.push_back(pos);
        pos = whole_genome.find(repeat, pos + 1);
    }
    size_t tokens_total_size = 0;
    size_t current_repeat_ind = 0;
    std::vector<genome_token> genome_tokens_res;
    for(size_t token_ind =0 ;token_ind < genome_tokens.size();) {
        if(tokens_total_size + genome_tokens[token_ind].size() < positions[current_repeat_ind]) {
            genome_tokens_res.push_back(genome_tokens[token_ind]);
            token_ind++;
            tokens_total_size+=genome_tokens[token_ind].size();
        }
        else if(tokens_total_size== positions[current_repeat_ind]) {
            auto [tok1, tok2 ] = genome_tokens[token_ind].split(repeat.size());

            repeat_poses_tok.push_back(token_ind);
            genome_tokens.push_back(tok1);
            current_repeat_ind++;
            tokens_total_size+=tok1.size();

            if(tok2.size()==0) {
                token_ind++;
           }
            else {
                genome_tokens[token_ind] = tok2;
            }

        }
        else {
            auto [tok1, tok2 ] = genome_tokens[token_ind].split(positions[current_repeat_ind]-tokens_total_size);
            genome_tokens_res.push_back(tok1);
            tokens_total_size+=tok1.size();
        }
    }

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
        std::cout << "Error";
        return 0;
    }

    std::vector<size_t> poses;


    bool read_pos_from_file_res = read_repeat_pos_from_file("/home/androposh/CLionProjects/RepeatScout/poses.txt", poses);
    if(!read_pos_from_file_res) {
        std::cout << "Error";
        return 0;
    }

    build_repeat_pos("f", genome, poses);


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
        std::cout << e.to_string();
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
            std::cout << rs.genome[p].to_string();
        }
        std::cout << '\n';
    }

}
