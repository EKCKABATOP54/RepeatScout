#include <iostream>
#include "repeat_scout.h"
#include <fstream>
#include <numeric>
#include <cassert>
#include <chrono>

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
    std::vector<bool> delims(std::accumulate(genome_tokens.begin(), genome_tokens.end(), static_cast<size_t>(0), [](const size_t res, const genome_token&t2) {
        return res + t2.size();
    }), false);


    size_t tokens_total_size = 0;
    for(const auto&tok :genome_tokens) {
        delims[tokens_total_size]=true;
        tokens_total_size+=tok.size();
        whole_genome += tok.to_string();
    }
    assert(tokens_total_size==whole_genome.size());
    std::vector<size_t> positions;
    size_t pos = whole_genome.find(repeat, 0);
    while (pos != std::string::npos) {
        positions.push_back(pos);
        delims[pos] = true;
        delims[pos+repeat.size()] = true;
        pos = whole_genome.find(repeat, pos + 1);
    }
    genome_tokens.clear();
    repeat_poses_tok.clear();
    std::string token;
    token.push_back(whole_genome[0]);
    size_t rep_ind=0;
    for(size_t i=1;i<whole_genome.size();i++) {
        if(delims[i]) {
            genome_tokens.emplace_back(token);
            token.clear();
            if(rep_ind < positions.size() && i == positions[rep_ind]) {
                rep_ind+=1;
                repeat_poses_tok.push_back(genome_tokens.size());
            }
        }
        token.push_back(whole_genome[i]);

    }
    positions.shrink_to_fit();
    genome_tokens.shrink_to_fit();
    genome_token::read_all_tokens(genome_tokens);
}


int main(int argc, char* argv[]) {///tmp/tmp.gM43pQaPIS/cmake-build-debug-remote-host/bioinformatics
    std::cout << argv[0];
    std::ofstream outFile("aboba4_left_right_small.txt");
    auto start_time = std::chrono::high_resolution_clock::now();
    std::cout.rdbuf(outFile.rdbuf());
    //std::iostream::sync_with_stdio(false);
    config c;
    c.max_offset = 5;
    c.gap_penalty = -5;
    c.cap_penalty = -10;
    c.match = 1;
    c.mismatch = -1;
    c.max_extend_len = 50;

    std::cout << "Starting..."  << std::endl;

    std::vector<genome_token> genome;
    //bool genome_read_res = read_genome_from_file("/home/androposh/CLionProjects/RepeatScout/humanFormattedUppercase.fa", genome);
    //bool genome_read_res = read_genome_from_file("/home/androposh/CLionProjects/RepeatScout/GCF_000002985.6_WBcel235_genomic.fna_modified", genome);
    //const std::string source_file = "/bioserver/aposhtarenko/hm/tokenizer/tokenized.txt";
    const std::string source_file = "/home/aposhtarenko/hm/tokenizer/tokenized.txt";
    bool genome_read_res = read_genome_from_file(source_file, genome);
    if(!genome_read_res) {
        std::cout << "Error";
        return 0;
    }

    std::vector<size_t> poses;


    //bool read_pos_from_file_res = read_repeat_pos_from_file("/home/androposh/CLionProjects/RepeatScout/poses.txt", poses);
    /*
    bool read_pos_from_file_res = read_repeat_pos_from_file("/home/aposhtarenko/hm/tests/poses.txt", poses);
    if(!read_pos_from_file_res) {
        std::cout << "Error";
        return 0;
    }
    poses.resize(poses.size(), 2);*/


    std::cout << "File read" << std::endl;
    build_repeat_pos("AAATTCAAGATAAAT", genome, poses);
    std::cout << "Position built" << std::endl;

    //genome_token::read_all_tokens_from_file("/bioserver/tokenizer/tokenized.txt");r
    //genome_token::read_all_tokens_from_file(source_file);

    //std::cout << "Tokens read" << std::endl;


//AAATTCAAGATAAAT
    repeat_scout rs(15, c, genome, poses);
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
    std::cout << "Printing left part of Q and seqs\n";
    //std::cout << std::string(2*(max_left_extention_len-(rs.left_q.size()-rs.exact_repeat_size)), '\t');
    //for(auto e:total_q) {
    for(const auto& e:rs.left_q) {
        std::cout << e.to_string() << ",";
    }
    std::cout << '\n';
    /*
    //for(auto p=rs.best_left_border[n];p<=rs.best_right_border[n];p++) {
            std::cout << rs.genome[p].to_string() << ",";
        }
        std::cout << '\n';
*/

    for(int n=0;n<rs.exact_repeat_pos.size();n++) {
        //std::size_t space_cnt = max_left_extention_len - (rs.exact_repeat_pos[n]-rs.best_left_border[n]);
        //std::cout << std::string(2*space_cnt, '\t');


        for(auto p= rs.best_left_border[n];p<=rs.exact_repeat_pos[n]+rs.exact_repeat_size;p++) {
            std::cout << rs.genome[p].to_string() << ",";
        }
        std::cout << '\n';
    }

    std::cout << "Printing right part of Q and seqs\n";

    for(const auto& e:rs.right_q) {
        std::cout << e.to_string() << ",";
    }
    std::cout << '\n';
    for(int n=0;n<rs.exact_repeat_pos.size();n++) {
        //std::size_t space_cnt = max_left_extention_len - (rs.exact_repeat_pos[n]-rs.best_left_border[n]);
        //std::cout << std::string(2*space_cnt, '\t');

        //for(auto p=rs.best_left_border[n];p<=rs.best_right_border[n];p++) {
        for(auto p= rs.exact_repeat_pos[n];p<=rs.best_right_border[n];p++) {
            std::cout << rs.genome[p].to_string() << ",";
        }
        std::cout << '\n';
    }


    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();

    std::cout << "Execution time: " << duration_seconds << " seconds." << std::endl;


}
