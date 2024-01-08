#ifndef REPEAT_SCOUT_H
#define REPEAT_SCOUT_H

#include <cstdint>
#include <vector>
#include <algorithm>
#include <iostream>
#include <sys/types.h>


struct config {
    int64_t max_offset;
    int64_t gap_penalty;
    int64_t cap_penalty;
    int64_t match;
    int64_t mismatch;
    std::size_t max_extend_len;
};

struct genome_token {
    char c = 0;

    [[nodiscard]] int64_t calc_mu(const genome_token& t, const int64_t match, const int64_t mismatch) const {
        if (c == t.c) {
            return match;
        } else {
            return mismatch;
        }
    }


    static std::vector<genome_token> all_tokens;
};

struct repeat_scout {
    repeat_scout(const int64_t exact_repeat_size_, const config& c_, std::vector<genome_token> genome_,
                 std::vector<std::size_t> exact_repeat_pos_) : conf(c_), exact_repeat_pos(std::move(exact_repeat_pos_)),
                                                               genome(std::move(genome_)),
                                                               exact_repeat_size(exact_repeat_size_) {
    }

    void init_right() {right_q = std::vector<genome_token>(genome.begin() + exact_repeat_pos.front(),
                                            genome.begin() + exact_repeat_pos.front() + exact_repeat_size);

        g.resize(exact_repeat_pos.size());
        best_seq_score.resize(exact_repeat_pos.size());
        best_right_border.resize(exact_repeat_pos.size());
        for (std::size_t i = 0; i < exact_repeat_pos.size(); i++) {
            g[i].resize(2);

            best_right_border[i] = exact_repeat_pos[i] + exact_repeat_size;

            g[i][0].resize(2 * conf.max_offset + 1);
            for (int64_t offset = -conf.max_offset; offset <= 0; offset++) {
                g[i][0][offset + conf.max_offset] = -conf.gap_penalty * offset + exact_repeat_size * conf.match;
                best_seq_score[i] = std::max(best_seq_score[i], g[i][0][offset + conf.max_offset]);
            }
            for (int64_t offset = 0; offset <= conf.max_offset; offset++) {
                g[i][0][offset + conf.max_offset] = conf.gap_penalty * offset + exact_repeat_size * conf.match;
                best_seq_score[i] = std::max(best_seq_score[i], g[i][0][offset + conf.max_offset]);
            }
            g[i][1] = g[i][0];
        }

    }

    void extend_right() {
        init_right();
        for (auto y = 0; y < conf.max_extend_len; y++) {
            auto best_new_token = genome_token::all_tokens.front();
            int64_t best_total_score = INT64_MIN;

            for (const auto& new_token: genome_token::all_tokens) {
                int64_t total_score = 0;
                for (auto n = 0; n < exact_repeat_pos.size(); n++) {
                    auto best_n_score = best_seq_score[n] + conf.cap_penalty;
                    for (int64_t offset = -conf.max_offset; offset <= conf.max_offset; offset++) {
                        auto temp_score = calc_score_right(y, n, offset, new_token);
                        best_n_score = std::max(best_n_score, temp_score);
                    }
                    total_score += best_n_score;
                }
                if (total_score > best_total_score) {
                    best_total_score = total_score;
                    best_new_token = new_token;
                }
            }
            right_q.push_back(best_new_token);

            //recalc best_seq_score
            for (size_t n = 0; n < exact_repeat_pos.size(); n++) {
                for (int64_t offset = -conf.max_offset; offset <= conf.max_offset; offset++) {
                    auto temp_score = calc_score_right(y, n, offset, best_new_token);
                    if (temp_score > best_seq_score[n]) {
                        best_seq_score[n] = temp_score;
                        best_right_border[n] = exact_repeat_pos[n] + exact_repeat_size + y - offset;
                    }
                }
            }
            //rcalc g

            for (auto n = 0; n < exact_repeat_pos.size(); n++) {
                for (int64_t offset = -conf.max_offset; offset <= conf.max_offset; offset++) {
                    g[n][y % 2][offset + conf.max_offset] = calc_score_right(y, n, offset, best_new_token);
                }
            }
        }

        if(0) {
            std::cout << "Printing answer...\n";
            for(auto x: right_q) {
                std::cout << x.c;
            }
            std::cout << std::endl;


            std::vector<genome_token> tst(genome.begin()+1711408, genome.begin()+1711408+100);
            for(auto x: tst) {
                std::cout << x.c;
            }
            std::cout << std::endl;

            for (int i = 0; i < exact_repeat_pos.size(); i++) {
                std::cout << i << ": ";
                for (int64_t x = exact_repeat_pos[i]; x <= best_right_border[i]; x++) {
                    std::cout << genome[x].c;
                }
                std::cout << '\n';
                //std::cout << exact_repeat_pos[i] << '\t' << best_right_border[i] << '\n';
            }
        }
    }

    void extend_left() {

        auto const old_best_right_border = best_right_border;
        auto const old_right_q = right_q;

        std::reverse(genome.begin(), genome.end());

        for (auto& pos: exact_repeat_pos) {
            pos = genome.size() - (pos + exact_repeat_size);
        }
        extend_right();

        best_left_border = best_right_border;
        for(auto& pos :best_left_border) {
            pos = genome.size()-pos-1;
        }
        left_q = right_q;
        std::reverse(left_q.begin(), left_q.end());

        right_q = old_right_q;
        best_right_border = old_best_right_border;

        for (auto& pos: exact_repeat_pos) {
            pos = genome.size() - (pos + exact_repeat_size);
        }

        std::reverse(genome.begin(), genome.end());
    }

    [[nodiscard]] int64_t calc_score_right(const std::size_t y, const std::size_t n, const int64_t offset,
                                           const genome_token new_token) const {
        int64_t max_g = g[n][(y + 1) % 2][offset + conf.max_offset] + new_token.calc_mu(
                            genome[exact_repeat_pos[n] + exact_repeat_size + y - offset], conf.match, conf.mismatch);

        for (int64_t k = 0; offset + k <= conf.max_offset; k++) {
            //k in range [1; b - offset]
            max_g = std::max(
                max_g, g[n][(y + 1) % 2][offset + k + conf.max_offset] + k * conf.gap_penalty + new_token.calc_mu(
                           genome[exact_repeat_pos[n] + exact_repeat_size + y - offset], conf.match,
                           conf.mismatch));

            max_g = std::max(max_g, g[n][(y + 1) % 2][offset + k - 1 + conf.max_offset] + (k + 1) * conf.gap_penalty);
        }
        //k = conf.max_offset-offset+1
            max_g = std::max(max_g, g[n][(y + 1) % 2][ conf.max_offset] + ((conf.max_offset - offset +  1) + 1) * conf.gap_penalty);

        return max_g;
    }

    const config conf;
    std::vector<std::size_t> exact_repeat_pos;
    std::vector<genome_token> genome;
    std::vector<std::size_t> best_right_border;
    std::vector<std::size_t> best_left_border;
    std::vector<genome_token> right_q;
    std::vector<genome_token> left_q;
    int64_t exact_repeat_size;
    std::vector<std::vector<std::vector<int64_t>>> g;
    std::vector<int64_t> best_seq_score;
};


#endif //REPEAT_SCOUT_H
