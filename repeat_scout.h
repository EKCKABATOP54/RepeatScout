#ifndef REPEAT_SCOUT_H
#define REPEAT_SCOUT_H

#include <cstdint>
#include <vector>
#include <algorithm>
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
        right_q = std::vector<genome_token>(genome.begin() + exact_repeat_pos.front(),
                                            genome.begin() + exact_repeat_pos.front() + exact_repeat_size);

        g.resize(exact_repeat_pos.size());
        for (std::size_t i = 0; i < exact_repeat_pos.size(); i++) {
            g[i].resize(2);

            g[i][0].resize(2 * conf.max_offset + 1);
            for (int64_t offset = -conf.max_offset; offset <= 0; offset++) {
                g[i][0][offset + conf.max_offset] = -conf.gap_penalty * offset + exact_repeat_size * conf.match;
                best_seq_score[i] = std::max(best_seq_score[i], g[i][0][offset]);
            }
            for (int64_t offset = 0; offset <= conf.max_offset; offset++) {
                g[i][0][offset + conf.max_offset] = conf.gap_penalty * offset + exact_repeat_size * conf.match;
                best_seq_score[i] = std::max(best_seq_score[i], g[i][0][offset]);
            }
            g[i][1] = g[i][0];
        }
    }

    void extend_right() {
        for (auto y = 0; y < conf.max_extend_len; y++) {
            auto best_new_token = genome_token::all_tokens.front();
            int64_t best_total_score = INT64_MIN;

            for (const auto& new_token: genome_token::all_tokens) {
                int64_t total_score = INT64_MIN;
                for (auto n = 0; n < exact_repeat_pos.size(); n++) {
                    auto best_n_score = best_seq_score[n] + conf.cap_penalty;
                    for (int64_t offset = -conf.max_offset; offset <= conf.max_offset; offset++) {
                        auto temp_score = calc_score_right(y, n, offset, new_token);
                        best_n_score += std::max(best_n_score, temp_score);
                    }
                    total_score += best_n_score;
                }
                if (total_score > best_total_score) {
                    best_total_score = total_score;
                    best_new_token = new_token;
                }
            }
            right_q.push_back(best_new_token);
        }
    }

    void extend_left() {
        std::reverse(right_q.begin(), right_q.end());
        std::reverse(genome.begin(), genome.end());



        std::reverse(right_q.begin(), right_q.end());
        std::reverse(genome.begin(), genome.end());
    }

    [[nodiscard]] int64_t calc_score_right(const std::size_t y, const std::size_t n, const int64_t offset,
                             const genome_token new_token) const  {
        int64_t max_g = g[n][(y - 1) % 2][offset] + new_token.calc_mu(
                            genome[exact_repeat_pos[n] + exact_repeat_size + y - offset], conf.match, conf.mismatch);

        for (int64_t k = 1; offset + k <= conf.max_offset; k++) {
            //k in range [1; b - offset]
            max_g = std::max(max_g, g[n][(y - 1) % 2][offset + k]) - k * conf.gap_penalty + new_token.calc_mu(
                        genome[exact_repeat_pos[n] + exact_repeat_size + y - offset], conf.match,
                        conf.mismatch);

            max_g = std::max(max_g, g[n][(y - 1) % 2][offset + k - 1] - (k + 1) * conf.gap_penalty);
        }
        return max_g;
    }

    const config conf;
    std::vector<std::size_t> exact_repeat_pos;
    std::vector<genome_token> genome;
    std::vector<genome_token> right_q;
    int64_t exact_repeat_size;
    std::vector<std::vector<std::vector<int64_t>>> g;
    std::vector<int64_t> best_seq_score;
};


#endif //REPEAT_SCOUT_H
