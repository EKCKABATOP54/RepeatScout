
#include <vector>
#include "repeat_scout.h"

std::map<std::pair<genome_token, genome_token>, int64_t> genome_token::cache{};
std::vector<genome_token> genome_token::all_tokens_{};