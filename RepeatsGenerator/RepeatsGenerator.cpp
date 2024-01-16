#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <list>
#include <random>
#include <set>

std::string gen_rand_exact_repeat(int64_t repeat_size) {
  std::random_device dev;
  std::mt19937 rng(dev());
  std::uniform_int_distribution<std::mt19937::result_type> dist(0, 3);
  char nucl[4] = {'A', 'C', 'T', 'G'};
  std::vector<char> exact_repeat(repeat_size);
  for (int64_t i = 0; i < repeat_size; i++) {
    exact_repeat[i] = nucl[dist(rng)];
  }
  return std::string(exact_repeat.begin(), exact_repeat.end());
}

std::string mutilate_exact_repeat(const std::string &s, int64_t max_gaps,
                                  int64_t max_mismatch) {
  std::random_device dev;
  std::mt19937_64 rng(dev());
  std::uniform_int_distribution<std::mt19937_64::result_type> dist_gaps(
      0, 2*max_gaps);
  std::uniform_int_distribution<std::mt19937_64::result_type> dist_mismatch(
      0, max_mismatch);
  std::uniform_int_distribution<std::mt19937_64::result_type> dist_nucl(0, 3);
  std::uniform_int_distribution<std::mt19937_64::result_type> dist_general(
      0, std::numeric_limits<int64_t>::max()-1);
  char nucl[4] = {'A', 'C', 'T', 'G'};
  const int64_t offset = dist_gaps(rng) - max_gaps;
  std::vector<char> genome(s.begin(), s.end());

  std::vector<int64_t> remove_add_pos;
  {
    std::set<int64_t> st;
    int64_t abs_offset = offset;
    if (abs_offset < 0) {
      abs_offset *= -1;
    }
    for (int64_t i = 0; i < abs_offset; i++) {
      int64_t new_pos = dist_general(rng)%s.size();
      while (st.contains(new_pos)) {
        new_pos = dist_general(rng)%s.size();
      }
      st.insert(new_pos);
    }
    remove_add_pos = std::vector(st.begin(), st.end());
  }

  if (offset < 0) {
    for (int64_t i = 0; i < -offset; i++) {
      genome[remove_add_pos[i]] = 'N';
    }
    genome.erase(std::remove(genome.begin(), genome.end(), 'N'), genome.end());
  } else {
    for (size_t i = 0; i < remove_add_pos.size(); i++) {
      remove_add_pos[i] += i; // position of new symbols after insertion
    }
    genome.resize(s.size() + offset);
    int64_t inserted_nucleotids = 0;
    for (int64_t x = s.size() + offset - 1; x >= 0; x--) {
      if (!remove_add_pos.empty() && x == remove_add_pos.back()) {
        genome[x] = nucl[dist_nucl(rng)];
        inserted_nucleotids += 1;
        remove_add_pos.pop_back();
      } else {
        genome[x] = genome[x - offset + inserted_nucleotids];
      }
    }
  }

  return std::string(genome.begin(), genome.end());
}

int main(int argc, char *argv[]) {
  if (argc != 7) {
    std::cerr << "Usage: " << argv[0]
              << " <repeat_size> <max_gaps> <repeats_cnt> <max_mismatch> "
                 "<non_repetitive_cnt> <file>"
              << std::endl;
    return 1;
  }
  int64_t repeat_size = std::atoll(argv[1]);
  int64_t max_gaps = std::atoll(argv[2]);
  int64_t repeats_cnt = std::atoll(argv[3]);
  int64_t max_mismatch = std::atoll(argv[4]);
  int64_t non_repetitive_cnt = std::atoll(argv[5]);
  const std::string filename(argv[6]);
  std::ofstream fs(filename);
  if (!fs) {
    std::cout << "Can't open file";
    return 1;
  }

  std::random_device dev;
  std::mt19937_64 rng(dev());
  std::uniform_int_distribution<std::mt19937_64::result_type> dist;

  std::vector<int64_t> rep_rand_poses;
  for (int64_t i = 0; i < repeats_cnt; i++) {
    rep_rand_poses.push_back(dist(rng) % non_repetitive_cnt);
  }
  std::sort(rep_rand_poses.begin(), rep_rand_poses.end());

  std::string exact_repeat = gen_rand_exact_repeat(repeat_size);

  int64_t rep_ind = 0;
  std::string result;
  for (int64_t gpos = 0; gpos < non_repetitive_cnt;) {
    if (rep_ind!=rep_rand_poses.size() &&  gpos == rep_rand_poses[rep_ind]) {
      result += mutilate_exact_repeat(exact_repeat, max_gaps, max_mismatch);
      rep_ind += 1;
    } else {
      char nucl[4] = {'A', 'C', 'T', 'G'};
      std::uniform_int_distribution<std::mt19937_64::result_type> dist_nucl(0,
                                                                            3);
      result += std::string(1, nucl[dist_nucl(rng)]);
      gpos += 1;
    }
  }

  fs << result;
  std::cout << "Exact repeat: " << exact_repeat;
}