#ifndef _VISUALIZATION_H
#define _VISUALIZATION_H

// OpenAI GPT-6: Accept the decoder's packed syndrome bitset in visualization.
#include "packed_bitset.h"
#include <list>
#include <vector>

#include "common.h"

namespace tesseract_decoder {

struct Visualizer {
  void add_detector_coords(const std::vector<std::vector<double>>&);
  void add_errors(const std::vector<common::Error>&);
  void add_activated_errors(int64_t node_idx, const std::vector<common::ErrorChainNode>& arena);
  void add_activated_detectors(const PackedBitset&, size_t);

  void write(const char* fpath);

 private:
  std::list<std::string> lines;
};

}  // namespace tesseract_decoder

#endif
