#ifndef TESSERACT_PACKED_BITSET_H
#define TESSERACT_PACKED_BITSET_H

// Author: OpenAI GPT-6-Luna.
// Purpose: Replace the dynamic bitset operations used by this vendored
// Tesseract snapshot with a compact std::vector<uint64_t> representation.

#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace tesseract_decoder {

class PackedBitset {
 public:
  static constexpr size_t npos = std::numeric_limits<size_t>::max();

  PackedBitset() = default;
  explicit PackedBitset(size_t bits, bool value = false)
      : num_bits_(bits), words_((bits + 63) / 64, value ? ~uint64_t{0} : 0) {
    mask_tail();
  }

  class Reference {
   public:
    Reference(PackedBitset& owner, size_t bit) : owner_(owner), bit_(bit) {}
    Reference& operator=(bool value) {
      owner_.set(bit_, value);
      return *this;
    }
    Reference& operator=(const Reference& value) { return *this = bool(value); }
    operator bool() const { return owner_.test(bit_); }
   private:
    PackedBitset& owner_;
    size_t bit_;
  };

  bool operator[](size_t bit) const { return test(bit); }
  Reference operator[](size_t bit) { return Reference(*this, bit); }
  bool test(size_t bit) const { return (words_[bit / 64] >> (bit % 64)) & 1; }
  void set(size_t bit, bool value) {
    uint64_t mask = uint64_t{1} << (bit % 64);
    if (value) words_[bit / 64] |= mask;
    else words_[bit / 64] &= ~mask;
  }
  PackedBitset& operator|=(const PackedBitset& rhs) {
    for (size_t i = 0; i < words_.size(); ++i) words_[i] |= rhs.words_[i];
    mask_tail();
    return *this;
  }
  PackedBitset& operator&=(const PackedBitset& rhs) {
    for (size_t i = 0; i < words_.size(); ++i) words_[i] &= rhs.words_[i];
    mask_tail();
    return *this;
  }
  PackedBitset operator~() const {
    PackedBitset result(*this);
    for (auto& word : result.words_) word = ~word;
    result.mask_tail();
    return result;
  }
  size_t find_first() const { return find_from(0); }
  size_t find_next(size_t bit) const { return bit == npos ? npos : find_from(bit + 1); }
  bool operator==(const PackedBitset&) const = default;
  size_t hash() const {
    size_t h = static_cast<size_t>(1469598103934665603ULL);
    for (uint64_t word : words_) {
      h ^= static_cast<size_t>(word);
      h *= static_cast<size_t>(1099511628211ULL);
    }
    h ^= num_bits_;
    return h;
  }

 private:
  size_t find_from(size_t bit) const {
    if (bit >= num_bits_) return npos;
    size_t wi = bit / 64;
    uint64_t word = words_[wi] & (~uint64_t{0} << (bit % 64));
    while (true) {
      if (word) return wi * 64 + __builtin_ctzll(word);
      if (++wi == words_.size()) return npos;
      word = words_[wi];
    }
  }
  void mask_tail() {
    if (num_bits_ % 64 && !words_.empty())
      words_.back() &= (uint64_t{1} << (num_bits_ % 64)) - 1;
  }

  size_t num_bits_ = 0;
  std::vector<uint64_t> words_;
};

struct PackedBitsetHash {
  size_t operator()(const PackedBitset& bits) const { return bits.hash(); }
};

}  // namespace tesseract_decoder

#endif
