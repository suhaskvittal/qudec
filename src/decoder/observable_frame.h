#ifndef DECODER_OBSERVABLE_FRAME_H
#define DECODER_OBSERVABLE_FRAME_H

#include <algorithm>
#include <bit>
#include <cstddef>
#include <concepts>
#include <cstdint>
#include <iterator>
#include <ostream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace dec {

// Sparse parity set of logical observables. Most physical faults cross only a
// few logical supports, even when the code encodes hundreds of qubits.
class ObservableFrame {
public:
    ObservableFrame(uint64_t low_bits = 0) {
        while (low_bits) {
            ids_.push_back(std::countr_zero(low_bits));
            low_bits &= low_bits - 1;
        }
    }
    static ObservableFrame single(uint32_t id) {
        ObservableFrame result;
        result.ids_.push_back(id);
        return result;
    }
    void toggle(uint32_t id) {
        auto it = std::lower_bound(ids_.begin(), ids_.end(), id);
        if (it != ids_.end() && *it == id) ids_.erase(it);
        else ids_.insert(it, id);
    }
    bool contains(uint32_t id) const {
        return std::binary_search(ids_.begin(), ids_.end(), id);
    }
    explicit operator bool() const { return !ids_.empty(); }
    // Checked conversion supports existing one-observable diagnostic callers.
    explicit operator uint64_t() const {
        uint64_t out = 0;
        for (auto id : ids_) {
            if (id >= 64) throw std::out_of_range("Frame does not fit in uint64_t");
            out |= uint64_t{1} << id;
        }
        return out;
    }
    const std::vector<uint32_t>& ids() const { return ids_; }
    ObservableFrame& operator^=(const ObservableFrame& other) {
        std::vector<uint32_t> out;
        std::set_symmetric_difference(ids_.begin(), ids_.end(),
            other.ids_.begin(), other.ids_.end(), std::back_inserter(out));
        ids_ = std::move(out);
        return *this;
    }
    friend ObservableFrame operator^(ObservableFrame lhs, const ObservableFrame& rhs) {
        lhs ^= rhs;
        return lhs;
    }
    friend bool operator==(const ObservableFrame&, const ObservableFrame&) = default;
    template<std::integral I>
    friend bool operator==(const ObservableFrame& lhs, I rhs) {
        return lhs == ObservableFrame(uint64_t(rhs));
    }
    template<std::integral I>
    friend bool operator==(I lhs, const ObservableFrame& rhs) {
        return ObservableFrame(uint64_t(lhs)) == rhs;
    }
    friend bool operator<(const ObservableFrame& lhs, const ObservableFrame& rhs) {
        size_t i = lhs.ids_.size(), j = rhs.ids_.size();
        while (i && j) {
            if (lhs.ids_[i - 1] != rhs.ids_[j - 1])
                return lhs.ids_[i - 1] < rhs.ids_[j - 1];
            --i;
            --j;
        }
        return i == 0 && j != 0;
    }
    friend std::ostream& operator<<(std::ostream& out, const ObservableFrame& frame) {
        out << '{';
        for (auto id : frame.ids_) out << 'L' << id << ' ';
        return out << '}';
    }
private:
    std::vector<uint32_t> ids_;
};

} // namespace dec
#endif
