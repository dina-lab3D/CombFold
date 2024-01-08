#ifndef BITID_H
#define BITID_H

#include <ostream>
#include <bitset>


template <size_t N>
class CustomBitset : public std::bitset<N> {
public:
    CustomBitset() : std::bitset<N>() {}

    CustomBitset(const std::bitset<N>& b) : std::bitset<N>(b) {}

    // assert val < n and then set the val bit, this is different from std::bitset
    CustomBitset(unsigned long long val) : std::bitset<N>(0) {
        assert(val < N);
        this->set(val);
    }

    // don't print leading zeros
    std::string to_string() const {
        std::string str = std::bitset<N>::to_string();
        size_t first_one = str.find('1');
        return first_one == std::string::npos ? "0" : str.substr(first_one);
    }
};

template <size_t N>
std::ostream& operator<<(std::ostream& os, const CustomBitset<N>& b) {
    os << b.to_string();
    return os;
}


namespace std {
    template <size_t N>
    struct hash<CustomBitset<N>> {
        size_t operator()(const CustomBitset<N>& b) const {
            return std::hash<std::string>()(b.to_string());
        }
    };
}


typedef CustomBitset<128> BitId;

#endif /* BITID_H */