#include "datf_middle.hpp"
#include <iostream>
using namespace gdl;
int main() {
    MiddleTables tables;
    int rounds;
    while (std::cin >> rounds) {
        State delta{}, gamma{};
        uint64_t rc;
        for (auto &w : delta) std::cin >> std::hex >> w;
        for (auto &w : gamma) std::cin >> std::hex >> w;
        std::cin >> std::hex >> rc >> std::dec;
        if (!std::cin) return 2;
        auto c = rounds == 1 ? middle_one(delta, gamma) : MiddleTwo(delta, tables).evaluate(gamma);
        if (__builtin_parityll(c.rc_mask & rc)) c.sign = -c.sign;
        std::cout << c.sign << ' ' << c.weight << ' ' << std::hex << c.rc_mask << std::dec << '\n';
    }
}
