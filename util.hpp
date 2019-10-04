#ifndef UTIL_HPP
#define UTIL_HPP

#include <vector>
#include <cstdio>
#include <cstring>
#include <string>

namespace util {
template<char SEP=';'>
std::vector<std::string> parse_dsv(const char* s) {
    std::vector<std::string> ss;
    const char *p1 = s, *p2 = std::strchr(p1, ';');
    for(;;) {
        ss.emplace_back(p1, p2 - p1);
        p1 = p2 + 1;
        if((p2 = std::strchr(p1, ';')) == nullptr) {
            p2 = std::strchr(p1, '\0');
            ss.emplace_back(p1, p2 - p1);
            break;
        }
    }
    return ss;
}
std::vector<std::string> parse_ids(const char* s) {
    return parse_dsv<';'>(s);
}
} // util

#endif
