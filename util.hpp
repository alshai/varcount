#ifndef UTIL_HPP
#define UTIL_HPP

#include <vector>
#include <cstdio>
#include <cstring>
#include <string>


std::vector<std::string> parse_ids(const char* s) {
    std::vector<std::string> ss;
#if TAHER
    char* s__ = (char*) malloc(sizeof(char) * (strlen(s) + 1));
    strcpy(s__, s);
    char* s_ = std::strtok(s__, ";");
    while (s_ != NULL) {
        ss.push_back(s_);
        s_ = std::strtok(NULL, ";");
    }
    free(s__);
#else
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
#endif
    return ss;
}

#endif
