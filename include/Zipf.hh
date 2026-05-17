#pragma once  // Avoid include header multiple times

#include <iostream> 
#include <fstream>    // for file handling
#include <string>     // for string handling
#include <utility>    // for std::pair
#include <map>        // for std::map   word -> frequency
#include <vector>     // for std::vector
#include <algorithm>  // for std::sort

// std::string const & w : for passing by reference without copying, and to prevent modification of the argument (more efficient for large objects like std::string)
// std::map<std::string, int> & m : for passing by reference to allow modification of the map 
void add_word(std::string const & w, 
              std::map<std::string, int> & m);

bool is_letter(char c);
char to_lowercase(char c);
void write_csv(const std::vector<std::pair<std::string, int>> & ranking, const std::string & file_name);