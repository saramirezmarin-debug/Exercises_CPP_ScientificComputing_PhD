#include "Zipf.hh"


// Add a word into the frequency map.
// If the word already exists, increment its frequency.
// Otherwise insert the word with frequency 1.
void add_word(std::string const & w, std::map<std::string, int> & m)
{
    auto res = m.find(w); // find the word in the map

    if (res != m.end())
    {
        ++res->second;    // if found, increment the frequency
    }
    else
    {
        m[w] = 1;         // if not found, add the word with frequency 1
    }

     std::cout << "TEST\n";

}

