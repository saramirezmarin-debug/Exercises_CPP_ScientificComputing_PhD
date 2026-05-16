#include "Zipf.hh"

void add_word(std::string const & w, std::map<std::string, int> & m)
{
    auto res = m.find(w);

    if( res != m.end() )
    {
        ++res->second; // word already exists, increment frequency
    }
    else
    {
        m[w] = 1; // word does not exist, add it with frequency 1
    }
}