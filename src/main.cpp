#include "Zipf.hh"

int main() {

  std::map<std::string,int> words;

  add_word( "hello", words );
  add_word( "hello", words );
  add_word( "world", words );

  for ( auto const & item : words ) {
    std::cout << item.first << ": " << item.second << '\n';
  }

  return 0;
}