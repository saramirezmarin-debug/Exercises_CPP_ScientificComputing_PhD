#include "Zipf.hh"

// ./exe/Zipf.exe test.txt
// argc is the argument count  -> 2
// argv is the argument vector -> argv[0] = "./exe/Zipf.exe", argv[1] = "test.txt"
int main(int argc, char *argv[])
{

    /*
    std::map<std::string,int> words; // map to store word frequencies

    add_word( "hello", words );
    add_word( "hello", words );
    add_word( "world", words );

    for ( auto const & item : words ) {
      std::cout << item.first << ": " << item.second << '\n';
    }
    */

    if (argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " <file>\n";
        return 1;
    }

    std::string file_name = argv[1]; // Copy the file name argument into a string variable.


    std::ifstream file(file_name);     // Open the file for reading
    std::string line; // Variable to store line read from the file

    // Check if the file was opened successfully
    if (!file)
    {
        std::cerr << "Error: Could not open file " << file_name << "\n";
        return 1;
    }

    // Read the file line by line
    while (std::getline(file, line))
    {
        std::cout << line << "\n";
        // update words
        std::string word; // Variable to build the current word

        for (auto c : line)
        {
            if (is_letter(c)) // If the character is a letter, add it to the current word
            {
                word += c; 
            }
            else              // If not is a letter, is a separator
            {
                if (!word.empty())
                {
                   std::cout << word << '\n';
                   word.clear();
                }
            }
        }

        // After the loop, check if there is a last word that was being built

        if (!word.empty())
        {
            std::cout << word << '\n';
        }
    }

    return 0;
}