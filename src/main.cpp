#include "Zipf.hh"

// ./exe/Zipf.exe test.txt
// argc is the argument count  -> 2
// argv is the argument vector -> argv[0] = "./exe/Zipf.exe", argv[1] = "test.txt"
int main(int argc, char *argv[])
{
    std::cout << "Program started\n";
    
    std::map<std::string,int> words; // map to store word frequencies
    std::vector<std::pair<std::string, int>> ranking; // vector to store word frequencies for ranking

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
        std::cout << "[DEBUG] New line read: \"" << line << "\"\n";
        // update words
        std::string word; // Variable to build the current word

        for (auto c : line)
        {
            if (is_letter(c)) // If the character is a letter, add it to the current word
            {
                //std::cout << c <<" -> is letter\n";
                word += to_lowercase(c);
            }
            else              // If not is a letter, is a separator
            {
                //std::cout << c << " -> separator\n";
                if (!word.empty())
                {
                   add_word(word, words); // Add the word to the frequency map
                   word.clear();
                }
            }
        }

        // After the loop, check if there is a last word that was being built

        if (!word.empty())
        {
            //std::cout << "Last word: " << word << '\n';
            add_word(word, words); // Add the last word to the frequency map
        }
    }

    // Copy map content into vector.
    for (auto const & item : words)
    {
        ranking.push_back(item);   // similar to append()
    }

    // begin, end, comparator (lambda)
    std::sort(ranking.begin(), ranking.end(), [](auto &left, auto &right)
    {
        return left.second > right.second; // Sort in descending order of frequency
    });

    // Print vector sorting frequencies.
    std::cout << "\n\n\nWord frequencies:\n";

    for (auto const & item : ranking)
    {
        std::cout << item.first
                  << " -> "
                  << item.second
                  << '\n';
    }

    return 0;
}


// map is used for search, insert and update -> O(log n)
// vector is used for ranking  -> O(n)