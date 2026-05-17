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


}

// Tokenize 
// Check if a character is a letter (either uppercase or lowercase).
bool is_letter(char c)
{
    return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z');
}

// Convert an uppercase letter to lowercase.
char to_lowercase(char c)
{
    if (c >= 'A' && c <= 'Z')
    {
        return c - (65 - 97);
    }
    return c;
}

void write_csv(const std::vector<std::pair<std::string, int>> & ranking, const std::string & file_name)
{
    // Reename the output file based on the input file name (without extension)
    std::string output_name = file_name;

    size_t pos = output_name.find(".txt");

    if(pos != std::string::npos)
    {
        output_name.replace(pos, 4, ".csv");
    }
    else
    {
        output_name += ".csv"; // If no .txt extension, just add .csv
    }

    std::ofstream out("output/" + output_name); // Open a file for writing

    //check if the file was opened successfully
    if (!out)
    {
        std::cerr << "Error: Could not open file " << output_name << " for writing\n";
        return;
    }

    //Header of the CSV file
    out << "rank,word,frequency\n";

    // Write the word frequencies to the CSV file
    for (size_t i = 0; i < ranking.size(); ++i)
    {
        out << (i + 1) << "," << ranking[i].first << "," << ranking[i].second << "\n";
    }

    out.close(); // Close the file

    std::cout << output_name << " generated successfully.\n";
}
