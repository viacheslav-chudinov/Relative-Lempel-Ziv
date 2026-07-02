#include <fstream>
#include <stdexcept>
#include <string>

std::string read_file(const std::string& file_name){
	// Reading a single string from a file into std::string
    std::ifstream file(file_name);

    if (!file.is_open()){
        throw std::runtime_error("Cannot open file: " + file_name);
    }

    std::string string;
    std::getline(file, string);
    return string;
}