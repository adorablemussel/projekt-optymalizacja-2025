#include "dataMatrix.h"
#include <fstream>

matrix dataMatrix(int rows, int cols, std::string filename)
{
    matrix result;
    result = matrix(rows, cols);

    std::ifstream file(filename);
    if (!file.is_open()) {
        throw string("Nie mozna otworzyć pliku: " + filename);
    }
    std::string line;
    int i = 0;
    while (std::getline(file, line)){
        int j = 0;
        std::istringstream lineStream(line);
        std::string value;
        while (lineStream >> value) {
            result(i, j) = std::stod(value);
            ++j;
        }
        ++i;
    }

    return result;
}