#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <algorithm>
#include <vector>
#include <zlib.h>

// Function to read chromosome sizes from a file
std::map<std::string, int> readChromSizes(const std::string& filepath) {
    std::map<std::string, int> chromSizes;
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Unable to open chromosome sizes file at " + filepath);
    }

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream ss(line);
        std::string chrom;
        int size;
        ss >> chrom >> size;
        chromSizes[chrom] = size;
    }
    file.close();
    return chromSizes;
}


// Function to read pairs and write regions
void readPairsAndWriteRegions(const std::string& directory, const std::string& pairsFile, const std::map<std::string, int>& chromSizes, const std::string& libid, int extbp, int selfbp) {
    std::string outputFile = directory + "/" + libid + ".ext" + std::to_string(extbp) + "bp.g" + std::to_string(selfbp) + "bp.region";
    std::ofstream fout(outputFile, std::ios::app);
    if (!fout.is_open()) {
        throw std::runtime_error("Unable to open output file at " + outputFile);
    }

    gzFile gz = gzopen(pairsFile.c_str(), "rb");
    if (!gz) {
        throw std::runtime_error("Unable to open pairs file at " + pairsFile);
    }

    char buffer[8192];
    std::string line;
    int i = 100000000;
    while (gzgets(gz, buffer, sizeof(buffer)) != Z_NULL) {
        line = buffer;
        if (line[0] != '#') {
            std::istringstream ss(line);
            std::vector<std::string> fields;
            std::string field;
            while (std::getline(ss, field, '\t')) {
                fields.push_back(field);
            }

            std::string chrom1 = fields[1];
            int pos1 = std::stoi(fields[2]);
            std::string chrom2 = fields[3];
            int pos2 = std::stoi(fields[4]);

            if (chrom1 == chrom2 && chrom1 != "chrM" && pos2 - pos1 > selfbp) {
                int start1 = std::max(0, pos1 - extbp);
                int end1 = std::min(pos1 + extbp, chromSizes.at(chrom1));
                int start2 = std::max(0, pos2 - extbp);
                int end2 = std::min(pos2 + extbp, chromSizes.at(chrom1));
                std::string gemid = libid + "-100-" + std::to_string(i) + "-HEA-7-4-sub-1-1";

                fout << chrom1 << "\t" << start1 << "\t" << end1 << "\t2\t" << gemid << "\n";
                fout << chrom1 << "\t" << start2 << "\t" << end2 << "\t2\t" << gemid << "\n";
                ++i;
            }
        }
    }
    gzclose(gz);
    fout.close();
}

int main(int argc, char* argv[]) {
    if (argc < 4 || argc > 6) {
        std::cerr << "Usage: " << argv[0] << " <directory> <pairs_file> <chrom_sizes_file> [<extbp> [<selfbp>]]" << std::endl;
        return 1;
    }

    std::string directory = argv[1];
    std::string pairsFile = argv[2];
    std::string chromSizesFile = argv[3];

    // User can specify these values
    int extbp = (argc > 4) ? std::stoi(argv[4]) : 250;
    int selfbp = (argc > 5) ? std::stoi(argv[5]) : 8000;

    std::string libid = pairsFile.substr(pairsFile.find_last_of("/") + 1, pairsFile.find(".bsorted.pairs.gz") - pairsFile.find_last_of("/") - 1);

    try {
        std::map<std::string, int> chromSizes = readChromSizes(chromSizesFile);
        readPairsAndWriteRegions(directory, pairsFile, chromSizes, libid, extbp, selfbp);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
