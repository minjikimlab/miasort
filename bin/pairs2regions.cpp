#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <zlib.h>
#include <stdexcept>
#include <algorithm>
#include <vector>

std::unordered_map<std::string, int> readChromSizes(const std::string& filepath) {
    std::unordered_map<std::string, int> chromSizes;
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

void processLine(const std::string& line, const std::unordered_map<std::string, int>& chromSizes, const std::string& libid, int extbp, int selfbp, std::ofstream& fout, int& i) {
    if (line.empty() || line[0] == '#') {
        return;
    }

    std::istringstream ss(line);
    std::vector<std::string> fields;
    std::string field;
    while (std::getline(ss, field, '\t')) {
        fields.push_back(field);
    }

    if (fields.size() < 5) {
        std::cerr << "Skipping malformed line: " << line << std::endl;
        return;
    }

    std::string gemid = fields[0];
    std::string chrom1 = fields[1];
    std::string chrom2 = fields[3];

    int pos1, pos2;
    try {
        pos1 = std::stoi(fields[2]);
        pos2 = std::stoi(fields[4]);
    } catch (const std::invalid_argument& e) {
        std::cerr << "Invalid position value in line: " << line << std::endl;
        return;
    }

    if (chrom1 == chrom2 && chrom1 != "chrM" && pos2 - pos1 > selfbp) {
        auto it = chromSizes.find(chrom1);
        if (it == chromSizes.end()) {
            return; // Skip if chromosome not found
        }
        int chromSize = it->second;

        int start1 = std::max(0, pos1 - extbp);
        int end1 = std::min(pos1 + extbp, chromSize);
        int start2 = std::max(0, pos2 - extbp);
        int end2 = std::min(pos2 + extbp, chromSize);
        ++i;

        fout << chrom1 << "\t" << start1 << "\t" << end1 << "\t2\t" << gemid << "\n";
        fout << chrom1 << "\t" << start2 << "\t" << end2 << "\t2\t" << gemid << "\n";
    }
}

void readPairsAndWriteRegions(const std::string& directory, const std::string& pairsFile, const std::unordered_map<std::string, int>& chromSizes, const std::string& libid, int extbp, int selfbp) {
    std::string outputFile = directory + "/" + libid + ".ext" + std::to_string(extbp) + "bp.g" + std::to_string(selfbp) + "bp.region";
    std::ofstream fout(outputFile);
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
        // Remove newline character if it exists
        if (!line.empty() && line.back() == '\n') {
            line.pop_back();
        }
        processLine(line, chromSizes, libid, extbp, selfbp, fout, i);
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

    int extbp = (argc > 4) ? std::stoi(argv[4]) : 250;
    int selfbp = (argc > 5) ? std::stoi(argv[5]) : 8000;

    std::string libid = pairsFile.substr(pairsFile.find_last_of("/") + 1, pairsFile.find(".bsorted.pairs.gz") - pairsFile.find_last_of("/") - 1);

    try {
        std::unordered_map<std::string, int> chromSizes = readChromSizes(chromSizesFile);
        readPairsAndWriteRegions(directory, pairsFile, chromSizes, libid, extbp, selfbp);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
