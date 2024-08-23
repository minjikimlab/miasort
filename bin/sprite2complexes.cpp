#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <zlib.h>

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

void processClusterLine(const std::string& line, const std::unordered_map<std::string, int>& chromSizes, const std::string& libid, int extbp, int selfbp, std::ofstream& fout, int& i) {
    std::istringstream ss(line);
    std::string id;
    ss >> id;

    std::unordered_map<std::string, std::vector<int>> chromPositions;
    std::string location;
    while (ss >> location) {
        size_t colonPos = location.find(':');
        if (colonPos == std::string::npos || colonPos == location.length() - 1) {
            std::cerr << "Error: Malformed location string: " << location << std::endl;
            continue; // Skip malformed locations
        }

        std::string chrom = location.substr(0, colonPos);
        std::string posStr = location.substr(colonPos + 1);

        // Skip malformed position strings
        if (posStr.empty() || !std::all_of(posStr.begin(), posStr.end(), ::isdigit)) {
            std::cerr << "Error: Invalid position string: " << posStr << " in line: " << line << std::endl;
            continue;
        }

        try {
            int pos = std::stoi(posStr);
            chromPositions[chrom].push_back(pos);
        } catch (const std::invalid_argument& e) {
            std::cerr << "Error: Invalid position string: " << posStr << " in line: " << line << std::endl;
            continue; // Skip this location
        } catch (const std::out_of_range& e) {
            std::cerr << "Error: Position out of range: " << posStr << " in line: " << line << std::endl;
            continue; // Skip this location
        }
    }

    // Process valid chromosomal positions
    for (auto& [chrom, positions] : chromPositions) {
        if (positions.size() == 1) {
            int pos1 = positions[0];
            auto it = chromSizes.find(chrom);
            if (it == chromSizes.end()) {
                continue;
            }
            int chromSize = it->second;
            int start1 = std::max(0, pos1 - extbp);
            int end1 = std::min(pos1 + extbp, chromSize);
            ++i;
            fout << chrom << "\t" << start1 << "\t" << end1 << "\t1\t" << id << "\n";
        } else if (positions.size() > 1) {
            std::sort(positions.begin(), positions.end());

            std::vector<int> validPositions;
            validPositions.push_back(positions.front());  // Always keep the smallest position

            for (size_t j = 1; j < positions.size(); ++j) {
                int pos1 = positions[0];
                int pos2 = positions[j];
                if (pos2 - pos1 > selfbp) {
                    validPositions.push_back(positions[j]);
                }
            }

            if (validPositions.size() == 1) {
                validPositions.clear(); // Discard the smallest position both if no other valid positions
            }

            int len = validPositions.size();
            for (int pos : validPositions) {
                auto it = chromSizes.find(chrom);
                if (it == chromSizes.end()) {
                    continue;
                }
                int chromSize = it->second;
                int start1 = std::max(0, pos - extbp);
                int end1 = std::min(pos + extbp, chromSize);
                ++i;
                fout << chrom << "\t" << start1 << "\t" << end1 << "\t" << len << "\t" << id << "\n";
            }
        }
    }
}

void readSpriteAndWriteRegions(const std::string& directory, const std::string& spriteFile, const std::unordered_map<std::string, int>& chromSizes, const std::string& libid, int extbp, int selfbp) {
    std::string outputFile = directory + "/" + libid + ".ext" + std::to_string(extbp) + "bp.g" + std::to_string(selfbp) + "bp.complexes";
    std::ofstream fout(outputFile);
    if (!fout.is_open()) {
        throw std::runtime_error("Unable to open output file at " + outputFile);
    }

    gzFile gz = gzopen(spriteFile.c_str(), "rb");
    if (!gz) {
        throw std::runtime_error("Unable to open sprite file at " + spriteFile);
    }

    std::string line;
    int i = 100000000;
    char buffer[8192];

    // Loop through the file, reading chunks until the end of the file
    while (gzgets(gz, buffer, sizeof(buffer)) != Z_NULL) {
        std::string part(buffer);
        while (part.back() != '\n' && gzgets(gz, buffer, sizeof(buffer)) != Z_NULL) {
            part.append(buffer);
        }
        line = part;

        processClusterLine(line, chromSizes, libid, extbp, selfbp, fout, i);
    }

    gzclose(gz);
    fout.close();
}

int main(int argc, char* argv[]) {
    if (argc < 4 || argc > 6) {
        std::cerr << "Usage: " << argv[0] << " <directory> <sprite_file> <chrom_sizes_file> [<extbp> [<selfbp>]]" << std::endl;
        return 1;
    }

    std::string directory = argv[1];
    std::string spriteFile = argv[2];
    std::string chromSizesFile = argv[3];

    int extbp = (argc > 4) ? std::stoi(argv[4]) : 250;
    int selfbp = (argc > 5) ? std::stoi(argv[5]) : 8000;

    std::string libid = spriteFile.substr(spriteFile.find_last_of("/") + 1, spriteFile.find(".clusters") - spriteFile.find_last_of("/") - 1);

    try {
        std::unordered_map<std::string, int> chromSizes = readChromSizes(chromSizesFile);
        readSpriteAndWriteRegions(directory, spriteFile, chromSizes, libid, extbp, selfbp);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
