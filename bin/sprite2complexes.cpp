#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <zlib.h>

std::string get_current_time() {
    // Get the current time
    auto now = std::chrono::system_clock::now();

    // Convert to time_t to get a time that we can format
    std::time_t currentTime = std::chrono::system_clock::to_time_t(now);

    // Convert to local time
    std::tm* localTime = std::localtime(&currentTime);

    // Create a string stream to hold the formatted time
    std::ostringstream oss;
    oss << std::put_time(localTime, "%Y-%m-%d %H:%M:%S");

    // Return the formatted string
    return oss.str();
}

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

void processClusterLine(const std::string& line, const std::unordered_map<std::string, int>& chromSizes,
                        const std::string& libid, int extbp, int selfbp, std::ofstream& fout, int& i,
                        int& num_frag, int& num_filtered, std::unordered_map<int, int>& histogram, int& max_frag, int& min_frag) {
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

    int num_frag_in_complex = 0;
    // Process valid chromosomal positions
    for (auto& [chrom, positions] : chromPositions) {
        num_frag += positions.size();
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
            num_frag_in_complex += 1;
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
                else {
                    num_filtered++;
                }
            }

            if (validPositions.size() == 1) {
                validPositions.clear(); // Discard the smallest position both if no other valid positions
            }

            int len = validPositions.size();
            num_frag_in_complex += len;
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

    if (num_frag_in_complex >= 6 && num_frag_in_complex <= 10) {
        histogram[6]++;
    } else if (num_frag_in_complex >= 11 && num_frag_in_complex <= 50) {
        histogram[11]++;
    } else if (num_frag_in_complex > 50) {
        histogram[51]++;
    } else {
        histogram[num_frag_in_complex]++;
    }

    if (num_frag_in_complex > max_frag) {
        max_frag = num_frag_in_complex;
    }

    if (num_frag_in_complex < min_frag && num_frag_in_complex != 0) {
        min_frag = num_frag_in_complex;
    }
}

void readSpriteAndWriteRegions(const std::string& directory, const std::string& spriteFile,
                                const std::unordered_map<std::string, int>& chromSizes, const std::string& libid,
                                int extbp, int selfbp, const std::string& logFile) {
    std::string outputFile = directory + "/" + libid + ".complexes";
    std::ofstream fout(outputFile);
    if (!fout.is_open()) {
        throw std::runtime_error("Unable to open output file at " + outputFile);
    }

    std::ofstream logfout(logFile);
    if (!logfout.is_open()) {
        throw std::runtime_error("Unable to open log file at " + logFile);
    }

    // Log the version number
    std::string version = "MIA-Sort sprite2complexes Version 0.1.1\n-------------------------------";
    logfout << version << std::endl;

    gzFile gz = gzopen(spriteFile.c_str(), "rb");
    if (!gz) {
        throw std::runtime_error("Unable to open sprite file at " + spriteFile);
    }

    std::string line;
    int i = 100000000;
    char buffer[8192];

    std::unordered_map<int, int> histogram;

    int num_frag = 0;
    int num_filtered = 0;
    int max_frag = INT_MIN;
    int min_frag = INT_MAX;
    int num_lines = 0;

    logfout << get_current_time() << " sprite2complexes starts\n" << std::endl;
    // Loop through the file, reading chunks until the end of the file
    while (gzgets(gz, buffer, sizeof(buffer)) != Z_NULL) {
        std::string part(buffer);
        while (part.back() != '\n' && gzgets(gz, buffer, sizeof(buffer)) != Z_NULL) {
            part.append(buffer);
        }
        line = part;
        num_lines++;
        processClusterLine(line, chromSizes, libid, extbp, selfbp, fout, i, num_frag, num_filtered, histogram, max_frag, min_frag);
    }

    int num_complexes = 0;
    // Iterate through the map and sum all the values
    for (const auto& pair : histogram) {
        num_complexes += pair.second;
    }

    logfout << "The number of fragments in the clusters file: " << num_frag << std::endl;
    logfout << "The total number of processed lines in the clusters file: " << num_lines << std::endl;
    logfout << "The number of lines filtered out due to selfbp: " << num_filtered << std::endl;
    logfout << "Total number of complexes written in the output complexes file: " << num_complexes << std::endl;

    logfout << "Histogram of fragments/complex in the output complexes file: " << std::endl;
    logfout << "1 frag: " << histogram[1] << std::endl;
    logfout << "2 frags: " << histogram[2] << std::endl;
    logfout << "3 frags: " << histogram[3] << std::endl;
    logfout << "4 frags: " << histogram[4] << std::endl;
    logfout << "5 frags: " << histogram[5] << std::endl;
    logfout << "6-10 frags: " << histogram[6] << std::endl;
    logfout << "11-50 frags: " << histogram[11] << std::endl;
    logfout << ">50 frags: " << histogram[51] << std::endl;
    logfout << "min frag: " << min_frag << std::endl;
    logfout << "max frag " << max_frag << std::endl;

    logfout << "\n" << get_current_time() << " sprite2complexes ends" << std::endl;

    gzclose(gz);
    fout.close();
    logfout.close();
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

    std::string logFile = directory + "/" + libid + ".log";

    try {
        std::unordered_map<std::string, int> chromSizes = readChromSizes(chromSizesFile);
        readSpriteAndWriteRegions(directory, spriteFile, chromSizes, libid, extbp, selfbp, logFile);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
