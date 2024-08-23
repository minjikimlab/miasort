#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <zlib.h>
#include <stdexcept>
#include <algorithm>
#include <vector>
#include <chrono>
#include <iomanip>
#include <ctime>

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

void processLine(const std::string& line, const std::unordered_map<std::string, int>& chromSizes, const std::string& libid,
                int extbp, int selfbp, std::ofstream& fout, int& i, int& num_filtered, int& num_complexes, int& num_lines) {
    if (line.empty() || line[0] == '#') {
        return;
    }

    std::istringstream ss(line);
    std::vector<std::string> fields;
    std::string field;
    while (std::getline(ss, field, '\t')) {
        fields.push_back(field);
    }

    num_lines++;

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
        num_complexes++;
    }
    else {
        num_filtered++;
    }
}

void readPairsAndWriteRegions(const std::string& directory, const std::string& pairsFile,
                            const std::unordered_map<std::string, int>& chromSizes, const std::string& libid,
                            int extbp, int selfbp, const std::string& logFile, int argc, char* argv[]) {
    std::string outputFile = directory + "/" + libid + ".complexes";

    int num_filtered = 0;
    int num_complexes = 0;
    int num_lines = 0;

    std::ofstream fout(outputFile);
    if (!fout.is_open()) {
        throw std::runtime_error("Unable to open output file at " + outputFile);
    }

    std::ofstream logfout(logFile);
    if (!logfout.is_open()) {
        throw std::runtime_error("Unable to open log file at " + logFile);
    }

    // Log the version number
    std::string version = "MIA-Sort pairs2complexes Version 0.1.1\n-----------------------------------";
    logfout << version << std::endl;

    std::string command;
    // Loop through all arguments and concatenate them into a single string
    for (int i = 0; i < argc; ++i) {
        command += argv[i];
        if (i < argc - 1) {
            command += " "; // Add a space between arguments
        }
    }

    // Output the exact command the user ran
    logfout << "User Command: " << command << "\n\n";

    gzFile gz = gzopen(pairsFile.c_str(), "rb");
    if (!gz) {
        throw std::runtime_error("Unable to open pairs file at " + pairsFile);
    }

    char buffer[8192];
    std::string line;
    int i = 100000000;

    logfout << get_current_time() << " pairs2complexes starts\n" << std::endl;
    while (gzgets(gz, buffer, sizeof(buffer)) != Z_NULL) {
        line = buffer;
        // Remove newline character if it exists
        if (!line.empty() && line.back() == '\n') {
            line.pop_back();
        }
        processLine(line, chromSizes, libid, extbp, selfbp, fout,
                    i, num_filtered, num_complexes, num_lines);
    }

    logfout << "The total number of processed lines in the pairs file: " << num_lines << std::endl;
    logfout << "The number of pairs lines filtered out due to selfbp: " << num_filtered << std::endl;
    logfout << "The number of complexes written in the output file: " << num_complexes << std::endl;

    logfout << "\n" << get_current_time() << " pairs2complexes ends" << std::endl;

    gzclose(gz);
    fout.close();
    logfout.close();
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

    std::string logFile = directory + "/" + libid + ".log";

    try {
        std::unordered_map<std::string, int> chromSizes = readChromSizes(chromSizesFile);
        readPairsAndWriteRegions(directory, pairsFile, chromSizes, libid, extbp, selfbp, logFile, argc, argv);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
