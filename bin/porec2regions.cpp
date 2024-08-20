#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <zlib.h>

void processLine(const std::string& line, std::vector<std::string>& result, std::string& currentReadName, int& fragmentCount) {
    std::istringstream ss(line);
    std::string field;
    std::vector<std::string> fields;
    while (std::getline(ss, field, ',')) {
        fields.push_back(field);
    }

    if (fields[17] == "True") {  // Column 18: pass_filter == "True"
        std::string chrom = fields[3];  // Column 4: chrom
        std::string start = fields[4];  // Column 5: start
        std::string end = fields[5];  // Column 6: end
        std::string readName = fields[7];  // Column 8: read_name

        if (readName != currentReadName) {
            if (!currentReadName.empty()) {
                // Write the previous group of fragments to result
                std::ostringstream oss;
                oss << chrom << "\t" << start << "\t" << end << "\t" << fragmentCount << "\t" << currentReadName;
                result.push_back(oss.str());
            }
            // Reset for the new group
            currentReadName = readName;
            fragmentCount = 1;
        } else {
            fragmentCount++;
        }
    }
}

void readCSVAndWriteRegions(const std::string& directory, const std::string& csvFile) {
    std::string outputFile = directory + "/GM12878_Pore-C_GSM4490689.hg38.region";  // TODO: revise hand-crafted filename
    std::ofstream fout(outputFile);
    if (!fout.is_open()) {
        throw std::runtime_error("Unable to open output file at " + outputFile);
    }

    gzFile gz = gzopen(csvFile.c_str(), "rb");
    if (!gz) {
        throw std::runtime_error("Unable to open gzipped CSV file at " + csvFile);
    }

    std::string line;
    std::vector<std::string> result;
    std::string currentReadName;
    int fragmentCount = 0;

    char buffer[8192];
    while (gzgets(gz, buffer, sizeof(buffer)) != Z_NULL) {
        line = buffer;

        // Continue reading the line until the newline character is found
        while (line.back() != '\n' && gzgets(gz, buffer, sizeof(buffer)) != Z_NULL) {
            line.append(buffer);
        }

        processLine(line, result, currentReadName, fragmentCount);
    }

    // Add the last group
    if (!currentReadName.empty()) {
        std::ostringstream oss;
        oss << result.back(); // last line
        fout << oss.str() << "\n";
    }

    for (const auto& res : result) {
        fout << res << "\n";
    }

    gzclose(gz);
    fout.close();
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <directory> <csv_file>" << std::endl;
        return 1;
    }

    std::string directory = argv[1];
    std::string csvFile = argv[2];

    try {
        readCSVAndWriteRegions(directory, csvFile);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
