#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <zlib.h>

struct Fragment {
    std::string chrom;
    std::string start;
    std::string end;
    std::string readName;
};

void processLine(const std::string& line, std::vector<Fragment>& fragments, std::string& currentReadName, int& fragmentCount, std::ofstream& fout) {
    std::istringstream ss(line);
    std::string field;
    std::vector<std::string> fields;

    while (std::getline(ss, field, ',')) {
        fields.push_back(field);
    }

    if (fields.size() > 17 && fields[17] == "True") {  // Column 18: pass_filter == "True"
        std::string chrom = fields[3];  // Column 4: chrom
        std::string start = fields[4];  // Column 5: start
        std::string end = fields[5];  // Column 6: end
        std::string readName = fields[7];  // Column 8: read_name

        if (readName != currentReadName) {
            if (!currentReadName.empty()) {
                // Output all stored fragments for the previous readName
                for (const auto& fragment : fragments) {
                    fout << fragment.chrom << "\t" << fragment.start << "\t" << fragment.end << "\t" << fragmentCount << "\t" << fragment.readName << "\n";
                }
                fragments.clear();  // Clear fragments for the next readName
            }
            currentReadName = readName;  // Update currentReadName to the new readName
            fragmentCount = 0;  // Reset fragment count for the new group
        }

        // Store the current fragment
        fragments.push_back({chrom, start, end, readName});
        fragmentCount++;
    }
}

void writeFragments(std::ofstream& fout, const std::vector<Fragment>& fragments, int fragmentCount) {
    for (const auto& fragment : fragments) {
        fout << fragment.chrom << "\t" << fragment.start << "\t" << fragment.end << "\t" << fragmentCount << "\t" << fragment.readName << "\n";
    }
}

void readCSVAndWriteRegions(const std::string& directory, const std::string& csvFile) {
    std::string outputFile = directory + "/GM12878_Pore-C_GSM4490689.hg38.region";
    std::ofstream fout(outputFile);
    if (!fout.is_open()) {
        throw std::runtime_error("Unable to open output file at " + outputFile);
    }

    gzFile gz = gzopen(csvFile.c_str(), "rb");
    if (!gz) {
        throw std::runtime_error("Unable to open gzipped CSV file at " + csvFile);
    }

    std::string line;
    std::string currentReadName;
    int fragmentCount = 0;
    std::vector<Fragment> fragments;

    char buffer[8192];
    while (gzgets(gz, buffer, sizeof(buffer)) != Z_NULL) {
        line = buffer;

        // Continue reading the line until the newline character is found
        while (!line.empty() && line.back() != '\n' && gzgets(gz, buffer, sizeof(buffer)) != Z_NULL) {
            line.append(buffer);
        }

        if (!line.empty()) {
            processLine(line, fragments, currentReadName, fragmentCount, fout);
        }
    }

    // Write the last group of fragments to the output file
    if (!fragments.empty()) {
        writeFragments(fout, fragments, fragmentCount);
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
