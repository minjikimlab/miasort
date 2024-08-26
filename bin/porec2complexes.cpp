#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <zlib.h>
#include <chrono>
#include <ctime>
#include <unordered_map>
#include <iomanip>

struct Fragment {
    std::string chrom;
    std::string start;
    std::string end;
    std::string readName;
};

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

void processLine(const std::string& line, std::vector<Fragment>& fragments,
                std::string& currentReadName, long long int& fragmentCount, std::ofstream& fout,
                long long int& num_frag, long long int& num_filtered, std::unordered_map<int, long long int>& histogram,
                long long int& max_frag, long long int& min_frag) {
    std::istringstream ss(line);
    std::string field;
    std::vector<std::string> fields;

    while (std::getline(ss, field, ',')) {
        fields.push_back(field);
    }

    num_frag++;

    if (fields.size() > 17 && fields[17] == "True") {  // Column 18: pass_filter == "True"
        std::string chrom = fields[3];  // Column 4: chrom
        std::string start = fields[4];  // Column 5: start
        std::string end = fields[5];  // Column 6: end
        std::string readName = fields[7];  // Column 8: read_name

        if (readName != currentReadName) {
            if (!currentReadName.empty()) {
                int num_frag_in_complex = fragments.size();
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
    } else {
        num_filtered++;
    }
}

void writeFragments(std::ofstream& fout, const std::vector<Fragment>& fragments,
                    int fragmentCount, std::unordered_map<int, long long int>& histogram,
                    long long int& min_frag, long long int& max_frag) {
    int num_frag_in_complex = fragments.size();
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
    for (const auto& fragment : fragments) {
        fout << fragment.chrom << "\t" << fragment.start << "\t" << fragment.end << "\t" << fragmentCount << "\t" << fragment.readName << "\n";
    }
}

void readCSVAndWriteRegions(const std::string& directory, const std::string& csvFile, const std::string& logFile) {
    // Extract the filename without extension
    std::string baseName = csvFile.substr(csvFile.find_last_of("/") + 1);
    baseName = baseName.substr(0, baseName.find("."));

    std::string outputFile = directory + "/" + baseName + ".complexes";
    std::ofstream fout(outputFile);
    if (!fout.is_open()) {
        throw std::runtime_error("Unable to open output file at " + outputFile);
    }

    std::ofstream logfout(logFile);
    if (!logfout.is_open()) {
        throw std::runtime_error("Unable to open log file at " + logFile);
    }

    // Log the version number and timestamp
    std::string version = "MIA-Sort porec2complexes Version 0.1.1\n-------------------------------";
    logfout << version << "\n";

    gzFile gz = gzopen(csvFile.c_str(), "rb");
    if (!gz) {
        throw std::runtime_error("Unable to open gzipped CSV file at " + csvFile);
    }

    std::string line;
    std::string currentReadName;
    long long int fragmentCount = 0;
    std::vector<Fragment> fragments;

    std::unordered_map<int, long long int> histogram;

    long long int num_frag = 0;
    long long int num_filtered = 0;
    long long int max_frag = LLONG_MIN;
    long long int min_frag = LLONG_MAX;
    long long int num_lines = 0;

    logfout << get_current_time() << " porec2complexes starts\n" << std::endl;

    char buffer[8192];
    while (gzgets(gz, buffer, sizeof(buffer)) != Z_NULL) {
        line = buffer;

        // Continue reading the line until the newline character is found
        while (!line.empty() && line.back() != '\n' && gzgets(gz, buffer, sizeof(buffer)) != Z_NULL) {
            line.append(buffer);
        }

        if (!line.empty()) {
            num_lines++;
            processLine(line, fragments, currentReadName, fragmentCount, fout, num_frag, num_filtered, histogram, max_frag, min_frag);
        }
    }

    // Write the last group of fragments to the output file
    if (!fragments.empty()) {
        writeFragments(fout, fragments, fragmentCount, histogram, min_frag, max_frag);
    }

    long long int num_complexes = 0;
    // Iterate through the map and sum all the values
    for (const auto& pair : histogram) {
        num_complexes += pair.second;
    }

    logfout << "The number of fragments in the pore-c file: " << num_frag << std::endl;
    logfout << "The total number of processed lines in the pore-c file: " << num_lines << std::endl;
    logfout << "The number of lines filtered out due to pass_filter: " << num_filtered << std::endl;
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

    logfout << "\n" << get_current_time() << " porec2complexes ends" << std::endl;

    gzclose(gz);
    fout.close();
    logfout.close();
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <directory> <csv_file>" << std::endl;
        return 1;
    }

    std::string directory = argv[1];
    std::string csvFile = argv[2];

    // Construct log file name based on the input file name
    std::string logFile = directory + "/" + csvFile.substr(csvFile.find_last_of("/") + 1, csvFile.find(".") - csvFile.find_last_of("/") - 1) + ".log";

    try {
        readCSVAndWriteRegions(directory, csvFile, logFile);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
