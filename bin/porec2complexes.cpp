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
    auto now = std::chrono::system_clock::now();
    std::time_t currentTime = std::chrono::system_clock::to_time_t(now);
    std::tm* localTime = std::localtime(&currentTime);

    std::ostringstream oss;
    oss << std::put_time(localTime, "%Y-%m-%d %H:%M:%S");

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

                for (const auto& fragment : fragments) {
                    fout << fragment.chrom << "\t" << fragment.start << "\t" << fragment.end << "\t" << fragmentCount << "\t" << fragment.readName << "\n";
                }
                fragments.clear();  // Clear fragments for the next readName
            }
            currentReadName = readName;  // Update currentReadName to the new readName
            fragmentCount = 0;  // Reset fragment count for the new group
        }

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

void readCSVAndWriteRegions(const std::string& csvFile, const std::string& outputFile, const std::string& logFile, int argc, char* argv[]) {
    std::ofstream fout(outputFile);
    if (!fout.is_open()) {
        throw std::runtime_error("Unable to open output file at " + outputFile);
    }

    std::ofstream logfout(logFile);
    if (!logfout.is_open()) {
        throw std::runtime_error("Unable to open log file at " + logFile);
    }

    std::string version = "MIA-Sort porec2complexes Version 0.1.2\n-------------------------------";
    logfout << version << "\n";

    std::string command;
    for (int i = 0; i < argc; ++i) {
        command += argv[i];
        if (i < argc - 1) {
            command += " ";
        }
    }
    logfout << "User Command: " << command << "\n\n";

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

        while (!line.empty() && line.back() != '\n' && gzgets(gz, buffer, sizeof(buffer)) != Z_NULL) {
            line.append(buffer);
        }

        if (!line.empty()) {
            num_lines++;
            processLine(line, fragments, currentReadName, fragmentCount, fout, num_frag, num_filtered, histogram, max_frag, min_frag);
        }
    }

    if (!fragments.empty()) {
        writeFragments(fout, fragments, fragmentCount, histogram, min_frag, max_frag);
    }

    long long int num_complexes = 0;
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
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <directory> <csv_file> <output_file>" << std::endl;
        return 1;
    }

    std::string directory = argv[1];
    std::string csvFile = argv[2];
    std::string outputFile = argv[3];

    std::string inputFileName = csvFile.substr(csvFile.find_last_of("/") + 1);
    inputFileName = inputFileName.substr(0, inputFileName.find(".csv.gz"));
    inputFileName = inputFileName.substr(0, inputFileName.find(".csv"));
    std::string logFile = directory + "/porec2complexes_" + inputFileName + ".log";

    try {
        readCSVAndWriteRegions(csvFile, outputFile, logFile, argc, argv);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
