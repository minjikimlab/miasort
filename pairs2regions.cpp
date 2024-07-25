#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <zlib.h>
#include <thread>
#include <mutex>
#include <queue>
#include <future>
#include <stdexcept>
#include <condition_variable>

class ThreadPool {
public:
    ThreadPool(size_t threads);
    ~ThreadPool();

    template<class F, class... Args>
    auto enqueue(F&& f, Args&&... args) -> std::future<typename std::result_of<F(Args...)>::type>;

private:
    // Workers
    std::vector<std::thread> workers;
    // Task queue
    std::queue<std::function<void()>> tasks;

    // Synchronization
    std::mutex queue_mutex;
    std::condition_variable condition;
    bool stop;
};


inline ThreadPool::ThreadPool(size_t threads) : stop(false) {
    for (size_t i = 0; i < threads; ++i)
        workers.emplace_back([this] {
            for (;;) {
                std::function<void()> task;

                {
                    std::unique_lock<std::mutex> lock(this->queue_mutex);
                    this->condition.wait(lock, [this] { return this->stop || !this->tasks.empty(); });
                    if (this->stop && this->tasks.empty())
                        return;
                    task = std::move(this->tasks.front());
                    this->tasks.pop();
                }

                task();
            }
        });
}


inline ThreadPool::~ThreadPool() {
    {
        std::unique_lock<std::mutex> lock(queue_mutex);
        stop = true;
    }
    condition.notify_all();
    for (std::thread &worker : workers)
        worker.join();
}


template<class F, class... Args>
auto ThreadPool::enqueue(F&& f, Args&&... args) -> std::future<typename std::result_of<F(Args...)>::type> {
    using return_type = typename std::result_of<F(Args...)>::type;

    auto task = std::make_shared<std::packaged_task<return_type()>>(std::bind(std::forward<F>(f), std::forward<Args>(args)...));

    std::future<return_type> res = task->get_future();
    {
        std::unique_lock<std::mutex> lock(queue_mutex);

        // Don't allow enqueueing after stopping the pool
        if (stop)
            throw std::runtime_error("enqueue on stopped ThreadPool");

        tasks.emplace([task]() { (*task)(); });
    }
    condition.notify_one();
    return res;
}


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


// Function to process a line and return the results
std::vector<std::string> processLine(const std::string& line, const std::map<std::string, int>& chromSizes, const std::string& libid, int extbp, int selfbp, int& i) {
    std::vector<std::string> result;
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
            ++i;

            result.push_back(chrom1 + "\t" + std::to_string(start1) + "\t" + std::to_string(end1) + "\t2\t" + gemid);
            result.push_back(chrom1 + "\t" + std::to_string(start2) + "\t" + std::to_string(end2) + "\t2\t" + gemid);
        }
    }
    return result;
}


// Function to read pairs and write regions using a thread pool
void readPairsAndWriteRegions(const std::string& directory, const std::string& pairsFile,
                             const std::map<std::string, int>& chromSizes,
                             const std::string& libid, int extbp, int selfbp) {
    std::string outputFile = directory + "/" + libid + ".ext" + std::to_string(extbp) +
                            "bp.g" + std::to_string(selfbp) + "bp.region";
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
    std::mutex mtx;
    ThreadPool pool(std::thread::hardware_concurrency());

    std::vector<std::future<std::vector<std::string>>> futures;

    while (gzgets(gz, buffer, sizeof(buffer)) != Z_NULL) {
        line = buffer;
        // Enqueue a task to process the line
        futures.emplace_back(pool.enqueue([&chromSizes, &libid, extbp, selfbp, line, &i]() mutable {
            return processLine(line, chromSizes, libid, extbp, selfbp, i);
        }));
    }

    for (auto& fut : futures) {
        std::vector<std::string> results = fut.get();
        std::lock_guard<std::mutex> lock(mtx);
        for (const auto& res : results) {
            fout << res << "\n";
        }
    }

    gzclose(gz);
    fout.close();
}


int main(int argc, char* argv[]) {
    if (argc < 4 || argc > 6) {
        std::cerr << "Usage: " << argv[0] <<
        " <directory> <pairs_file> <chrom_sizes_file> [<extbp> [<selfbp>]]" << std::endl;
        return 1;
    }

    std::string directory = argv[1];
    std::string pairsFile = argv[2];
    std::string chromSizesFile = argv[3];

    // User can specify these values
    int extbp = (argc > 4) ? std::stoi(argv[4]) : 250;
    int selfbp = (argc > 5) ? std::stoi(argv[5]) : 8000;

    std::string libid = pairsFile.substr(pairsFile.find_last_of("/") + 1,
                                        pairsFile.find(".bsorted.pairs.gz") - pairsFile.find_last_of("/") - 1);

    try {
        std::map<std::string, int> chromSizes = readChromSizes(chromSizesFile);
        readPairsAndWriteRegions(directory, pairsFile, chromSizes, libid, extbp, selfbp);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
