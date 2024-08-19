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

class ThreadPool {
public:
    ThreadPool(size_t threads);
    ~ThreadPool();

    template<class F, class... Args>
    auto enqueue(F&& f, Args&&... args) -> std::future<typename std::result_of<F(Args...)>::type>;

private:
    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;
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
        if (stop)
            throw std::runtime_error("enqueue on stopped ThreadPool");
        tasks.emplace([task]() { (*task)(); });
    }
    condition.notify_one();
    return res;
}

std::vector<std::string> processPoreCDataBatch(const std::vector<std::string>& lines, int& read_idx) {
    std::vector<std::string> result;
    for (const auto& line : lines) {
        if (!line.empty() && line[0] != '#') {
            std::istringstream ss(line);
            std::vector<std::string> fields;
            std::string field;

            while (std::getline(ss, field, ',')) {  // Adjusted to parse CSV
                fields.push_back(field);
            }

            // Correct the column indices based on your file structure
            std::string pass_filter = fields[17];             // 'pass_filter' column index
            std::string chrom = fields[3];                    // 'chrom' column index
            std::string start = fields[4];                    // 'start' column index
            std::string end = fields[5];                      // 'end' column index
            std::string readidx = fields[0];                  // 'readidx' column index
            std::string num_fragments = fields[20];           // 'num_overlapping_fragments' column index

            if (pass_filter == "True") {
                std::string output = chrom + "\t" + start + "\t" + end + "\t" + num_fragments + "\t" + readidx + "-" + std::to_string(read_idx);
                result.push_back(output);
                ++read_idx;
            }
        }
    }
    return result;
}


void readPoreCAndWriteRegions(const std::string& directory, const std::string& poreCFile) {
    std::string outputFile = directory + "/PoreC_regions.region";
    std::ofstream fout(outputFile);
    if (!fout.is_open()) {
        throw std::runtime_error("Unable to open output file at " + outputFile);
    }

    gzFile gz = gzopen(poreCFile.c_str(), "rb");
    if (!gz) {
        throw std::runtime_error("Unable to open Pore-C file at " + poreCFile);
    }

    char buffer[8192];
    std::string line;
    int read_idx = 1;
    std::mutex mtx;
    ThreadPool pool(std::thread::hardware_concurrency());
    std::vector<std::future<std::vector<std::string>>> futures;
    std::vector<std::string> batch;
    const size_t BATCH_SIZE = 100000;  // Adjust batch size as needed

    while (gzgets(gz, buffer, sizeof(buffer)) != Z_NULL) {
        line = buffer;
        batch.push_back(line);
        if (batch.size() >= BATCH_SIZE) {
            futures.emplace_back(pool.enqueue([batch, &read_idx]() mutable {
                return processPoreCDataBatch(batch, read_idx);
            }));
            batch.clear();
        }
    }

    if (!batch.empty()) {
        futures.emplace_back(pool.enqueue([batch, &read_idx]() mutable {
            return processPoreCDataBatch(batch, read_idx);
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
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <directory> <porec_file>" << std::endl;
        return 1;
    }

    std::string directory = argv[1];
    std::string poreCFile = argv[2];

    try {
        readPoreCAndWriteRegions(directory, poreCFile);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
