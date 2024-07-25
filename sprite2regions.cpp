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


std::vector<std::string> processSpriteBatch(const std::vector<std::string>& lines, const std::map<std::string, int>& chromSizes, int extbp, int& cluster_id) {
    std::vector<std::string> result;
    for (const auto& line : lines) {
        if (!line.empty() && line[0] != '#') {
            std::istringstream ss(line);
            std::string barcodes, chrom;
            int pos;
            std::vector<std::string> fields;
            while (std::getline(ss, barcodes, '\t')) {
                fields.push_back(barcodes);
            }
            for (size_t j = 1; j < fields.size(); ++j) {
                std::istringstream loc(fields[j]);
                loc >> chrom >> pos;
                int start = std::max(0, pos - extbp);
                int end = std::min(pos + extbp, chromSizes.at(chrom));
                std::string region_id = "SPRITE-" + std::to_string(cluster_id) + "-R" + std::to_string(j);
                result.push_back(chrom + "\t" + std::to_string(start) + "\t" + std::to_string(end) + "\t" + region_id);
            }
            ++cluster_id;
        }
    }
    return result;
}


void readSpriteAndWriteRegions(const std::string& directory, const std::string& spriteFile, const std::map<std::string, int>& chromSizes, int extbp) {
    std::string outputFile = directory + "/SPRITE_regions.ext" + std::to_string(extbp) + "bp.region";
    std::ofstream fout(outputFile);
    if (!fout.is_open()) {
        throw std::runtime_error("Unable to open output file at " + outputFile);
    }

    gzFile gz = gzopen(spriteFile.c_str(), "rb");
    if (!gz) {
        throw std::runtime_error("Unable to open SPRITE file at " + spriteFile);
    }

    char buffer[8192];
    std::string line;
    int cluster_id = 1;
    std::mutex mtx;
    ThreadPool pool(std::thread::hardware_concurrency());
    std::vector<std::future<std::vector<std::string>>> futures;
    std::vector<std::string> batch;
    const size_t BATCH_SIZE = 100000;  // hyperparam, TODO:

    while (gzgets(gz, buffer, sizeof(buffer)) != Z_NULL) {
        line = buffer;
        batch.push_back(line);
        if (batch.size() >= BATCH_SIZE) {
            futures.emplace_back(pool.enqueue([&chromSizes, extbp, batch, &cluster_id]() mutable {
                return processSpriteBatch(batch, chromSizes, extbp, cluster_id);
            }));
            batch.clear();
        }
    }

    if (!batch.empty()) {
        futures.emplace_back(pool.enqueue([&chromSizes, extbp, batch, &cluster_id]() mutable {
            return processSpriteBatch(batch, chromSizes, extbp, cluster_id);
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
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <directory> <sprite_file> <chrom_sizes_file> [<extbp>]" << std::endl;
        return 1;
    }

    std::string directory = argv[1];
    std::string spriteFile = argv[2];
    std::string chromSizesFile = argv[3];

    int extbp = (argc > 4) ? std::stoi(argv[4]) : 250;

    try {
        std::map<std::string, int> chromSizes = readChromSizes(chromSizesFile);
        readSpriteAndWriteRegions(directory, spriteFile, chromSizes, extbp);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
