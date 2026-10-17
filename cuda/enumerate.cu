// GPU off-target enumeration. Same hit definition as cpp/enumerate.cpp:
// 1-based 23 bp + strand window, mismatches on the 20 bp protospacer,
// PAM classes such as NGG / NAG / NGA, both strands.
//
// Build: nvcc -O3 -o enumerate_gpu enumerate.cu
// CLI matches the CPU binary so the benchmark can diff the TSV.

#include <cuda_runtime.h>

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

struct DeviceHit {
    int pos;
    int mm;
    int strand;  // 0 plus, 1 minus
    int pam_index;
    char site[24];
};

static char complement_host(char b) {
    switch (b) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default: return 'N';
    }
}

__device__ char dcomp(char b) {
    if (b == 'A') return 'T';
    if (b == 'C') return 'G';
    if (b == 'G') return 'C';
    if (b == 'T') return 'A';
    return 'N';
}

__device__ bool pam_ok(const char* pam, const char* cls) {
    for (int i = 0; i < 3; ++i) {
        if (cls[i] != 'N' && cls[i] != pam[i]) return false;
    }
    return true;
}

__global__ void enumerate_kernel(const char* seq, int n, const char* guide,
                                 const char* pams, int n_pams, int max_mm,
                                 DeviceHit* hits, int capacity, int* count) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int last = n - 22;
    if (i >= last) return;
    char window[23];
    bool clean = true;
    for (int k = 0; k < 23; ++k) {
        window[k] = seq[i + k];
        if (window[k] != 'A' && window[k] != 'C' && window[k] != 'G' && window[k] != 'T') clean = false;
    }
    if (!clean) return;

    for (int strand = 0; strand < 2; ++strand) {
        char site[23];
        if (strand == 0) {
            for (int k = 0; k < 23; ++k) site[k] = window[k];
        } else {
            for (int k = 0; k < 23; ++k) site[k] = dcomp(window[22 - k]);
        }
        int pam_index = -1;
        for (int p = 0; p < n_pams; ++p) {
            if (pam_ok(site + 20, pams + p * 3)) {
                pam_index = p;
                break;
            }
        }
        if (pam_index < 0) continue;
        int mm = 0;
        for (int k = 0; k < 20; ++k) mm += guide[k] != site[k];
        if (mm > max_mm) continue;
        int slot = atomicAdd(count, 1);
        if (slot >= capacity) continue;
        DeviceHit hit;
        hit.pos = i + 1;
        hit.mm = mm;
        hit.strand = strand;
        hit.pam_index = pam_index;
        for (int k = 0; k < 23; ++k) hit.site[k] = site[k];
        hit.site[23] = 0;
        hits[slot] = hit;
    }
}

static std::vector<std::pair<std::string, std::string>> load_fasta(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    std::vector<std::pair<std::string, std::string>> chroms;
    std::string line, name, seq;
    auto flush = [&]() {
        if (!name.empty()) chroms.emplace_back(name, std::move(seq));
        seq.clear();
    };
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            flush();
            auto id = line.substr(1);
            auto sp = id.find_first_of(" \t");
            name = sp == std::string::npos ? id : id.substr(0, sp);
        } else {
            for (char c : line) {
                if (c == '\r') continue;
                seq.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(c))));
            }
        }
    }
    flush();
    return chroms;
}

int main(int argc, char** argv) {
    std::string fasta, guide, pam_arg = "NGG";
    int max_mm = 4;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        auto need = [&]() { return std::string(argv[++i]); };
        if (arg == "--fasta") fasta = need();
        else if (arg == "--guide") guide = need();
        else if (arg == "--pam") pam_arg = need();
        else if (arg == "--max-mm") max_mm = std::stoi(need());
        else if (arg == "--algo" || arg == "--threads") { need(); }
        else { std::cerr << "unknown arg " << arg << "\n"; return 2; }
    }
    if (fasta.empty() || guide.size() != 20) {
        std::cerr << "usage: enumerate_gpu --fasta F --guide 20mer --pam NGG,NAG --max-mm 4\n";
        return 2;
    }
    std::vector<std::string> pams;
    std::stringstream ss(pam_arg);
    std::string item;
    while (std::getline(ss, item, ',')) if (!item.empty()) pams.push_back(item);
    std::string pamflat;
    for (const auto& p : pams) pamflat += p;

    auto chroms = load_fasta(fasta);
    const int capacity = 200000;
    std::cout << "chrom\tpos\tstrand\tmm\tpam_class\tsite23\n";
    struct Row { std::string chrom; int pos; char strand; int mm; std::string pam; std::string site; };
    std::vector<Row> rows;

    for (const auto& [chrom, seq] : chroms) {
        if (seq.size() < 23) continue;
        char* d_seq = nullptr;
        char* d_guide = nullptr;
        char* d_pams = nullptr;
        DeviceHit* d_hits = nullptr;
        int* d_count = nullptr;
        cudaMalloc(&d_seq, seq.size());
        cudaMalloc(&d_guide, 20);
        cudaMalloc(&d_pams, pamflat.size());
        cudaMalloc(&d_hits, sizeof(DeviceHit) * capacity);
        cudaMalloc(&d_count, sizeof(int));
        cudaMemcpy(d_seq, seq.data(), seq.size(), cudaMemcpyHostToDevice);
        cudaMemcpy(d_guide, guide.data(), 20, cudaMemcpyHostToDevice);
        cudaMemcpy(d_pams, pamflat.data(), pamflat.size(), cudaMemcpyHostToDevice);
        cudaMemset(d_count, 0, sizeof(int));
        int threads = 256;
        int blocks = (static_cast<int>(seq.size()) + threads - 1) / threads;
        enumerate_kernel<<<blocks, threads>>>(d_seq, static_cast<int>(seq.size()), d_guide, d_pams,
                                              static_cast<int>(pams.size()), max_mm, d_hits, capacity, d_count);
        cudaDeviceSynchronize();
        int count = 0;
        cudaMemcpy(&count, d_count, sizeof(int), cudaMemcpyDeviceToHost);
        int keep = std::min(count, capacity);
        std::vector<DeviceHit> hits(keep);
        if (keep) cudaMemcpy(hits.data(), d_hits, sizeof(DeviceHit) * keep, cudaMemcpyDeviceToHost);
        for (const auto& hit : hits) {
            rows.push_back(Row{chrom, hit.pos, hit.strand ? '-' : '+', hit.mm, pams[hit.pam_index],
                               std::string(hit.site, hit.site + 23)});
        }
        cudaFree(d_seq);
        cudaFree(d_guide);
        cudaFree(d_pams);
        cudaFree(d_hits);
        cudaFree(d_count);
        if (count > capacity) {
            std::cerr << "hit capacity exceeded on " << chrom << "\n";
            return 3;
        }
    }
    std::sort(rows.begin(), rows.end(), [](const Row& a, const Row& b) {
        if (a.chrom != b.chrom) return a.chrom < b.chrom;
        if (a.pos != b.pos) return a.pos < b.pos;
        if (a.strand != b.strand) return a.strand < b.strand;
        return a.pam < b.pam;
    });
    for (const auto& row : rows) {
        std::cout << row.chrom << '\t' << row.pos << '\t' << row.strand << '\t' << row.mm << '\t'
                  << row.pam << '\t' << row.site << '\n';
    }
    return 0;
}
