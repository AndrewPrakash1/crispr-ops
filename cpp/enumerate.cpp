// Bit-exact off-target enumeration for SpCas9 spacers.
// Two algorithms (full scan and PAM-first) share one window function so a
// disagreement between them is a bug, not a hidden filter.
//
// Coordinates: pos is the 1-based start of the 23 bp + strand window
// (protospacer then PAM on +, or PAM-side then protospacer on -).
// Mismatches are counted on the 20 bp protospacer only.

#include <algorithm>
#include <cctype>
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>
#include <utility>
#include <vector>

struct Hit {
    std::string chrom;
    int pos = 0;
    char strand = '+';
    int mm = 0;
    std::string site23;
    std::string pam_class;
};

static char complement(char b) {
    switch (b) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default: return 'N';
    }
}

static std::string revcomp(const std::string& seq) {
    std::string out(seq.size(), 'N');
    for (size_t i = 0; i < seq.size(); ++i) {
        out[seq.size() - 1 - i] = complement(seq[i]);
    }
    return out;
}

static bool acgt23(const std::string& w) {
    for (char b : w) {
        if (b != 'A' && b != 'C' && b != 'G' && b != 'T') return false;
    }
    return true;
}

static int hamming20(const std::string& guide, const char* site) {
    int mm = 0;
    for (int i = 0; i < 20; ++i) mm += guide[i] != site[i];
    return mm;
}

static bool pam_ok(const char* pam, const std::string& cls) {
    for (int i = 0; i < 3; ++i) {
        if (cls[i] != 'N' && cls[i] != pam[i]) return false;
    }
    return true;
}

static bool any_pam(const char* pam, const std::vector<std::string>& classes) {
    for (const auto& cls : classes) {
        if (pam_ok(pam, cls)) return true;
    }
    return false;
}

static int class_index(const char* pam, const std::vector<std::string>& classes) {
    for (int i = 0; i < static_cast<int>(classes.size()); ++i) {
        if (pam_ok(pam, classes[i])) return i;
    }
    return -1;
}

struct Task {
    std::string fasta;
    std::string guide;
    std::vector<std::string> pams;
    int max_mm = 4;
    std::string algo = "pam-skip";
    int threads = 1;
};

static bool pam_region_possible(const std::string& seq, int i, const std::vector<std::string>& pams) {
    char forward[3] = {seq[i + 20], seq[i + 21], seq[i + 22]};
    if (any_pam(forward, pams)) return true;
    char reverse[3] = {complement(seq[i + 2]), complement(seq[i + 1]), complement(seq[i])};
    return any_pam(reverse, pams);
}

static bool is_acgt(char b) {
    return b == 'A' || b == 'C' || b == 'G' || b == 'T';
}

static void record_hit(std::vector<Hit>& out, const std::string& chrom, int pos, char strand,
                       int mm, const char* site, const std::string& pam_class) {
    out.push_back(Hit{chrom, pos, strand, mm, std::string(site, 23), pam_class});
}

static void consider_at(const std::string& chrom, const std::string& seq, int i,
                        const Task& task, std::vector<Hit>& out) {
    const char* guide = task.guide.data();
    char fpam[3] = {seq[i + 20], seq[i + 21], seq[i + 22]};
    if (is_acgt(fpam[0]) && is_acgt(fpam[1]) && is_acgt(fpam[2])) {
        int which = class_index(fpam, task.pams);
        if (which >= 0) {
            int mm = 0;
            bool ok = true;
            for (int k = 0; k < 20; ++k) {
                char b = seq[i + k];
                if (!is_acgt(b)) { ok = false; break; }
                mm += b != guide[k];
                if (mm > task.max_mm) { ok = false; break; }
            }
            if (ok) {
                char site[23];
                for (int k = 0; k < 20; ++k) site[k] = seq[i + k];
                site[20] = fpam[0]; site[21] = fpam[1]; site[22] = fpam[2];
                record_hit(out, chrom, i + 1, '+', mm, site, task.pams[which]);
            }
        }
    }
    if (is_acgt(seq[i]) && is_acgt(seq[i + 1]) && is_acgt(seq[i + 2])) {
        char rpam[3] = {complement(seq[i + 2]), complement(seq[i + 1]), complement(seq[i])};
        int which = class_index(rpam, task.pams);
        if (which >= 0) {
            int mm = 0;
            bool ok = true;
            for (int k = 0; k < 20; ++k) {
                char raw = seq[i + 22 - k];
                if (!is_acgt(raw)) { ok = false; break; }
                mm += complement(raw) != guide[k];
                if (mm > task.max_mm) { ok = false; break; }
            }
            if (ok) {
                char site[23];
                for (int k = 0; k < 20; ++k) site[k] = complement(seq[i + 22 - k]);
                site[20] = rpam[0]; site[21] = rpam[1]; site[22] = rpam[2];
                record_hit(out, chrom, i + 1, '-', mm, site, task.pams[which]);
            }
        }
    }
}

static std::vector<std::pair<std::string, std::string>> load_fasta(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("cannot open " + path);
    in.seekg(0, std::ios::end);
    auto bytes = static_cast<size_t>(in.tellg());
    in.seekg(0);
    std::string file(bytes, '\0');
    in.read(file.data(), static_cast<std::streamsize>(bytes));
    std::vector<std::pair<std::string, std::string>> chroms;
    std::string name, seq;
    seq.reserve(bytes);
    auto flush = [&]() {
        if (!name.empty()) chroms.emplace_back(name, std::move(seq));
        seq.clear();
        seq.reserve(bytes);
    };
    size_t i = 0;
    while (i < file.size()) {
        if (file[i] == '>') {
            flush();
            size_t end = file.find('\n', i);
            if (end == std::string::npos) end = file.size();
            auto id = file.substr(i + 1, end - (i + 1));
            if (!id.empty() && id.back() == '\r') id.pop_back();
            auto sp = id.find_first_of(" \t");
            name = sp == std::string::npos ? id : id.substr(0, sp);
            i = end + 1;
            continue;
        }
        size_t end = file.find('\n', i);
        if (end == std::string::npos) end = file.size();
        for (size_t j = i; j < end; ++j) {
            unsigned char c = static_cast<unsigned char>(file[j]);
            if (c == '\r') continue;
            seq.push_back(static_cast<char>(std::toupper(c)));
        }
        i = end + 1;
    }
    flush();
    return chroms;
}

static void scan_range(const std::string& chrom, const std::string& seq, int begin, int end,
                       const Task& task, std::vector<Hit>& out) {
    bool ngg = false, nag = false, nga = false, other = false;
    for (const auto& cls : task.pams) {
        if (cls == "NGG") ngg = true;
        else if (cls == "NAG") nag = true;
        else if (cls == "NGA") nga = true;
        else other = true;
    }
    const bool standard = task.algo == "pam-skip" && !other;
    const char* s = seq.data();
    for (int i = begin; i < end; ++i) {
        if (standard) {
            const char b0 = s[i], b1 = s[i + 1], b21 = s[i + 21], b22 = s[i + 22];
            const bool forward = (ngg && b21 == 'G' && b22 == 'G') || (nag && b21 == 'A' && b22 == 'G') ||
                                 (nga && b21 == 'G' && b22 == 'A');
            const bool reverse = (ngg && b0 == 'C' && b1 == 'C') || (nag && b0 == 'C' && b1 == 'T') ||
                                 (nga && b0 == 'T' && b1 == 'C');
            if (!forward && !reverse) continue;
        }
        consider_at(chrom, seq, i, task, out);
    }
}

static std::vector<Hit> enumerate_chrom(const std::string& chrom, const std::string& seq, const Task& task) {
    if (seq.size() < 23) return {};
    int last = static_cast<int>(seq.size()) - 22;
    int workers = std::max(1, task.threads);
    if (workers == 1 || last < 100000) {
        std::vector<Hit> hits;
        scan_range(chrom, seq, 0, last, task, hits);
        return hits;
    }
    std::vector<std::vector<Hit>> parts(workers);
    std::vector<std::thread> pool;
    int chunk = (last + workers - 1) / workers;
    for (int t = 0; t < workers; ++t) {
        int begin = t * chunk;
        int end = std::min(last, begin + chunk);
        if (begin >= end) continue;
        pool.emplace_back([&, t, begin, end]() {
            scan_range(chrom, seq, begin, end, task, parts[t]);
        });
    }
    for (auto& th : pool) th.join();
    std::vector<Hit> hits;
    for (auto& part : parts) {
        hits.insert(hits.end(), part.begin(), part.end());
    }
    return hits;
}

static void usage() {
    std::cerr << "usage: enumerate_cpu --fasta F --guide 20mer --pam NGG,NAG --max-mm 4 "
                 "[--algo pam-skip|scan] [--threads N]\n";
}

int main(int argc, char** argv) {
    Task task;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        auto need = [&](const std::string& flag) {
            if (i + 1 >= argc) throw std::runtime_error("missing value for " + flag);
            return std::string(argv[++i]);
        };
        if (arg == "--fasta") task.fasta = need(arg);
        else if (arg == "--guide") task.guide = need(arg);
        else if (arg == "--pam") {
            std::stringstream ss(need(arg));
            std::string item;
            while (std::getline(ss, item, ',')) {
                if (!item.empty()) task.pams.push_back(item);
            }
        } else if (arg == "--max-mm") task.max_mm = std::stoi(need(arg));
        else if (arg == "--algo") task.algo = need(arg);
        else if (arg == "--threads") task.threads = std::stoi(need(arg));
        else if (arg == "--help") { usage(); return 0; }
        else { usage(); return 2; }
    }
    if (task.fasta.empty() || task.guide.size() != 20 || task.pams.empty()) {
        usage();
        return 2;
    }
    for (char b : task.guide) {
        if (b != 'A' && b != 'C' && b != 'G' && b != 'T') {
            std::cerr << "guide must be ACGT\n";
            return 2;
        }
    }
    using clock = std::chrono::steady_clock;
    auto t0 = clock::now();
    auto chroms = load_fasta(task.fasta);
    auto t1 = clock::now();
    std::vector<Hit> hits;
    for (const auto& [chrom, seq] : chroms) {
        auto part = enumerate_chrom(chrom, seq, task);
        hits.insert(hits.end(), part.begin(), part.end());
    }
    auto t2 = clock::now();
    std::cerr << "load_s\t" << std::chrono::duration<double>(t1 - t0).count()
              << "\tscan_s\t" << std::chrono::duration<double>(t2 - t1).count() << "\n";
    std::sort(hits.begin(), hits.end(), [](const Hit& a, const Hit& b) {
        if (a.chrom != b.chrom) return a.chrom < b.chrom;
        if (a.pos != b.pos) return a.pos < b.pos;
        if (a.strand != b.strand) return a.strand < b.strand;
        return a.pam_class < b.pam_class;
    });
    std::cout << "chrom\tpos\tstrand\tmm\tpam_class\tsite23\n";
    for (const auto& hit : hits) {
        std::cout << hit.chrom << '\t' << hit.pos << '\t' << hit.strand << '\t'
                  << hit.mm << '\t' << hit.pam_class << '\t' << hit.site23 << '\n';
    }
    return 0;
}
