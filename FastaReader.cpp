/**
 * FastaReader.cpp - FASTA file parser implementation (C++17)
 */

#include "FastaReader.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cctype>
#include <filesystem>

namespace fs = std::filesystem;

namespace mlmaws {

std::string FastaReader::cleanSequence(const std::string& raw) {
    std::string clean;
    clean.reserve(raw.size());
    for (char c : raw) {
        char upper = std::toupper(c);
        if (upper == 'A' || upper == 'C' || upper == 'G' || upper == 'T') {
            clean.push_back(upper);
        }
    }
    return clean;
}

std::vector<FastaRecord> FastaReader::readFile(const std::string& filepath) {
    std::vector<FastaRecord> records;
    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "[ERROR] Cannot open FASTA file: " << filepath << std::endl;
        return records;
    }

    std::string line;
    FastaRecord current;
    bool hasRecord = false;

    while (std::getline(file, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;

        if (line[0] == '>') {
            if (hasRecord && !current.sequence.empty()) {
                current.sequence = cleanSequence(current.sequence);
                records.push_back(std::move(current));
                current = FastaRecord();
            }
            current.name = line.substr(1);
            size_t spacePos = current.name.find_first_of(" \t");
            if (spacePos != std::string::npos) {
                current.name = current.name.substr(0, spacePos);
            }
            hasRecord = true;
        } else {
            current.sequence += line;
        }
    }

    if (hasRecord && !current.sequence.empty()) {
        current.sequence = cleanSequence(current.sequence);
        records.push_back(std::move(current));
    }

    file.close();
    return records;
}

std::vector<FastaRecord> FastaReader::readDirectory(const std::string& dirpath) {
    std::vector<FastaRecord> allRecords;

    if (!fs::exists(dirpath) || !fs::is_directory(dirpath)) {
        std::cerr << "[ERROR] Directory not found: " << dirpath << std::endl;
        return allRecords;
    }

    std::vector<fs::path> files;
    for (const auto& entry : fs::directory_iterator(dirpath)) {
        if (entry.is_regular_file()) {
            std::string ext = entry.path().extension().string();
            std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
            if (ext == ".fasta" || ext == ".fa" || ext == ".fna" || ext == ".fas") {
                files.push_back(entry.path());
            }
        }
    }
    std::sort(files.begin(), files.end());

    for (const auto& fpath : files) {
        std::vector<FastaRecord> recs = readFile(fpath.string());
        if (recs.empty()) continue;

        if (recs.size() == 1) {
            allRecords.push_back(std::move(recs[0]));
        } else {
            FastaRecord merged;
            merged.name = fpath.stem().string();
            for (auto& r : recs) {
                merged.sequence += r.sequence;
            }
            allRecords.push_back(std::move(merged));
        }
    }

    std::cerr << "[INFO] Read " << allRecords.size()
              << " species from directory: " << dirpath << std::endl;
    return allRecords;
}

} // namespace mlmaws
