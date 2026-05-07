/**
 * FastaReader.h - FASTA file parser
 * 
 * Reads single or multi-FASTA files and extracts species names + sequences.
 */

#ifndef FASTA_READER_H
#define FASTA_READER_H

#include <string>
#include <vector>
#include <utility>

namespace mlmaws {

/**
 * Represents a single FASTA record (species name + DNA sequence).
 */
struct FastaRecord {
    std::string name;     // species/taxon name (from header line)
    std::string sequence; // concatenated DNA sequence (uppercase, no whitespace)
};

/**
 * Reads FASTA files. Supports both single-species and multi-species files.
 */
class FastaReader {
public:
    /**
     * Read a single FASTA file containing one or more sequences.
     * Returns vector of FastaRecord.
     */
    static std::vector<FastaRecord> readFile(const std::string& filepath);

    /**
     * Read all FASTA files from a directory.
     * Each file is treated as containing one species.
     * Returns vector of FastaRecord (one per file).
     */
    static std::vector<FastaRecord> readDirectory(const std::string& dirpath);

    /**
     * Clean a DNA sequence: uppercase, remove non-ACGT characters.
     */
    static std::string cleanSequence(const std::string& raw);
};

} // namespace mlmaws

#endif // FASTA_READER_H
