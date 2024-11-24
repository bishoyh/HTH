
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <future>
#include <limits>
#include <cmath>
#include <algorithm>
#include <cctype>
#include <mutex>

/*
Based on the code from C Code for Helix-Turn-Helix Motif Identification by  Conrad Halling 
And the paper Dodd, I. B., and J. B. Egan. 1990.  Improved detection of helix-turn-helix DNA-binding motifs in protein sequences. Nucleic Acids Res.18:5019-5026.
Added multithreading, C++ features, and error handling, and improved output.

-BSK 2024
*/

constexpr int WINDOW_SIZE = 22;
constexpr double NON_HTH_MEAN_SCORE = 238.71;
constexpr double NON_HTH_STD_DEV = 293.61;
constexpr int AMINO_ACIDS_COUNT = 20;

const int weightMatrix[AMINO_ACIDS_COUNT][WINDOW_SIZE] = {
    /* A (alanine) */
    {-125, -194,  -84,   70,   36,   54,  238,  -15,   77,   26, -194,
     -194,  -56,  -84,   14,   77,  -56,  -56,  -56,   46, -195,   36},
    /* C (cysteine) */
    { -64,  -64,  -63,  -63,  -64,  -64,  -64,  -64,  -64,   47,   47,
      -63,  -63,  -64,  -64,  -64,  -64,  -64,  -64,  -63,  -64,   47},
    /* D (aspartate) */
    {-156, -154, -156, -154,  109, -156, -156,  109, -154, -156,    6,
     -156, -154,  -85, -156, -156, -156, -156, -154, -154, -156,  -85},
    /* E (glutamate) */
    { -31,   -9, -171,   70,  156, -171, -171,  107,   50,  -60,  -60,
     -171,  -60,   78,   86, -171, -171, -101, -171, -170,   86,    9},
    /* F (phenylalanine) */
    {  10, -130,   10, -130, -130,   10, -130, -129, -130,  102, -130,
     -130, -130, -130, -129, -130, -130, -129, -129, -129,  180, -130},
    /* G (glycine) */
    {  30,    5, -190,  -51, -191, -191,   18, -191, -191, -191,  202,
     -191,  -10, -191,    5, -190, -191,  -80, -190, -190, -191,  -51},
    /* H (histidine ) */
    {  62,   33,  -76,  -76,   -7,  -78,  -78,   33,   -7,  -78,   84,
      -78,   33,   33,  -78,   -7,  -78,   -7,   62,   84,  -78,   -7},
    /* I (isoleucine ) */
    {  75, -156,  101,  -45,  -86,  116, -156,  -16,   65,  -16, -156,
     128, -156,  -86, -156, -155,  188, -155,  -16,   53,  122, -155},
    /* K (lysine ) */
    { -31,  -31,   10,   70,   79, -170, -170,   94,   70, -171,   -9,
     -100, -100,   25, -100, -170, -171,   -9,   38,  -31,   -9,  101},
    /* L (leucine) */
    {  66, -212,   72, -213, -212,  144, -213, -102,   37,  132, -213,
      97, -213, -142, -212, -212,   97, -212, -212,   37,   88, -213},
    /* M (methionine) */
    { 122,  -74,   -3,  -73,  -73,  -73,  -74,  -74,   88,  122,  -73,
     158,  -74,  -74,   -3,  -74,   -3,  -74,  -74,   -3,  -74,  -73},
    /* N (asparagine) */
    {-137,   72, -137, -136, -137, -137, -137,  -67, -136, -136,  128,
     -137,   72, -136,    2,  -67, -137,    2,   84, -137, -137,  104},
    /* P (proline) */
    {-156,   23, -157, -156, -157, -157, -157, -157, -157, -157, -157,
     -157,  -46,  101,   39, -157, -157, -157, -157, -157, -157,  -46},
    /* Q (glutamine) */
    { -60, -130,  175,   90,  110, -131, -131,   90,   78, -131, -131,
     -60, -130,  154,   65,  119, -131,  -20,  119,   31,  -20,   90},
    /* R (arginine) */
    {  65,   76,  110,   65,    7, -155, -154,  123,   76, -155, -154,
     -155, -154,  129,   54,   40, -155,  129,  179,  -45, -155,  123},
    /* S (serine) */
    {-118,   96, -188,   21,  -48, -187,   -8, -187, -118, -118,  -77,
     -188,  174, -187,  135,  -26, -188,  150,  -77, -188, -187,  -26},
    /* T (threonine) */
    {  11,  149,  -59,   80, -169,   -8, -170,  -99,  -99,  -30, -170,
      -8,  131,  -30,  -59,  198, -170,  -30, -169, -170,  -30,  -59},
    /* V (valine) */
    {  17,  -67, -177, -177, -108,  100, -178, -178, -108,   71, -178,
     160, -178,  -16,  -67,  -67,  169, -178,   17,   31,   17, -178},
    /* W (tryptophan) */
    {  44,  -26,  -26,  -26,  -26,  -26,  -26,  -26,  -26,  -26,  -26,
     -26,  -26,  -25,  -26,  -25,  -26,  -26,   44,  279,  -26,  -26},
    /* Y (tyrosine) */
    { -40, -110,   30,    1, -110, -109, -110, -110,  -40,   30, -110,
     -109, -109,  -40, -110,  -40, -110,  162,   52,   86, -110, -110}
};

std::unordered_map<char, int> createAminoAcidMap() {
    std::unordered_map<char, int> aminoAcidMap;
    const std::string aminoAcids = "ACDEFGHIKLMNPQRSTVWY";
    for (int i = 0; i < aminoAcids.size(); ++i) {
        char aa = aminoAcids[i];
        aminoAcidMap[aa] = i;
        aminoAcidMap[std::tolower(aa)] = i;
    }
    return aminoAcidMap;
}

struct Sequence {
    std::string id;
    std::string sequence;
};

struct Result {
    std::string id;
    double convertedScore = 0.0;
    size_t maxScorePosition = 0; // 1-based index
    std::string maxScoreSequence;
    std::string interpretation;
    std::string errorMessage; // In case of errors
};

std::vector<Sequence> readSequences(const std::string& filename) {
    std::vector<Sequence> sequences;
    std::ifstream infile(filename);
    if (!infile) {
        throw std::runtime_error("Could not open input file.");
    }
    std::string line;
    Sequence seq;
    while (std::getline(infile, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!seq.id.empty()) {
                sequences.push_back(seq);
                seq = Sequence();
            }
            seq.id = line.substr(1); // Remove '>'
        } else {
            // Remove whitespace and append to sequence
            line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());
            seq.sequence += line;
        }
    }
    // Add the last sequence
    if (!seq.id.empty()) {
        sequences.push_back(seq);
    }
    return sequences;
}

Result processSequence(const Sequence& seq, const std::unordered_map<char, int>& aminoAcidMap) {
    Result result;
    result.id = seq.id;
    const std::string& sequence = seq.sequence;
    size_t sequenceLength = sequence.size();
    if (sequenceLength < WINDOW_SIZE) {
        result.errorMessage = "Sequence too short to analyze.";
        return result;
    }

    int maxScore = std::numeric_limits<int>::min();
    size_t maxScorePosition = 0;
    double convertedScore = 0.0;
    bool invalidCharacter = false;

    size_t maxWindowPosition = sequenceLength - WINDOW_SIZE + 1;

    for (size_t i = 0; i < maxWindowPosition; ++i) {
        int tempScore = 0;
        for (int j = 0; j < WINDOW_SIZE; ++j) {
            char residue = sequence[i + j];
            auto it = aminoAcidMap.find(residue);
            if (it == aminoAcidMap.end()) {
                result.errorMessage = "Invalid character in sequence.";
                invalidCharacter = true;
                break;
            }
            int aminoAcidIndex = it->second;
            tempScore += weightMatrix[aminoAcidIndex][j];
        }
        if (invalidCharacter) {
            break;
        }
        if (tempScore > maxScore) {
            maxScore = tempScore;
            maxScorePosition = i;
            convertedScore = (maxScore - NON_HTH_MEAN_SCORE) / NON_HTH_STD_DEV;
        }
    }
    if (invalidCharacter) {
        return result;
    }

    // Get the max score sequence
    result.maxScoreSequence = sequence.substr(maxScorePosition, WINDOW_SIZE);
    result.maxScorePosition = maxScorePosition + 1; // 1-based index
    result.convertedScore = convertedScore;

    // Determine interpretation
    int percentage = 0;
    if (convertedScore >= 4.5)
        percentage = 100;
    else if (convertedScore >= 4.0)
        percentage = 90;
    else if (convertedScore >= 3.5)
        percentage = 71;
    else if (convertedScore >= 3.0)
        percentage = 50;
    else if (convertedScore >= 2.5)
        percentage = 25;

    if (convertedScore < 2.5) {
        result.interpretation = "Not significant.";
    } else {
        result.interpretation = "Approximately " + std::to_string(percentage) + "% probability that this protein contains a helix-turn-helix motif.";
    }

    return result;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: hth <input_fasta_file> <output_tsv_file>" << std::endl;
        return 1;
    }
    std::string inputFile = argv[1];
    std::string outputFile = argv[2];

    // Create amino acid mapping
    auto aminoAcidMap = createAminoAcidMap();

    // Read sequences
    std::vector<Sequence> sequences;
    try {
        sequences = readSequences(inputFile);
    } catch (const std::exception& ex) {
        std::cerr << "Error reading sequences: " << ex.what() << std::endl;
        return 1;
    }

    // Process sequences in parallel
    std::vector<std::future<Result>> futures;
    for (const auto& seq : sequences) {
        futures.push_back(std::async(std::launch::async, processSequence, seq, std::ref(aminoAcidMap)));
    }

    // Collect results
    std::vector<Result> results;
    for (auto& f : futures) {
        results.push_back(f.get());
    }

    // Write results to TSV file
    std::ofstream outfile(outputFile);
    if (!outfile) {
        std::cerr << "Could not open output file for writing." << std::endl;
        return 1;
    }

    // Write header
    outfile << "Sequence_ID\tConverted_Score\tMax_Score_Position\tMax_Score_Sequence\tInterpretation\tError" << std::endl;

    for (const auto& res : results) {
        outfile << res.id << '\t';
        if (!res.errorMessage.empty()) {
            outfile << "\t\t\t\t" << res.errorMessage << std::endl;
        } else {
            outfile << res.convertedScore << '\t';
            outfile << res.maxScorePosition << '\t';
            outfile << res.maxScoreSequence << '\t';
            outfile << res.interpretation << '\t' << std::endl;
        }
    }

    return 0;
}
