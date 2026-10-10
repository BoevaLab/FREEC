#include "GenomeCopyNumber.h"
#include "SeekSubclones.h"
#include <cassert>
#include <fstream>
#include <string>

// Globals normally supplied by the command-line entry point.
int verbose = 0;
double minMappabilityPerWindow = 0.85;
bool uniqueMatch = false;
char *ULetters[90] = {};

int main(int argc, char** argv) {
    assert(argc == 3);
    const std::string stem = argv[1];
    const bool segmented = std::string(argv[2]) == "segmented";
    const int bins = segmented ? 200 : 30;
    {
        std::ofstream counts((stem + ".cnp").c_str());
        for (int i = 0; i < bins; ++i)
            counts << "1\t" << i * 1000 << "\t" << (segmented && i >= bins / 2 ? 1.5 : 1.0) << "\n";
    }
    ThreadPoolManager::init(1);
    GenomeCopyNumber sample;
    sample.setWESanalysis(false);
    sample.setIfLogged(false);
    sample.readCopyNumber(stem + ".cnp");
    sample.setPloidy(2);
    ChrCopyNumber& chr = sample.getChrCopyNumberAt(0);
    chr.fillInRatio(false);
    if (segmented)
        chr.calculateBreakpoints(0.8, bins, 2);
    chr.calculateCopyNumberMedian(2, 1, false, false, false);
    assert(chr.getLength() == bins);
    assert(segmented ? !chr.getBreakPoints().empty() : chr.getBreakPoints().empty());
    const std::vector<float> original = chr.getRatio();
    sample.setSeekSubclones(true);
    SeekSubclones subclones(sample, 2, stem, 20);
    assert(chr.getRatio() == original);
    assert(chr.getLength() == bins);
    std::cout << "Subclone evaluation completed for " << (segmented ? "multiple segments" : "one segment") << "\n";
}
