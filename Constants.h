#ifndef CONSTANTS_H
#define CONSTANTS_H

constexpr bool kApplyEffCorrection = true;
const char* genRec = "";
const char* cutSets[] = {"", "_cpa0", "_cpa1", "_mass0", "_mass1", "_tpccls0", "_tpccls1", "_pid0", "_pid1", "_dca0"};

namespace constants{
constexpr double epsilon = 1e-3;
}

namespace fileIO{
const std::string inFileName = "./AnalysisResultsLHC15o18qrGridJob.root"; // AnalysisResults.root
const std::string outFileName = "out_run2_LHC15o18qr_grid_job.root";
}

namespace bins{
// small
constexpr int bins_small[3]{10, 100, 1};
constexpr double xmin_small[3]{0., 0., 0.};
constexpr double xmax_small[3]{10., 100., 0.8};

constexpr int bins[3]{10, 10, 1};
constexpr double xmin[3]{0., 0., 0.};
constexpr double xmax[3]{10., 100., 0.8};
}

#endif // CONSTANTS_H
