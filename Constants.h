#ifndef CONSTANTS_H
#define CONSTANTS_H

constexpr bool kApplyEffCorrection = true;
const char* genRec = "";
const char* cutSets[] = {"", "_cpa0", "_cpa1", "_dcadaugh0", "_dcadaugh1", "_dcatopv0", "_dcatopv1", "_mass0", "_mass1", "_tpccls0", "_tpccls1", "_pidtof0", "_pidtof1", "_pidtpc0", "_pidtpc1", "_dca0", "_zvtx0"};
const char* period = "18qr";

namespace constants{
constexpr double epsilon = 1e-3;
}

namespace fileIO{
const std::string inFileName = "./AnalysisResultsLHC18qrSysHyperloop2.root"; // AnalysisResults.root
const std::string inFileNameMC0 = "./AnalysisResultsLHC20g7abcSysHyperloop.root"; // AnalysisResults.root
const std::string inFileNameMC1 = "./AnalysisResultsLHC20e3abcSysHyperloop.root"; // AnalysisResults.root
const std::string outFileName = "out_run2_LHC18qr_grid_sys_hyperloop3.root";
}

namespace bins{
// small
constexpr int bins_small[3]{10, 100, 1};
constexpr double xmin_small[3]{0., 0., 0.};
constexpr double xmax_small[3]{10., 100., 0.8};

constexpr int bins[3]{10, 1, 1};
constexpr double xmin[3]{0., 0., 0.};
constexpr double xmax[3]{10., 80., 0.8};
}

#endif // CONSTANTS_H
