#ifndef CONSTANTS_H
#define CONSTANTS_H

constexpr bool kApplyEffCorrection = false;

namespace constants{
constexpr double epsilon = 1e-3;
}

namespace fileIO{
const std::string inFileName = "./AnalysisResults_LHC15o_20240420.root"; // AnalysisResults.root
const std::string outFileName = "out_run2_LHC15o_20240420_noEff.root";
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
