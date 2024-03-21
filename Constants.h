#ifndef CONSTANTS_H
#define CONSTANTS_H

namespace constants{
constexpr double epsilon = 1e-3;
}

namespace fileIO{
const std::string inFileName = "/data/mciacco/matchingqa/AnalysisResultsLHC15o.root"; // AnalysisResults.root
const std::string outFileName = "out_run2_evecuts.root";
}

namespace bins{
// small
constexpr int bins_small[3]{10, 100, 4};
constexpr double xmin_small[3]{0., 0., 0.};
constexpr double xmax_small[3]{10., 100., 0.8};

constexpr int bins[3]{10, 10, 4};
constexpr double xmin[3]{0., 0., 0.};
constexpr double xmax[3]{10., 100., 0.8};
}

#endif // CONSTANTS_H
