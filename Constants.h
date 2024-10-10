#ifndef CONSTANTS_H
#define CONSTANTS_H

constexpr bool kApplyEffCorrection = true;
const char* genRec = "";
//const char* cutSets[] = {"", "_cpa0", "_cpa1", "_dcadaugh0", "_dcadaugh1", "_dcatopv0", "_dcatopv1", "_mass0", "_mass1", "_tpccls0", "_tpccls1", "_pidtof0", "_pidtof1", "_pidtpc0", "_pidtpc1", "_dca0", "_zvtx0"};
const char* cutSetsName[] = {"", "cos#theta_{p}^{low}", "cos#theta_{p}^{up}", "DCA_{tracks}^{low}", "DCA_{tracks}^{up}", "DCA_{V^{0}PV}^{low}", "DCA_{V^{0}PV}^{up}", "#it{M}^{low}", "#it{M}^{up}", "#it{n}_{TPCCls}^{low}", "#it{n}_{TPCCls}^{up}", "n#sigma_{TOF}^{low}", "n#sigma_{TOF}^{up}", "n#sigma_{TPC}^{low}", "n#sigma_{TPC}^{up}", "DCA^{low}", "V_{z}^{up}"};

//const char* cutSets[] = {"", "_sys1", "_sys2", "_sys3", "_sys4", "_sys5", "_sys6", "_sys7", "_sys8", "_sys9", "_sys10", "_sys11", "_sys12", "_sys13", "_sys14", "_sys15", "_sys16", "_sys17", "_sys18", "_sys19", "_sys20", "_sys21", "_sys22", "_sys23", "_sys24", "_sys25", "_sys26", "_sys27", "_sys28", "_sys29", "_sys30"};
const char* cutSets[] = {"", "_sys01", "_sys02", "_sys03", "_sys04", "_sys05", "_sys06", "_sys07", "_sys08", "_sys09", "_sys10", "_sys11", "_sys12", "_sys13", "_sys14", "_sys15", "_sys16", "_sys17", "_sys18", "_sys19", "_sys20", "_sys21", "_sys22", "_sys23", "_sys24", "_sys25", "_sys26", "_sys27", "_sys28", "_sys29", "_sys30", "_sys31", "_sys32", "_sys33", "_sys34", "_sys35", "_sys36", "_sys37", "_sys38", "_sys39", "_sys40", "_sys41", "_sys42", "_sys43", "_sys44", "_sys45", "_sys46", "_sys47", "_sys48", "_sys49", "_sys50"};

const char* period = "18qr";
//const char* cutSets[] = {"", "_sys1", "_sys2", "_sys3", "_sys4", "_sys5", "_sys6", "_sys7", "_sys8", "_sys9", "_sys10", "_sys11", "_sys12", "_sys13", "_sys14", "_sys15", "_sys16", "_sys17", "_sys18", "_sys19", "_sys20", "_sys21", "_sys22", "_sys23", "_sys24", "_sys25", "_sys26", "_sys27", "_sys28", "_sys29", "_sys30"};


namespace constants{
constexpr double epsilon = 1e-3;
}

namespace fileIO{
const std::string inFileName = "data/AnalysisResultsMultitrial_LHC18qr.root"; //"./AnalysisResultsLHC18qrSysHyperloop2.root"; // AnalysisResults.root
const std::string inFileNameMC0 = "data/AnalysisResultsMultitrial_MC_LHC20g7.root"; //"./AnalysisResultsLHC20g7abcSysHyperloop.root"; // AnalysisResults.root
const std::string inFileNameMC1 = "data/AnalysisResultsMultitrial_MC_LHC20e3abc.root"; //"./AnalysisResultsLHC20e3abcSysHyperloop.root"; // AnalysisResults.root
const std::string outFileName = "out_run2_Multitrial.root"; // "out_run2_LHC18qr_grid_sys_hyperloop3.root";
// const std::string inFileName = "./AnalysisResults_LHC18qr_sys_multitrial.root"; // AnalysisResults.root
// const std::string inFileNameMC0 = "./AnalysisResults_LHC20g7abc_sysTot.root"; // AnalysisResults.root
// const std::string inFileNameMC1 = "./AnalysisResultsLHC20e3abcSysTot.root"; // AnalysisResults.root
// const std::string outFileName = "out_run2_LHC18qr_multitrial_hyperloop_test.root";
}

constexpr double calib_param[]{1.2105821371078491, 0.19601094722747803, 0.9940611124038696, 2.7593915462493896, 0.0774163082242012, 1.0421967506408691};

namespace bins{
// small
constexpr int bins_small[3]{10, 100, 1};
constexpr double xmin_small[3]{0., 0., 0.};
constexpr double xmax_small[3]{10., 100., 0.8};

constexpr int bins[3]{10, 8, 1};
constexpr double xmin[3]{0., 0., 0.};
constexpr double xmax[3]{10., 80., 0.8};
}

const char* kDataDir = "data";
const char* kCalibDir = "calib";

#endif // CONSTANTS_H
