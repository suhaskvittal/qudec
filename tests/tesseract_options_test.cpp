// Author: OpenAI Codex. CLI parsing, default preservation and decoder handoff.
#include "../main/tesseract_options.h"
#include "decoder/tesseract.h"
#include <argparse/argparse.h>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
void require(bool value,const char* message)
{if (!value) throw std::runtime_error(message);}
TesseractOptions parse(std::vector<std::string> words)
{
    TesseractOptions options;
    ARGPARSE parser;
    options.add_arguments(parser);
    std::vector<char*> argv;
    for (auto& word:words) argv.push_back(word.data());
    parser.parse(argv.size(),argv.data());
    options.validate("tesseract",argv.size(),argv.data());
    return options;
}
}
int main()
{
    const stim::DetectorErrorModel dem(R"(
        error(0.1) D0
        error(0.1) D0 D1
        error(0.1) D1 L0
        detector(0, 0, 0) D0
        detector(1, 0, 0) D1
    )");
    const auto defaults=parse({"qudec"}).make_config(dem);
    require(defaults.det_beam==5 && defaults.pqlimit==200000 && !defaults.beam_climbing &&
        defaults.no_revisit_dets && defaults.merge_errors && !defaults.verbose && defaults.det_penalty==0 &&
        !defaults.sparsify_errors && defaults.sparsify_base_degree==-1 && defaults.sparsify_max_degree==-1 &&
        defaults.sparsify_reactivate_limit==-1 && defaults.detector_orders.size()==1 &&
        defaults.detector_orders[0].get_method()==tesseract_decoder::DetOrder::Index,"Changed Tesseract defaults");
    const auto custom=parse({"qudec","--tesseract-beam","20","--tesseract-pqlimit","1000000",
        "--tesseract-beam-climbing","--tesseract-num-det-orders","21","--tesseract-det-order-seed","518278944",
        "--tesseract-det-order-bfs","--tesseract-det-order-index","--tesseract-det-order-coordinate",
        "--tesseract-allow-revisit-dets","--tesseract-no-merge-errors","--tesseract-det-penalty","0.25",
        "--tesseract-sparsify-errors","--tesseract-sparsify-base-degree","1",
        "--tesseract-sparsify-max-degree","2","--tesseract-sparsify-reactivate-limit","3"});
    auto config=custom.make_config(dem);
    require(config.dem==dem && config.det_beam==20 && config.pqlimit==1000000 && config.beam_climbing &&
        !config.no_revisit_dets && !config.merge_errors && config.det_penalty==.25 &&
        config.sparsify_errors && config.sparsify_base_degree==1 && config.sparsify_max_degree==2 &&
        config.sparsify_reactivate_limit==3,"CLI settings were not passed to TesseractConfig");
    require(config.detector_orders.size()==21,"Mixed order count is not the requested total");
    const std::vector<tesseract_decoder::DetOrder> methods{
        tesseract_decoder::DetOrder::BFS,tesseract_decoder::DetOrder::Index,tesseract_decoder::DetOrder::Coordinate};
    auto repeat=custom.make_config(dem);
    for (size_t j=0;j<config.detector_orders.size();j++)
    {
        require(config.detector_orders[j].get_method()==methods[j%3],"Wrong detector-order mixture");
        config.detector_orders[j].resolve(dem);
        tesseract_decoder::DetectorOrder expected(methods[j%3],518278944+j);
        expected.resolve(dem);
        repeat.detector_orders[j].resolve(dem);
        require(config.detector_orders[j].get_order()==expected.get_order() &&
            config.detector_orders[j].get_order()==repeat.detector_orders[j].get_order(),"Detector-order seeds are not reproducible");
    }
    std::ostringstream description;
    custom.print(description);
    require(description.str().find("beam=20 beam_climbing=1 pqlimit=1000000")!=std::string::npos &&
        description.str().find("order_methods=bfs,index,coordinate")!=std::string::npos,"Missing effective configuration description");
    dec::Tesseract decoder(config),original(dem);
    dec::SyndromeType syndrome(2);syndrome[0]=1;syndrome[1]=1;
    dec::ObsType truth(1);
    require(decoder.decode(syndrome,truth).flipped_obs==original.decode(syndrome,truth).flipped_obs,
        "Configured Tesseract could not decode the graph");
    const auto unlimited=parse({"qudec","--tesseract-beam","-1","--tesseract-pqlimit","-1"}).make_config(dem);
    require(unlimited.det_beam==tesseract_decoder::INF_DET_BEAM &&
        unlimited.pqlimit==std::numeric_limits<size_t>::max(),"Unlimited sentinels mapped incorrectly");
    for (const auto& args:std::vector<std::vector<std::string>>{
        {"qudec","--tesseract-beam","-2"},{"qudec","--tesseract-beam","65536"},
        {"qudec","--tesseract-pqlimit","0"},{"qudec","--tesseract-num-det-orders","0"},
        {"qudec","--tesseract-det-order-seed","-1"},{"qudec","--tesseract-det-penalty","nan"},
        {"qudec","--tesseract-sparsify-errors"},{"qudec","--tesseract-sparsify-base-degree","1"},
        {"qudec","--tesseract-det-order-bfs","--tesseract-det-order-coordinate"}})
    {
        bool rejected=false;
        try {(void)parse(args);} catch (const std::invalid_argument&) {rejected=true;}
        require(rejected,"Invalid Tesseract CLI configuration accepted");
    }
    char program[]="qudec",flag[]="--tesseract-beam",value[]="5";
    char* argv[]{program,flag,value};
    bool rejected=false;
    try {TesseractOptions{}.validate("pymatching",3,argv);} catch (const std::invalid_argument&) {rejected=true;}
    require(rejected,"Explicit default Tesseract option was silently ignored for another decoder");
    std::cout<<"Tesseract options: defaults, parsing, mixed order seeds, decoder handoff and invalid inputs passed\n";
}
