#ifndef SCOREPARAMS_H
#define SCOREPARAMS_H

struct ScoreParams {
    ScoreParams()
        : small_interface_ratio(0.4), max_penetration(-2.5), ns_thr(0.6), rec_as_thr(0.0), lig_as_thr(0.0),
          patch_res_num(500), weights(5, 0) {}
    // bool add(const char* str);
    float small_interface_ratio;
    float max_penetration;
    float ns_thr;
    float rec_as_thr;
    float lig_as_thr;
    int patch_res_num;
    std::vector<int> weights; //[1.0-up],[-1.0,1.0], [-2.2,-1.0], [-3.6,-2.2], [-5.0,-3.6]
};

#endif // SCOREPARAMS_H
