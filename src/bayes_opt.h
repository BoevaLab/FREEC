#ifndef _BAYES_OPT_H
#define _BAYES_OPT_H

#include <iostream>
#include <omp.h>

#include <limbo/bayes_opt/boptimizer.hpp>
#include "limbo/acqui/ei.hpp"
#include "limbo/opt/parallel_repeater.hpp"
#include "limbo/stat.hpp"

#include "GenomeGMM.h"
#include "GenomeCopyNumber.h"

struct Params;


double concrete_eval_func(double alpha, const Eigen::MatrixXd &data);


double my_eval_func(int tau, double alpha, const std::vector<float>& ratio,
                  const std::vector<float>& maf, GenomeGMM sampleGMM);


template<typename Params>
struct eval_func;


template<typename Params>
struct MeanFWModel;


std::map<int, std::pair<double, double>> findBestPurity(GenomeCopyNumber &sampleCopyNumber, std::vector<int> ploidies, ConfigFile& cf, std::string outputDir);

#endif