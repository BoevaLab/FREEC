
#include "bayes_opt.h"
#include "limbo/opt/nlopt_base.hpp"
#include <fstream>

#ifndef USE_NLOPT
#warning No NLOpt
#else
#include <limbo/opt/nlopt_base.hpp>
#endif

using namespace limbo;

struct Params
{
    struct bayes_opt_boptimizer : public defaults::bayes_opt_boptimizer
    {
        BO_PARAM(int, hp_period, -1);
    };

    struct bayes_opt_bobase : public defaults::bayes_opt_bobase
    {
        BO_PARAM(int, stats_enabled, true);
        BO_PARAM(bool, bounded, true);
        BO_DYN_PARAM(std::string, base_dir);
        BO_DYN_PARAM(std::string, dir_with_results_prefix);
    };

    struct stop_maxiterations
    {
        BO_DYN_PARAM(int, iterations);
    };

    struct acqui_ei : public defaults::acqui_ei
    {
        BO_PARAM(double, jitter, 0.01);
    };

    struct init_randomsampling
    {
        BO_PARAM(int, samples, 15);
    };

    struct kernel : public defaults::kernel
    {
        BO_DYN_PARAM(double, noise);
    };

    struct kernel_maternfivehalves : public defaults::kernel_maternfivehalves
    {
    };

    struct opt_rprop : public defaults::opt_rprop
    {
    };
    struct opt_nloptnograd : public defaults::opt_nloptnograd
    {
    };
};

double optim_scoring_func(int tau, double alpha, GenomeCopyNumber genomeCopyNumber)
{
    auto start = std::chrono::high_resolution_clock::now();
    double score = genomeCopyNumber.processGMM("/tmp/res-" + std::to_string(tau) + "-" + std::to_string(alpha), tau, alpha, true, genomeCopyNumber.getDataSubsamplingRateInPuritySearch(), false);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);

    if (score > 1e8)
    {
        std::cout << "Clipping score to 1e8 from: " << score << std::endl;
        score = 1e8;
    }
    else if (score < -1e8 || std::isnan(score))
    {
        std::cout << "Clipping score to -1e8 from: " << score << std::endl;
        score = -1e8;
    }

    std::cout << "processGMM with tau/alpha:\t" << tau << " / " << alpha << "\ttook " << duration.count() << " s" << std::endl;
    std::cout << "SCORE with tau/alpha:\t" << tau << " / " << alpha << "\t-\t" << score << std::endl;
    return score;
}

template <typename Params>
struct eval_func
{
    BO_PARAM(size_t, dim_in, 1);

    BO_PARAM(size_t, dim_out, 1);
    int tau;
    GenomeCopyNumber genomeCopyNumber;

    explicit eval_func(int tau, GenomeCopyNumber genomeCopyNumber)
    {
        this->tau = tau;
        this->genomeCopyNumber = genomeCopyNumber;
    }

    Eigen::VectorXd operator()(const Eigen::VectorXd &x) const
    {
        double alpha = x(0);
        double score = optim_scoring_func(tau, alpha, genomeCopyNumber);
        return tools::make_vector(score);
    }
};

BO_DECLARE_DYN_PARAM(int, Params::stop_maxiterations, iterations);
BO_DECLARE_DYN_PARAM(double, Params::kernel, noise);
BO_DECLARE_DYN_PARAM(std::string, Params::bayes_opt_bobase, base_dir);
BO_DECLARE_DYN_PARAM(std::string, Params::bayes_opt_bobase, dir_with_results_prefix);

void misc_eval(GenomeCopyNumber genomeCopyNumber, int tau, double alpha)
{
    double score;
    score = genomeCopyNumber.processGMM("/tmp/res-trunc-" + std::to_string(tau) + "-" + std::to_string(alpha), tau, alpha, true, genomeCopyNumber.getDataSubsamplingRateInPuritySearch(), false);
#pragma omp critical
    std::cout << "TEST SCORE (trunc) with tau/alpha:\t" << tau << " / " << alpha << "\t-\t" << score << std::endl;
    score = genomeCopyNumber.processGMM("/tmp/res-full-" + std::to_string(tau) + "-" + std::to_string(alpha), tau, alpha, true, genomeCopyNumber.getDataSubsamplingRateInPloidyEvaluation(), true);
#pragma omp critical
    std::cout << "TEST SCORE (full) with tau/alpha:\t" << tau << " / " << alpha << "\t-\t" << score << std::endl;
}

std::map<int, std::pair<double, double>> findBestPurity(GenomeCopyNumber &sampleCopyNumber, std::vector<int> ploidies, ConfigFile &cf, std::string outputDir)
{
    int optimObjectiveCalls = (int)cf.Value("bayesopt", "optimObjectiveCalls", 7);
    double kernelNoise = (double)cf.Value("bayesopt", "kernelNoise", 1e-10);
    Params::stop_maxiterations::set_iterations(optimObjectiveCalls);
    Params::kernel::set_noise(kernelNoise);
    Params::bayes_opt_bobase::set_base_dir(outputDir);

    using kernel_t = kernel::MaternFiveHalves<Params>;
    using gp_t = model::GP<Params, kernel_t>;
    using acqui_t = acqui::EI<Params, gp_t>;
    using acqui_opt_t = opt::NLOptNoGrad<Params>;
    using init_t = init::RandomSampling<Params>;

    using stat_t = boost::fusion::vector<limbo::stat::ConsoleSummary<Params>,
                                         limbo::stat::AggregatedObservations<Params>,
                                         limbo::stat::BestAggregatedObservations<Params>,
                                         limbo::stat::BestObservations<Params>,
                                         limbo::stat::Samples<Params>,
                                         limbo::stat::BestSamples<Params>,
                                         limbo::stat::GPLikelihood<Params>,
                                         limbo::stat::GPMeanHParams<Params>,
                                         limbo::stat::GPPredictionDifferences<Params>,
                                         limbo::stat::Observations<Params>,
                                         limbo::stat::GPKernelHParams<Params>>;

    std::map<int, std::pair<double, double>> best_score_alpha_per_tau;
    std::map<int, std::pair<double, double>> best_result;

    omp_lock_t lck;
    omp_init_lock(&lck);

#pragma omp parallel num_threads(ploidies.size())
    {
#pragma omp for
        for (int tau : ploidies)
        {
            std::string dir_with_results_prefix = "ploidy_" + std::to_string(tau) + "_optim_results_";
            omp_set_lock(&lck);
            Params::bayes_opt_bobase::set_dir_with_results_prefix(dir_with_results_prefix);
            bayes_opt::BOptimizer<Params, modelfun<gp_t>, acquifun<acqui_t>,
                                acquiopt<acqui_opt_t>, initfun<init_t>, statsfun<stat_t>>
                boptimizer;
            omp_unset_lock(&lck);
            boptimizer.optimize(eval_func<Params>(tau, sampleCopyNumber));

            std::string filename = boptimizer.res_dir() + "/tau-" + std::to_string(tau) + "-alpha-" + std::to_string(boptimizer.best_sample()(0));
            std::ofstream(filename.c_str());

#pragma omp critical
            best_score_alpha_per_tau.insert({tau, std::make_pair(boptimizer.best_observation()(0), boptimizer.best_sample()(0))});
        }
    }
    omp_destroy_lock(&lck);
    
    float best_purity = -INFINITY;
    float best_score = -INFINITY;
    float best_tau = -1;
    std::string clusteringOutputDir = outputDir+"/clustering_results";
    boost::filesystem::create_directories(clusteringOutputDir);
#pragma omp parallel num_threads(ploidies.size())
    {
#pragma omp for
    for (int i = 0; i < best_score_alpha_per_tau.size(); i++)
    {
        auto x = best_score_alpha_per_tau.begin();
        std::advance(x, i);
        double score_on_trunc = x->second.first;
        double purity = x->second.second;
        int tau = x->first;

        auto start = std::chrono::high_resolution_clock::now();
        double score_on_full = sampleCopyNumber.processGMM(clusteringOutputDir + "/res-" + std::to_string(tau) + "-" + std::to_string(purity), tau, purity, true, sampleCopyNumber.getDataSubsamplingRateInPloidyEvaluation(), true);
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);

        #pragma omp critical
        {
            std::cout << "processGMM (FULL) with tau/alpha:\t" << tau << " / " << purity << "\ttook " << duration.count() << " s" << std::endl;
            std::cout << "TAU: " << tau << "\nBest sample: " << purity << " - Best observation (trunc): " << score_on_trunc << " - Best observation (full): " << score_on_full << std::endl;
            if (score_on_full > best_score)
            {
                best_purity = purity;
                best_score = score_on_full;
                best_tau = tau;
            }
        }
    }
    }

    best_result.insert({best_tau, std::make_pair(best_score, best_purity)});
    return best_result;
}
