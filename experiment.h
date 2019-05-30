/*
GNU GPL v2
Copyright (c) 2019 Hiroki Takizawa
*/

#ifndef RINTDWR_EXPERIMENT_H_
#define RINTDWR_EXPERIMENT_H_

namespace rintdwr {

void ComputationTimeExperiment1();
void ComputationTimeExperiment2();
void ComputationTimeExperiment3();
void ComputationTimeExperiment4();
void ComputationTimeExperiment5(const int W);

void AccuracyExperiment1(const int num, const bool perform_non_fourier);
void AccuracyExperiment2();

void SuppMaxDistExperiment1();
void SuppMaxDistExperiment2();
void SuppMaxDistExperiment3();

int ComputeBestMaxSpan(
	const std::string& sequence,
	const std::string param_file_name,
	const double temperature,
	const int max_loop,
	const int sample_amount,
	const double threshold_parameter);

void HeatResistanceExperiment(const std::string filename, double threshold);

void HeatResistanceExperimentPK(const std::string sequencefilename, const std::string structurefilename, double threshold);

void RrnaAccuracyExperiment(const std::string filename);
}



#endif//RINTDWR_EXPERIMENT_H_