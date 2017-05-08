#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <list>
#include <numeric>
#include <random>
#include <functional>
#include <utility>
#include <iomanip>
#include <set>
#include <cstdint>
#include <assert.h>
#include <unistd.h>
#include <boost/math/distributions/chi_squared.hpp>

#define Int int32_t
Int IDX_FEATURE;
double NUM_PATTERN;
bool VERBOSE = false;
bool WY = false;
double SIGMA = 0.2;
double ADMISSIBLE;
double FWER_ESTIMATE;

using namespace std;

struct cmp {
	bool operator ()(const pair<double, double> &a,
			const pair<double, double> &b) {
		// return a.second < b.second;
		return a.first > b.first;
	}
};

// output a 2D vector
template<typename T>
ostream &operator<<(ostream& out, const vector<vector<T>>& mat) {
	for (Int i = 0; i < mat.size() - 1; i++) {
		for (auto&& x : mat[i]) {
			out << x << " ";
		}
		out << endl;
	}
	for (auto&& x : mat.back()) {
		out << x << " ";
	}
	return out;
}
// output a vector
template<typename T>
ostream &operator<<(ostream& out, const vector<T>& vec) {
	if (vec.size() == 0)
		return out;
	for (Int i = 0; i < vec.size() - 1; ++i) {
		out << vec[i] << " ";
	}
	out << vec.back();
	return out;
}
// output a set
template<typename T>
ostream &operator<<(ostream& out, const set<T>& vec) {
	auto first = vec.begin();
	auto last = vec.empty() ? vec.end() : prev(vec.end()); // in case vec is empty
	for (auto it = first; it != last; ++it) {
		out << *it << ", ";
	}
	out << *(prev(vec.end()));
	return out;
}

// read a database file
void readFromCSV(ifstream& ifs, vector<vector<double>>& data) {
	string line;
	vector<vector<double>> tdata;
	while (getline(ifs, line)) {
		stringstream lineStream(line);
		string cell;
		vector<double> tmp;
		while (getline(lineStream, cell, ',')) {
			tmp.push_back(stod(cell));
		}
		tdata.push_back(tmp);
	}

	// TODO: cache efficient transposition
	int nTrans = tdata.size();
	int nFeatures = tdata[0].size();
	data.resize(nFeatures);
	for (int f = 0; f < nFeatures; ++f) {
		data[f].resize(nTrans);
	}
	for (int t = 0; t < nTrans; ++t) {
		for (int f = 0; f < nFeatures; ++f) {
			data[f][t] = tdata[t][f];
		}
	}
}

// read a database file
void readFromCSV(ifstream& ifs, vector<vector<double>>& data,
		int dim_limit,
		bool reverse) {
	string line;
	while (getline(ifs, line)) {
		stringstream lineStream(line);
		string cell;
		vector<double> tmp;
		while (getline(lineStream, cell, ',')
				&& tmp.size() <= dim_limit) {
			tmp.push_back(stod(cell));
		}
		// YJ: Add reversed (negated) values for all features so that the method can detect negative correlation.
		if (reverse) {
			int size = tmp.size();
			for (int i = 0; i < size; ++i) {
				tmp.push_back(-tmp[i]);
			}
		}
		data.push_back(tmp);
	}
}

void readClassFromCSV(ifstream& ifs, vector<Int>& cl) {
	string line;
	while (getline(ifs, line)) {
		cl.push_back(stoi(line));
	}
}

// compute frequency (not used)
double computeFreq(vector<vector<Int>>& rank, vector<Int>& fset,
		vector<double>& freq_current) {
	Int N = rank.size();
	double freq = 0.0;
	if (fset.size() == 0)
		return 1.0;
	for (Int i = 0; i < N; ++i) {
		double freq_each = 1.0;
		for (auto&& j : fset) {
			freq_each *= (double) rank[i][j] / (double) N;
		}
		freq_current[i] = freq_each;
		freq += freq_each;
	}
	freq /= (double) N;
	return freq;
}
// compute frequency using the previous result
double computeFreqUpdate(vector<vector<double>>& rankn, Int fnew,
		vector<double>& freq_current) {
	Int N = rankn[0].size();
	double freq = 0.0;
	// if (fset.size() == 0) return 1.0;
	for (Int i = 0; i < N; ++i) {
		freq_current[i] *= (double) rankn[fnew][i];
		freq += freq_current[i];
	}
	freq /= (double) N;
	return freq;
}

// compute the maximum frequency (eta)
double eta_max(Int k, Int N, Int N_class) {
	double sum = 0.0;
	for (Int i = 1; i <= N_class; ++i) {
		sum += pow((double) (N - i + 1) / (double) N, (double) k);
	}
	sum /= (double) N;
	return sum;
}
// compute the mininum frequency (eta)
double eta_min(Int k, Int N, Int N_class) {
	vector<Int> idx(N_class);
	iota(idx.begin(), idx.end(), 1);
	// ascending order vector
	vector<double> ascending_order(N_class);
	for (Int i = 0; i < N_class; ++i) {
		ascending_order[i] = (double) idx[i] / (double) N;
	}
	// current value
	vector<double> v;
	copy(ascending_order.begin(), ascending_order.end(),
			back_inserter(v));

	for (Int j = 0; j < k - 1; ++j) {
		sort(v.begin(), v.end(), greater<double>());
		for (Int i = 0; i < N_class; ++i) {
			v[i] *= ascending_order[i];
		}
	}
	double sum = accumulate(v.begin(), v.end(), 0.0);
	sum /= (double) N;
	return sum;
}
// compute the maximum achievable frequency for the smaller class
double eta_max_0(vector<double>& freq_current, Int N0, Int N) {
	vector<double> v;
	copy(freq_current.begin(), freq_current.end(), back_inserter(v));
	partial_sort(v.begin(), v.begin() + N0, v.end(),
			greater<double>());
	double freq0 = accumulate(v.begin(), v.begin() + N0, 0.0);
	freq0 /= (double) N;
	return freq0;
}
// compute the maximum achievable KL divergence
double kl_max(vector<double>& freq_current, Int N0, Int N) {
	vector<double> po;
	vector<double> pe;
	double r0 = (double) N0 / (double) N;
	double r1 = (double) (N - N0) / (double) N;
	double freq0 = eta_max_0(freq_current, N0, N);
	double freq = accumulate(freq_current.begin(), freq_current.end(),
			0.0);
	freq /= (double) N;

	// cout << endl << "freq =  " << freq << endl;
	// cout << "freq0 = " << freq0 << endl;

	po.push_back(freq0);
	po.push_back(freq - freq0);
	po.push_back(r0 - freq0);
	po.push_back(r1 - freq + freq0);

	pe.push_back(r0 * freq);
	pe.push_back(r1 * freq);
	pe.push_back(r0 - r0 * freq);
	pe.push_back(r1 - r1 * freq);

	// cout << "po: " << po << endl;
	// cout << "pe: " << pe << endl;

	double kl = 0.0;
	for (Int i = 0; i < po.size(); ++i) {
		kl += po[i] * log(po[i] / pe[i]);
	}
	return kl;
}

double kl_max_fast_upper(double freq, Int N0, Int N) {
	assert(0 <= freq && freq <= 1);
	assert(0 <= N0 && N0 <= N);
	double r0 = (double) N0 / (double) N;
	if (freq < r0) {
		double r = freq * log(1.0 / r0)
				+ (r0 - freq) * log((r0 - freq) / (r0 - r0 * freq))
				+ (1.0 - r0) * log(1.0 / (1.0 - freq));
		assert(0.0 <= r);
		return r;
	} else {
		double r = r0 * log(1 / r0) + (1 - r0) * log(1 / (1 - r0));
		assert(0.0 <= r);
		return r;
	}
}
double kl_max_fast(double freq, Int N0, Int N) {
	assert(0 <= freq && freq <= 1);
	assert(0 <= N0 && N0 <= N);
	double r0 = (double) N0 / (double) N;
	if (freq < r0) {
		double r = freq * log(1.0 / r0)
				+ (r0 - freq) * log((r0 - freq) / (r0 - r0 * freq))
				+ (1.0 - r0) * log(1.0 / (1.0 - freq));
		assert(0.0 <= r);
		return r;
	} else {
		double r = r0 * log(1 / freq)
				+ (freq - r0) * log((freq - r0) / (freq - freq * r0))
				+ (1 - freq) * log(1 / (1 - r0));
		assert(0.0 <= r);
		return r;
	}
}

double kl(vector<double>& freq_current, vector<Int>& cl, double freq,
Int N0,
Int N) {
	double r0 = (double) N0 / (double) N;
	double r1 = (double) (N - N0) / (double) N;
	double freq0 = 0.0;
	for (Int i = 0; i < N; ++i) {
		if (cl[i] == 0)
			freq0 += freq_current[i];
	}
	freq0 /= (double) N;
	double freq1 = freq - freq0;

	vector<double> po;
	po.push_back(freq0);
	po.push_back(freq1);
	po.push_back(r0 - freq0);
	po.push_back(r1 - freq1);

	vector<double> pe;
	pe.push_back(r0 * freq);
	pe.push_back(r1 * freq);
	pe.push_back(r0 - r0 * freq);
	pe.push_back(r1 - r1 * freq);

	double kl = 0.0;
	for (Int i = 0; i < po.size(); ++i) {
		kl += po[i] * log(po[i] / pe[i]);
	}
	return kl;
}
// compute p-value
double computePvalue(double kl, Int N) {
	boost::math::chi_squared chisq_dist(1);
	// else pval = 1 - boost::math::cdf(chisq_dist, 2 * (double)N * kl);
	// if (pval > 1) pval = 1.0;
	// if (VERBOSE) cout << "kl: " << kl << endl;
	double pval = 0.0;
	if (kl <= pow(10, -8))
		pval = 1.0;
	else
		pval = 1.0
				- boost::math::cdf(chisq_dist, 2.0 * (double) N * kl);
	return pval;
}
// compute a list of thresholds
void computeThresholds(vector<double>& freq_thrs, Int n, Int N) {
	freq_thrs.push_back(eta_max(1, N, N));
	for (Int k = n; k >= 2; --k) {
		freq_thrs.push_back(eta_max(k, N, N));
		freq_thrs.push_back(eta_min(k, N, N));
	}
	sort(freq_thrs.begin(), freq_thrs.end(), greater<double>());
}

// Frequent Pattern Mining
void runFPM(vector<vector<double>>& rankn, vector<Int>& fset,
		vector<double>& freq_current, Int i_prev, Int n,
		Int size_limit, ofstream& ofs) {
	Int N = rankn[0].size();
	for (Int i = i_prev + 1; i < n; i++) {
		fset.push_back(i);
		double freq = computeFreqUpdate(rankn, i, freq_current);
		if (freq > SIGMA && fset.size() <= size_limit) {
			if (VERBOSE)
				cout << "  Frequency of {" << fset << "} = " << freq
						<< endl;
			// ofs << fset << " (" << freq << ")" << endl;
			NUM_PATTERN += 1.0;
			runFPM(rankn, fset, freq_current, i, n, size_limit, ofs);
		}
		// extract the feature fset.back();
		for (Int i = 0; i < N; ++i) {
			freq_current[i] /= (double) rankn[i][fset.back()];
		}
		fset.pop_back();
	}
}

/**
 * Significant Pattern Mining for continuous variables
 * Depth-first search with recursive call.
 *
 * Arguments
 *	rankn			: Nxn dataset.
 *	fset           	: current pattern
 *	freq_current	: vector of frequencies for transactions
 *	i_prev			: iterator for branch (dfs)
 *	n				: number of types (items)
 *	size_limit		: INT_MAX is put -- I guess it is for debugging.
 *	ofs				: ???
 *	N0 				: number of transactions with class negative (fixed)
 *	freq_pval_list 	: list of testable pattern currently hold. pair<frequency, p-value>
 *	alpha			: error rate (fixed)
 */
void runSPM(vector<vector<double>>& rankn, vector<Int>& fset,
		vector<double>& freq_current, Int i_prev, Int n,
		Int size_limit, ofstream& ofs, Int N0,
		set<pair<double, double>, cmp>& freq_pval_list,
		double alpha) {
	Int N = rankn[0].size();
	for (Int i = i_prev + 1; i < n; i++) {
		fset.push_back(i);
		double freq = computeFreqUpdate(rankn, i, freq_current);
		if (freq > SIGMA && fset.size() <= size_limit) {
			NUM_PATTERN += 1.0;
			if (VERBOSE)
				cout << "  Frequency of {" << fset << "} = " << freq
						<< endl;
			// if (freq < (double)N0 / (double)N) {
			// double pval = computePvalue(kl_max(freq_current, N0, N), N);
			double pval = computePvalue(kl_max_fast(freq, N0, N), N);
			freq_pval_list.insert(make_pair(freq, pval));
			ADMISSIBLE = alpha / (*prev(freq_pval_list.end())).second;
			while (ADMISSIBLE < NUM_PATTERN) {
				// SIGMA = (*prev(freq_pval_list.end())).first;
				SIGMA = min((double) N0 / (double) N,
						(*prev(freq_pval_list.end())).first);
				// cout << "update sigma:" << SIGMA << endl;
				freq_pval_list.erase(prev(freq_pval_list.end()));
				NUM_PATTERN -= 1.0;
				ADMISSIBLE =
						freq_pval_list.empty() ?
								1e20 :
								alpha
										/ (*prev(freq_pval_list.end())).second;
			}
			// }
			// ofs << fset << " (" << freq << ")" << endl;
			if (VERBOSE)
				cout << "  Current freq threshold =    " << SIGMA
						<< endl;
			if (VERBOSE)
				cout << "  Current admissible number = " << ADMISSIBLE
						<< endl;
			runSPM(rankn, fset, freq_current, i, n, size_limit, ofs,
					N0, freq_pval_list, alpha);
		}
		// extract the feature fset.back();
		// TODO: this is the part it takes the most time.
		// ok DFS is better than BrFS i guess.
		for (Int i = 0; i < N; ++i) {
			freq_current[i] /= (double) rankn[fset.back()][i];
		}
		fset.pop_back();
	}
}

double GetLowerBoundOfPValue(double freq, int N0, int N) {
	double upperBound = kl_max_fast_upper(freq, N0, N);
	return computePvalue(upperBound, N);
}

/**
 * Initialize a table to store the count of each frequency.
 * N: number of transactions.
 * size: size of the table (thus deciding the space efficiency).
 * @return: thresholds for discretized frequency and
 *  		the minimal pvalue for each threshold.
 */
vector<pair<double, double> > InitializeThresholdTable(int N, int N0,
		int size, double alpha, double thre_ratio) {
	// TODO: the table should be more efficient with inversed.
	vector<double> thresholds(size);
	double max_freq = (double) (N + 1.0) / (double) (2.0 * N);
	// TODO: Current discretization is way too rough.
	//       Need to find a way to edit the granularity.
	for (int i = 0; i < size; ++i) {
		max_freq = max_freq * thre_ratio;
		thresholds[i] = max_freq;
		double pbound = GetLowerBoundOfPValue(thresholds[i], N0, N);
		if (pbound >= alpha) {
			thresholds.erase(thresholds.begin() + i,
					thresholds.end());
			break;
		}
	}
	sort(thresholds.begin(), thresholds.end());

	vector<pair<double, double> > table(thresholds.size());
	for (int i = 0; i < table.size(); ++i) {
		table[i].first = thresholds[i];
		table[i].second = GetLowerBoundOfPValue(table[i].first, N0,
				N);
	}
	printf("The domain of discrete Fr(X) = {0..%d}\n", table.size());
	for (int i = 0; i < table.size(); ++i) {
		printf("freq/minp = %.2f/%.6f\n", table[i].first,
				table[i].second);
	}
	return table;
}

/**
 * It takes the maximum frequency
 */
int GetDiscretizedFrequency(vector<pair<double, double> >& thresholds,
		double freq) {
	assert(0 <= freq);
	assert(freq <= 1.0);
	int i = 0;
	while (i < thresholds.size() && freq > thresholds[i].first) {
		++i;
	}
	if (i == thresholds.size()) {
//		printf("Fr_d(%.2f) = %d\n", freq, i);
		--i;
	} else {
		assert(thresholds[i].first > freq);
	}
	assert(0 <= i);
	assert(i < thresholds.size());
	return i;
}

/**
 * It takes the maximum frequency
 */
int GetDiscretizedFrequencyInadmissible(
		vector<pair<double, double> >& thresholds, double freq) {
	assert(0 <= freq);
	assert(freq <= 1.0);
	int i = 0;
	while (i < thresholds.size() && freq > thresholds[i].first) {
		++i;
	}
	if (i == thresholds.size()) {
//		printf("Fr_d(%.2f) = %d\n", freq, i);
		--i;
	} else {
		assert(thresholds[i].first > freq);
	}
	assert(0 <= i);
	assert(i < thresholds.size());
	return i;
}

void printCountTable(vector<int>& freq_count,
		vector<pair<double, double> >& thresholds) {
	assert(freq_count.size() == thresholds.size());
	printf("The domain of discrete Fr(X) = {0..%d}\n",
			thresholds.size() - 1);
	for (int i = 0; i < thresholds.size(); ++i) {
		if (i == thresholds.size() - 1) {
			printf(
					"[%2d] freq = %.2f, minp = %.6f, #>items = %4d (accum = %4d)\n",
					i, thresholds[i].first, thresholds[i].second,
					freq_count[i], freq_count[i]);
		} else {
			printf(
					"[%2d] freq = %.2f, minp = %.6f, #>items = %4d (accum = %4d)\n",
					i, thresholds[i].first, thresholds[i].second,
					freq_count[i], freq_count[i] - freq_count[i + 1]);
		}
	}
}

// TODO: not sure what it is doing...
void IncCsAccum(vector<int>* freq_count, int disFreq, int lambda) {
	assert(0 <= disFreq);
	assert(disFreq < freq_count->size());
	for (int i = 0; i < disFreq; i++) {
		freq_count->operator [](i)++;}
	}

// TODO: Implement a table to store the count of discretized frequency.
void runLSSPM(vector<vector<double>>& rankn, vector<Int>& fset,
		vector<double>& freq_current, Int i_prev, Int n,
		Int size_limit, ofstream& ofs, Int N0,
		vector<int>* freq_count,
		vector<pair<double, double> >& thresholds, int& curr_lambda,
		double alpha) {

	Int N = rankn[0].size();
	for (Int i = i_prev + 1; i < n; i++) {
		fset.push_back(i);
		double freq = computeFreqUpdate(rankn, i, freq_current);
		int disFreq = GetDiscretizedFrequency(thresholds, freq);
		if (disFreq > curr_lambda && fset.size() <= size_limit) {
			NUM_PATTERN += 1.0;
			if (VERBOSE)
				cout << "  Frequency of {" << fset << "} = " << freq
						<< endl;

			IncCsAccum(freq_count, disFreq, curr_lambda);

			if (freq_count->operator [](curr_lambda)
					* thresholds[curr_lambda].second >= alpha) {
				++curr_lambda;
				assert(curr_lambda < thresholds.size());
			}

			// }
			// ofs << fset << " (" << freq << ")" << endl;
			if (VERBOSE)
				cout << "  Current freq threshold =    " << SIGMA
						<< endl;
			if (VERBOSE)
				cout << "  Current admissible number = " << ADMISSIBLE
						<< endl;
			runLSSPM(rankn, fset, freq_current, i, n, size_limit, ofs,
					N0, freq_count, thresholds, curr_lambda, alpha);
		}
		// extract the feature fset.back();
		// TODO: this is the part it takes the most time.
		// ok DFS is better than BrFS i guess.
		for (Int i = 0; i < N; ++i) {
			freq_current[i] /= (double) rankn[fset.back()][i];
		}
		fset.pop_back();
	}
}

// Significant Pattern Mining with Westfall-Young permutation
void runSPMWY(vector<vector<double>>& rankn, vector<Int>& fset,
		vector<double>& freq_current, Int i_prev, Int n,
		Int size_limit, ofstream& ofs, Int N0,
		set<pair<double, double>, cmp>& freq_pval_list,
		vector<vector<Int>>& cl_perm, vector<double>& pmin_perm,
		double alpha) {
	assert(false && "rankn not calced yet");
	Int N = rankn.size();
	for (Int i = i_prev + 1; i < n; i++) {
		fset.push_back(i);
		double freq = computeFreqUpdate(rankn, i, freq_current);
		if (freq > SIGMA && fset.size() <= size_limit) {
			NUM_PATTERN += 1.0;
			if (VERBOSE)
				cout << "  Frequency of {" << fset << "} = " << freq
						<< endl;
			double pval = computePvalue(kl_max_fast(freq, N0, N), N);
			freq_pval_list.insert(make_pair(freq, pval));

			// compute p-values for permutations and estimate the FWER
			for (Int ip = 0; ip < cl_perm.size(); ++ip) {
				double p_tmp = computePvalue(
						kl(freq_current, cl_perm[ip], freq, N0, N),
						N);
				if (p_tmp < pmin_perm[ip])
					pmin_perm[ip] = p_tmp;
			}
			Int counter = 0;
			for (auto&& x : pmin_perm) {
				if (x < (*prev(freq_pval_list.end())).second)
					counter++;
			}
			FWER_ESTIMATE = (double) counter
					/ (double) cl_perm.size();
			// sort(pmin_perm.begin(), pmin_perm.end());
			// cout << (*prev(freq_pval_list.end())).second << ": ";
			// for (Int ii = 0; ii < 10; ii++) cout << pmin_perm[ii] << " ";
			// cout << "(" << FWER_ESTIMATE << ")" << endl;

			// cout << "Minimum achivable p-value: " << (*prev(freq_pval_list.end())).second << endl;
			// cout << "pmin_perm: " << pmin_perm << endl;
			// cout << "FWER estimation: " << FWER_estimate << endl;
			while (alpha < FWER_ESTIMATE) {
				SIGMA = min((double) N0 / (double) N,
						(*prev(freq_pval_list.end())).first);
				// cout << "update sigma:" << SIGMA << endl;
				freq_pval_list.erase(prev(freq_pval_list.end()));
				NUM_PATTERN -= 1.0;
				if (freq_pval_list.empty()) {
					FWER_ESTIMATE = 0.0;
				} else {
					counter = 0;
					for (auto&& x : pmin_perm) {
						if (x < (*prev(freq_pval_list.end())).second)
							counter++;
					}
					FWER_ESTIMATE = (double) counter
							/ (double) cl_perm.size();
				}
			}
			if (VERBOSE)
				cout << "  Current freq threshold =    " << SIGMA
						<< endl;
			if (VERBOSE)
				cout << "  Current admissible number = " << ADMISSIBLE
						<< endl;
			runSPMWY(rankn, fset, freq_current, i, n, size_limit, ofs,
					N0, freq_pval_list, cl_perm, pmin_perm, alpha);
		}
		// extract the feature fset.back();
		for (Int i = 0; i < N; ++i) {
			freq_current[i] /= (double) rankn[i][fset.back()];
		}
		fset.pop_back();
	}
}
// Frequent Pattern Mining with finding singificant patterns
void runSPM_sig(vector<vector<double>>& rankn, vector<Int>& fset,
		vector<double>& freq_current, Int i_prev, Int n,
		Int size_limit, Int N0, vector<Int>& cl,
		double alpha_corrected, ofstream& ofs) {
	Int N = rankn[0].size();
	for (Int i = i_prev + 1; i < n; i++) {
		fset.push_back(i);
		double freq = computeFreqUpdate(rankn, i, freq_current);
		if (freq > SIGMA && fset.size() <= size_limit) {
			// TODO: Need to calculate pmin first to assess if it is a testable pattern?
			double pval = computePvalue(
					kl(freq_current, cl, freq, N0, N), N);
			// if (VERBOSE) cout << "  corrected p-value of {" << fset << "} = " << pval * (double)num_testable << endl;
			if (pval < alpha_corrected) {
				NUM_PATTERN += 1.0;
				// ofs << fset << " (" << pval * (double)num_testable << ")" << endl;
				ofs << fset << " (" << pval << ")" << endl;
			}
			runSPM_sig(rankn, fset, freq_current, i, n, size_limit,
					N0, cl, alpha_corrected, ofs);
		}
		// extract the feature fset.back();
		for (Int i = 0; i < N; ++i) {
			freq_current[i] /= (double) rankn[fset.back()][i];
		}
		fset.pop_back();
	}
}

int main(int argc, char *argv[]) {
	bool flag_in = false;
	bool flag_class_in = false;
	bool flag_out = false;
	bool flag_stat = false;
	bool fpm = false;
	bool bonferroni = false;
	bool lscpm = false;
	double thre_ratio = 0.9;
	char *input_file = NULL;
	char *input_class_file = NULL;
	char *output_file = NULL;
	char *stat_file = NULL;
	Int size_limit = INT32_MAX;
	Int wy_repeat = 1000;
	double alpha = 0.05;
	clock_t ts, te;

	Int dim_limit = INT32_MAX;

// YJ: In order to detect negative correlation, current method needs to add revered rank.
	bool reverse = false;

// get arguments
	char opt;
	while ((opt = getopt(argc, argv, "i:c:o:t:s:a:k:fvwp:bl:d:r"))
			!= -1) {
		switch (opt) {
		case 'i':
			input_file = optarg;
			flag_in = true;
			break;
		case 'c':
			input_class_file = optarg;
			flag_class_in = true;
			break;
		case 'o':
			output_file = optarg;
			flag_out = true;
			break;
		case 't':
			stat_file = optarg;
			flag_stat = true;
			break;
		case 's':
			SIGMA = atof(optarg);
			break;
		case 'a':
			alpha = atof(optarg);
			break;
		case 'k':
			size_limit = atoi(optarg);
			break;
		case 'f':
			fpm = true;
			break;
		case 'v':
			VERBOSE = true;
			break;
		case 'w':
			WY = true;
			break;
		case 'p':
			wy_repeat = atoi(optarg);
			break;
		case 'b':
			bonferroni = true;
			break;
		case 'l':
			lscpm = true;
			thre_ratio = atof(optarg);
			break;
		case 'd':
			dim_limit = atoi(optarg);
			break;
		case 'r':
			reverse = true;
			break;
		}
	}

	if (!flag_in) {
		cerr << "> ERROR: Input file (-i [input_file]) is missing!"
				<< endl;
		exit(1);
	}
	if (!flag_out) {
		output_file = (char *) "output";
	}
	if (!flag_stat) {
		stat_file = (char *) "stat";
	}
	ofstream sfst(stat_file);

// --------------------------------- //
// ---------- Read a file ---------- //
// --------------------------------- //
	cout << "> Reading a database file \"" << input_file << "\" ... "
			<< flush;
	sfst << "> Reading a database file \"" << input_file << "\" ... ";
	ifstream ifs(input_file);
	assert(ifs.is_open());
	vector<vector<double>> data;
	vector<Int> cl;
	readFromCSV(ifs, data);
//	readFromCSV(ifs, data, dim_limit, reverse);
	cout << "end" << endl << flush;
	sfst << "end" << endl;
	Int n = data.size();
	Int N = data[0].size();
	Int N0 = 0;
	if (!fpm) {
		cout << "> Reading a class file    \"" << input_class_file
				<< "\" ... " << flush;
		sfst << "> Reading a class file    \"" << input_class_file
				<< "\" ... ";
		ifstream cfs(input_class_file);
		assert(cfs.is_open());
		readClassFromCSV(cfs, cl);
		cout << "end" << endl << flush;
		sfst << "end" << endl;
		for (auto&& x : cl)
			if (x == 0)
				N0++;
		if (N0 > N / 2) {
			cerr << "> ERROR: Class 0 must be a minor class!" << endl;
			exit(1);
		}
	}
	cout << "  # samples in total:     " << N << endl << flush;
	sfst << "  # samples in total:     " << N << endl;
	if (!fpm) {
		cout << "  # samples in class 0:   " << N0 << endl << flush;
		sfst << "  # samples in class 0:   " << N0 << endl;
	}
	cout << "  # features:             " << n << endl << flush;
	sfst << "  # features:             " << n << endl;
	if (fpm) {
		cout << "  frequency threshold:    " << SIGMA << endl
				<< flush;
		sfst << "  frequency threshold:    " << SIGMA << endl;
	}

// ----------------------------------- //
// ---------- Compute ranks ---------- //
// ----------------------------------- //
	vector<vector<Int>> rank;
	rank = vector<vector<Int>>(n, vector<Int>(N, 0));
	vector<vector<double>> rankn;
	rankn = vector<vector<double>>(n, vector<double>(N, 0.0));
	vector<size_t> idx(N);
	for (Int j = 0; j < n; ++j) {
		// initialize index vector
		iota(idx.begin(), idx.end(), 0);
		IDX_FEATURE = j;
		// sort indexes based on comparing values in v
		sort(idx.begin(), idx.end(),
				[&data](Int i1, Int i2) {return data[IDX_FEATURE][i1] < data[IDX_FEATURE][i2];});
		for (Int i = 0; i < N; i++) {
			rank[j][i] = idx[i] + 1;
			rankn[j][i] = (double) rank[j][i] / (double) N;
		}
	}

// --------------------------------- //
// ---------- Enumeration ---------- //
// --------------------------------- //
	if (VERBOSE) {
		cout << "> Start enumeration:" << endl << flush;
	} else {
		cout << "> Start enumeration ... " << flush;
		sfst << "> Start enumeration ... ";
	}
	vector<double> freq_current(N, 1.0); // track the current rank product for each sample
	vector<Int> fset; // pattern (a set of features)
	ofstream ofs(output_file);
	if (VERBOSE)
		cout << "  Frequency of (" << fset << ") = "
				<< computeFreq(rank, fset, freq_current) << endl;

	if (fpm) {
		ts = clock();
		runFPM(rankn, fset, freq_current, -1, n, size_limit, ofs);
		te = clock();
		if (!VERBOSE) {
			cout << "end" << endl << flush;
			sfst << "end" << endl << flush;
		}
		cout << "  # frequent patterns:    " << NUM_PATTERN << endl
				<< flush;
		sfst << "  # frequent patterns:    " << NUM_PATTERN << endl;
		cout << "  Running time:           "
				<< (float) (te - ts) / CLOCKS_PER_SEC << " [sec]"
				<< endl << flush;
		sfst << "  Running time:           "
				<< (float) (te - ts) / CLOCKS_PER_SEC << " [sec]"
				<< endl;
		exit(1);
	}

	set<pair<double, double>, cmp> freq_pval_list;
// freq_pval_list.insert(make_pair(0.0, 1e-20));
	SIGMA = 0.0;
	ADMISSIBLE = 1e20;
	ts = clock();
	double alpha_corrected;
	if (!WY && !bonferroni && !lscpm) {
		// =====
		// lamp
		// =====
		runSPM(rankn, fset, freq_current, -1, n, size_limit, ofs, N0,
				freq_pval_list, alpha);
		alpha_corrected = alpha / NUM_PATTERN;
	} else if (lscpm) {
		// =====
		// Linear Space Continuous Pattern Mining LAMP
		// =====
		int GRAN_DISCRETE = 1000;
		vector<pair<double, double>> thresholds =
				InitializeThresholdTable(N, N0, GRAN_DISCRETE, alpha, thre_ratio);
//		printf("threshold%d\n", thresholds.size());
		vector<int>* count = new vector<int>(thresholds.size(), 0);
		int lambda = 0;
		runLSSPM(rankn, fset, freq_current, -1, n, size_limit, ofs,
				N0, count, thresholds, lambda, alpha);
//		++lambda;
		int num_pat = count->operator [](lambda);
//		for (int i = 0; i < lambda; ++i) {
//			num_pat += count->operator [](i);
//		}
		NUM_PATTERN = num_pat;
		alpha_corrected = alpha / NUM_PATTERN;
		SIGMA = thresholds[lambda - 1].first;
		printCountTable(*count, thresholds);
		delete count;
	} else if (bonferroni) {
		// =====
		// Bonferroni
		// =====
		SIGMA = 0.0;
		runFPM(rankn, fset, freq_current, -1, n, size_limit, ofs);
		NUM_PATTERN += 1.0;
		alpha_corrected = alpha / NUM_PATTERN;
	} else {
		// =====
		// Westfall-Young permutation
		// =====
		// prepare random permutations
		vector<vector<Int>> cl_perm(wy_repeat);
		unsigned seed = 0;
		for (auto&& cl_each : cl_perm) {
			copy(cl.begin(), cl.end(), back_inserter(cl_each));
			// random_shuffle(cl_each.begin(), cl_each.end());
			shuffle(cl_each.begin(), cl_each.end(),
					default_random_engine(seed));
			seed++;
		}
		vector<double> pmin_perm(cl_perm.size(), 1.0);
		runSPMWY(rankn, fset, freq_current, -1, n, size_limit, ofs,
				N0, freq_pval_list, cl_perm, pmin_perm, alpha);
		sort(pmin_perm.begin(), pmin_perm.end());
		// cout << pmin_perm << endl;
		alpha_corrected = pmin_perm[wy_repeat * alpha];
	}
// Int num_testable = NUM_PATTERN;
	te = clock();
	if (!VERBOSE) {
		cout << "end" << endl << flush;
		sfst << "end" << endl << flush;
	}
	cout << "  # testable patterns:    " << NUM_PATTERN << endl
			<< flush;
	sfst << "  # testable patterns:    " << NUM_PATTERN << endl;
	cout << "  Corrected alpha:        " << alpha_corrected << endl
			<< flush;
	sfst << "  Corrected alpha:        " << alpha_corrected << endl;
	cout << "  Frequency threshold:    " << SIGMA << endl << flush;
	sfst << "  Frequency threshold:    " << SIGMA << endl;
	cout << "  Running time:           "
			<< (float) (te - ts) / CLOCKS_PER_SEC << " [sec]" << endl
			<< flush;
	sfst << "  Running time:           "
			<< (float) (te - ts) / CLOCKS_PER_SEC << " [sec]" << endl;

// ------------------------------------------------------------------------- //
// ---------- Enumerate significant patterns with corrected alpha ---------- //
// ------------------------------------------------------------------------- //
	cout << "> Run Frequent Pattern Mining with a threshold " << SIGMA
			<< endl << flush;
	sfst << "> Run Frequent Pattern Mining with a threshold " << SIGMA
			<< endl;
	NUM_PATTERN = 0.0;
	ts = clock();
	for (auto&& x : freq_current)
		x = 1.0;
	runSPM_sig(rankn, fset, freq_current, -1, n, size_limit, N0, cl,
			alpha_corrected, ofs);
	te = clock();
	cout << "  # significant patterns: " << NUM_PATTERN << endl
			<< flush;
	sfst << "  # significant patterns: " << NUM_PATTERN << endl;
	cout << "  Running time:           "
			<< (float) (te - ts) / CLOCKS_PER_SEC << " [sec]" << endl;
	sfst << "  Running time:           "
			<< (float) (te - ts) / CLOCKS_PER_SEC << " [sec]" << endl;

	exit(0);
}
