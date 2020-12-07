// parameters.h

#ifndef PARAMETERS_H
#define PARAMETERS_H

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppGSL)]]

#include <vector>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <ctime>
#include <limits>
#include <omp.h>
#include <Rcpp.h>
using namespace std;

#include "randomizer.h"
#include "distribution.h"
#include "helper.h"
#include "process_spec.h"
#include "changes.h"
#include "parameters.h"
#include "compartment.h"
#include "reporter.h"
#include "observer.h"

struct Matrix;
struct Change;

void ParamSet(Discrete& variable, Rcpp::RObject& value);
void ParamSet(vector<double>& variable, Rcpp::RObject& value);
void ParamSet(Matrix& variable, Rcpp::RObject& value);
void ParamSet(vector<Matrix>& variable, Rcpp::RObject& value);
void ParamSet(vector<ProcessSpec>& variable, Rcpp::RObject& value);

// Population-level parameters for batch of simulations
struct PopulationParameters
{
public:
    PopulationParameters() : needs_recalc(true) {}

    bool needs_recalc;  // does contact matrix need recalculation?
    Matrix cm;          // contact matrix

    Discrete dE;
    Discrete dEa;
    Discrete dIp;       // TODO: any need for these to be age-specific?
    Discrete dIa;
    Discrete dIs;
    Discrete dC;

    vector<double> size;
    vector<double> imm0;
    vector<Matrix> matrices;
    vector<double> contact;
    vector<double> contact_mult;
    vector<double> contact_lowerto;
    vector<double> u;
    vector<double> fIp;
    vector<double> fIa;
    vector<double> fIs;
    vector<double> y;
    vector<double> omega;
    vector<double> rho;
    vector<double> tau;
    vector<double> v;
    vector<double> v12;
    vector<double> v2;
    vector<double> ev;
    vector<double> ei_v;
    vector<double> ed_vi;
    vector<double> ev2;
    vector<double> ei_v2;
    vector<double> ed_vi2;
    vector<double> pi_r;
    vector<double> pd_ri;
    vector<double> wn;
    vector<double> wv;
    vector<double> wv2;
    vector<double> A;
    vector<double> B;
    vector<double> D;
    vector<double> season_A;    // note - size 1 vector
    vector<double> season_T;    // note - size 1 vector
    vector<double> season_phi;  // note - size 1 vector

    vector<double> seed_times;
    Discrete dist_seed_ages;

    ///Observer observer;
    ///vector<ScheduleEntry> schedule;

    string name;
    vector<string> group_names;

    bool Set(Parameters* parent, string& name, Rcpp::RObject& value);
    bool Set(Parameters* parent, string& name, vector<double>& value);

    void Recalculate()
    {
        if (needs_recalc)
        {
            if (matrices.empty())
                throw logic_error("No contact matrices defined.");
            if (matrices.size() != contact.size())
                throw logic_error("Number of contact components not equal to number of matrices.");

            unsigned int ncol = matrices[0].nc;
            unsigned int nrow = matrices[0].x.size() / ncol;

            auto c_mult = [&](int m) { if (!contact_mult.empty()) return contact_mult[m]; return 1.0; };
            auto c_lowerto = [&](int m) { if (!contact_lowerto.empty()) return contact_lowerto[m]; return std::numeric_limits<double>::max(); };

            cm = Matrix(0, nrow, ncol);

            for (unsigned int r = 0; r < nrow; ++r)
                for (unsigned int c = 0; c < ncol; ++c)
                    for (unsigned int m = 0; m < matrices.size(); ++m)
                        cm(r, c) += matrices[m](r, c) * min(contact[m] * c_mult(m), c_lowerto(m));

            needs_recalc = false;
        }
    }
};

struct Parameters
{
public:
    void FilterForRun(unsigned int r);

    string model;
    double time_step;
    double time0;
    double time1;
    unsigned int report_every;
    bool fast_multinomial;
    bool deterministic;

    vector<PopulationParameters> pop;

    vector<ProcessSpec> processes;
    ///Observer observer;
    Matrix travel;
    ChangeSet changes;
};

void SetParameters(Parameters& P, Rcpp::List list, Randomizer& Rand);

#endif
