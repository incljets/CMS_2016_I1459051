// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include <tuple> //< remove from R2.6

namespace Rivet {


  // Inclusive jet pT
  class CMS_2016_I1459051 : public Analysis {
  public:

    // Jet radii
    static constexpr array<double,6> JET_RADII = {0.3, 0.4, 0.5, 0.6, 0.7, 1.0};
    // |y| bin edges
    static constexpr array<double,8> RAP_BINEDGES = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.7};


    // Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2016_I1459051);


    // Book histograms and initialize projections:
    void init() {
      const FinalState fs;

      for (double R : JET_RADII) {

        // Initialize the projectors:
        const size_t iR = size_t(10 * R);
        declare(FastJets(fs, FastJets::ANTIKT, R), "JetsAK" + toString(iR));

        // Book histograms:
        for (size_t iy = 0; iy < RAP_BINEDGES.size()-1; ++iy) {
          const string hsuff = toString(iR) + "_y" + toString(10*RAP_BINEDGES[iy]) + "_" + toString(10*RAP_BINEDGES[iy]);
          const string hname = "pT_R" + hsuff;
          _hists[make_tuple(iR, iy)] = bookHisto1D(hname, refData(iy+1, 1, 1));
        }
      }

    }


    // Analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      for (double R : JET_RADII) {
        const size_t iR = size_t(10 * R);

        // Get jets
        const FastJets& fj = apply<FastJets>(event, "JetsAK" + toString(iR));
        const Jets& jets = fj.jets(Cuts::ptIn(114*GeV, 2200.0*GeV) && Cuts::absrap < 4.7);

        // Fill histograms
        for (const Jet& j : jets) {

          // pT spectra in |y| bins
          // Safe index decoding from Rivet 2.6.0:
          // const int iy = binIndex(j.absrap(), RAP_BINEDGES);
          // if (iy < 0) continue;
          // For now:
          if (!inRange(j.absrap(), RAP_BINEDGES.front(), RAP_BINEDGES.back())) continue; // outside range
          const size_t iy = std::distance(RAP_BINEDGES.begin(), std::upper_bound(RAP_BINEDGES.begin(), RAP_BINEDGES.end(), j.absrap()));
          _hists[make_tuple(iR, iy)]->fill(j.pT()/GeV, weight);

          // Angularities
          double scalar_pt = 0, sum1020 = 0, sum1010 = 0, sum1005 = 0, sum0000 = 0, sum2000 = 0;
          for (const Particle& p : jets.front().particles()) {
            scalar_pt += p.pT();
            sum1020 += pow(p.pT(), 2.0) * pow(deltaR(p, jets[0]), 1.0);
            sum1010 += pow(p.pT(), 1.0) * pow(deltaR(p, jets[0]), 1.0);
            sum1005 += pow(p.pT(), 0.5) * pow(deltaR(p, jets[0]), 1.0);
            sum0000 += pow(p.pT(), 0.0) * pow(deltaR(p, jets[0]), 0.0);
            sum2000 += pow(p.pT(), 0.0) * pow(deltaR(p, jets[0]), 2.0);
          }
          const double ga1020 = sum1020 / pow(scalar_pt, 2.0) * pow(R, 1.0);
          const double ga1010 = sum1010 / pow(scalar_pt, 1.0) * pow(R, 1.0);
          const double ga1005 = sum1005 / pow(scalar_pt, 0.5) * pow(R, 1.0);
          const double ga0000 = sum0000 / pow(scalar_pt, 0.0) * pow(R, 1.0);
          const double ga2000 = sum2000 / pow(scalar_pt, 0.0) * pow(R, 2.0);

        }
      }

    }


    // Finalize
    void finalize() {
      for (auto k_hptr : _hists)
        scale(k_hptr.second, crossSection()/sumOfWeights()/2.0);
    }


  private:

    map<tuple<size_t,size_t>, Histo1DPtr> _hists;

  };


  // Hook for the plugin system.
  DECLARE_RIVET_PLUGIN(CMS_2016_I1459051);


}
