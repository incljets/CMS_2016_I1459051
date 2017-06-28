// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {

  // Inclusive jet pT
  class CMS_2016_I1459051 : public Analysis {
  public:

    // Constructor
    CMS_2016_I1459051() : Analysis("CMS_2016_I1459051") {}


    // Book histograms and initialize projections:
    void init() {
      const FinalState fs;

      // Initialize the projectors:
      addProjection(FastJets(fs, FastJets::ANTIKT, 0.4),"JetsAK4");
      addProjection(FastJets(fs, FastJets::ANTIKT, 0.7),"JetsAK7");

      // Book histograms:
      _hist_sigmaAK4.addHistogram(0.0, 0.5, bookHisto1D(8, 1, 1));
      _hist_sigmaAK4.addHistogram(0.5, 1.0, bookHisto1D(9, 1, 1));
      _hist_sigmaAK4.addHistogram(1.0, 1.5, bookHisto1D(10, 1, 1));
      _hist_sigmaAK4.addHistogram(1.5, 2.0, bookHisto1D(11, 1, 1));
      _hist_sigmaAK4.addHistogram(2.0, 2.5, bookHisto1D(12, 1, 1));
      _hist_sigmaAK4.addHistogram(2.5, 3.0, bookHisto1D(13, 1, 1));
      
      _hist_sigmaAK7.addHistogram(0.0, 0.5, bookHisto1D(1, 1, 1));
      _hist_sigmaAK7.addHistogram(0.5, 1.0, bookHisto1D(2, 1, 1));
      _hist_sigmaAK7.addHistogram(1.0, 1.5, bookHisto1D(3, 1, 1));
      _hist_sigmaAK7.addHistogram(1.5, 2.0, bookHisto1D(4, 1, 1));
      _hist_sigmaAK7.addHistogram(2.0, 2.5, bookHisto1D(5, 1, 1));
      _hist_sigmaAK7.addHistogram(2.5, 3.0, bookHisto1D(6, 1, 1));
      
      _hist_sigmaAK7Forward = bookHisto1D(7, 1, 1);
      _hist_sigmaAK4Forward = bookHisto1D(14, 1, 1);
      
    }

    // Analysis
    void analyze(const Event &event) {
      const double weight = event.weight();
      const FastJets& fjAK4 = applyProjection<FastJets>(event,"JetsAK4");
      const FastJets& fjAK7 = applyProjection<FastJets>(event,"JetsAK7");
      const Jets& jetsAK4 = fjAK4.jets(Cuts::ptIn(114*GeV, 2200.0*GeV) && Cuts::absrap < 4.7);
      const Jets& jetsAK7 = fjAK7.jets(Cuts::ptIn(114*GeV, 2200.0*GeV) && Cuts::absrap < 4.7);

      // Fill the relevant histograms:
      foreach(const Jet& j, jetsAK4) {
        _hist_sigmaAK4.fill(j.absrap(), j.pT(), weight);
	if(j.absrap()>3.2 && j.absrap()<4.7) _hist_sigmaAK4Forward->fill(j.pT(), weight);
      }

      foreach(const Jet& j, jetsAK7) {
        _hist_sigmaAK7.fill(j.absrap(), j.pT(), weight);
	if(j.absrap()>3.2 && j.absrap()<4.7) _hist_sigmaAK7Forward->fill(j.pT(), weight);
      }

    }

    // Finalize
    void finalize() {
      _hist_sigmaAK4.scale(crossSection()/sumOfWeights()/2.0, this);
      _hist_sigmaAK7.scale(crossSection()/sumOfWeights()/2.0, this);
      scale(_hist_sigmaAK4Forward,crossSection()/sumOfWeights()/3.0);
      scale(_hist_sigmaAK7Forward,crossSection()/sumOfWeights()/3.0);
    }

  private:
    BinnedHistogram<double> _hist_sigmaAK4;
    BinnedHistogram<double> _hist_sigmaAK7;
    Histo1DPtr _hist_sigmaAK4Forward;
    Histo1DPtr _hist_sigmaAK7Forward;
  };

  // This global object acts as a hook for the plugin system.
  DECLARE_RIVET_PLUGIN(CMS_2016_I1459051);

}
