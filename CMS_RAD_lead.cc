// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  /// Inclusive jet pT at 13 TeV
  class CMS_RAD_lead : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_RAD_lead);


    /// Book histograms and initialize projections:
    void init() {

      // Initialize the projections
      const FinalState fs;
      declare(FastJets(fs, FastJets::ANTIKT, 0.2), "JetsAK2");
      declare(FastJets(fs, FastJets::ANTIKT, 0.3), "JetsAK3");
      declare(FastJets(fs, FastJets::ANTIKT, 0.4), "JetsAK4");
      declare(FastJets(fs, FastJets::ANTIKT, 0.5), "JetsAK5");
      declare(FastJets(fs, FastJets::ANTIKT, 0.6), "JetsAK6");
      declare(FastJets(fs, FastJets::ANTIKT, 0.7), "JetsAK7");
      declare(FastJets(fs, FastJets::ANTIKT, 0.8), "JetsAK8");
      declare(FastJets(fs, FastJets::ANTIKT, 0.9), "JetsAK9");
      declare(FastJets(fs, FastJets::ANTIKT, 1.0), "JetsAK10");
      declare(FastJets(fs, FastJets::ANTIKT, 1.1), "JetsAK11");


      // Book sets of histograms, binned in absolute rapidity
      _hist_sigmaAK7.addHistogram(0.0, 3.5, bookHisto1D("d01-x01-y01-AK7", refData(1,1,1)));
      _hist_sigmaAK2.addHistogram(0.0, 3.5, bookHisto1D("d01-x01-y01-AK2", refData(1,1,1)));
      _hist_sigmaAK3.addHistogram(0.0, 3.5, bookHisto1D("d01-x01-y01-AK3", refData(1,1,1)));
      _hist_sigmaAK4.addHistogram(0.0, 3.5, bookHisto1D("d01-x01-y01-AK4", refData(1,1,1)));
      _hist_sigmaAK5.addHistogram(0.0, 3.5, bookHisto1D("d01-x01-y01-AK5", refData(1,1,1)));
      _hist_sigmaAK6.addHistogram(0.0, 3.5, bookHisto1D("d01-x01-y01-AK6", refData(1,1,1)));
      _hist_sigmaAK8.addHistogram(0.0, 3.5, bookHisto1D("d01-x01-y01-AK8", refData(1,1,1)));
      _hist_sigmaAK9.addHistogram(0.0, 3.5, bookHisto1D("d01-x01-y01-AK9", refData(1,1,1)));
      _hist_sigmaAK10.addHistogram(0.0, 3.5, bookHisto1D("d01-x01-y01-AK10", refData(1,1,1)));
      _hist_sigmaAK11.addHistogram(0.0, 3.5, bookHisto1D("d01-x01-y01-AK11", refData(1,1,1)));


    }


    /// Per-event analysis
    void analyze(const Event &event) {

      const double weight = event.weight();

      // AK4 jets
      const FastJets& fjAK4 = applyProjection<FastJets>(event, "JetsAK4");
      const Jets& jetsAK4 = fjAK4.jetsByPt(Cuts::ptIn(114*GeV, 2200.0*GeV) && Cuts::absrap < 3.5);
      for (const Jet& j : jetsAK4) {
        _hist_sigmaAK4.fill(j.absrap(), j.pT(), weight);
      break;}

      // AK7 jets
      const FastJets& fjAK7 = applyProjection<FastJets>(event, "JetsAK7");
      const Jets& jetsAK7 = fjAK7.jetsByPt(Cuts::ptIn(114*GeV, 2200.0*GeV) && Cuts::absrap < 3.5);
      for (const Jet& j : jetsAK7) {
        _hist_sigmaAK7.fill(j.absrap(), j.pT(), weight);
      break;}



      const FastJets& fjAK2 = applyProjection<FastJets>(event, "JetsAK2");
      const Jets& jetsAK2 = fjAK2.jetsByPt(Cuts::ptIn(114*GeV, 2200.0*GeV) && Cuts::absrap < 3.5);
      for (const Jet& j : jetsAK2) {
        _hist_sigmaAK2.fill(j.absrap(), j.pT(), weight);
      break;}


      const FastJets& fjAK3 = applyProjection<FastJets>(event, "JetsAK3");
      const Jets& jetsAK3 = fjAK3.jetsByPt(Cuts::ptIn(114*GeV, 2200.0*GeV) && Cuts::absrap < 3.5);
      for (const Jet& j : jetsAK3) {
        _hist_sigmaAK3.fill(j.absrap(), j.pT(), weight);
      break;}


      const FastJets& fjAK5 = applyProjection<FastJets>(event, "JetsAK5");
      const Jets& jetsAK5 = fjAK5.jetsByPt(Cuts::ptIn(114*GeV, 2200.0*GeV) && Cuts::absrap < 3.5);
      for (const Jet& j : jetsAK5) {
        _hist_sigmaAK5.fill(j.absrap(), j.pT(), weight);
      break;}


      const FastJets& fjAK6 = applyProjection<FastJets>(event, "JetsAK6");
      const Jets& jetsAK6 = fjAK6.jetsByPt(Cuts::ptIn(114*GeV, 2200.0*GeV) && Cuts::absrap < 3.5);
      for (const Jet& j : jetsAK6) {
        _hist_sigmaAK6.fill(j.absrap(), j.pT(), weight);
      break;}


      const FastJets& fjAK8 = applyProjection<FastJets>(event, "JetsAK8");
      const Jets& jetsAK8 = fjAK8.jetsByPt(Cuts::ptIn(114*GeV, 2200.0*GeV) && Cuts::absrap < 3.5);
      for (const Jet& j : jetsAK8) {
        _hist_sigmaAK8.fill(j.absrap(), j.pT(), weight);
      break;}


      const FastJets& fjAK9 = applyProjection<FastJets>(event, "JetsAK9");
      const Jets& jetsAK9 = fjAK9.jetsByPt(Cuts::ptIn(114*GeV, 2200.0*GeV) && Cuts::absrap < 3.5);
      for (const Jet& j : jetsAK9) {
        _hist_sigmaAK9.fill(j.absrap(), j.pT(), weight);
        break;
      }


      const FastJets& fjAK10 = applyProjection<FastJets>(event, "JetsAK10");
      const Jets& jetsAK10 = fjAK10.jetsByPt(Cuts::ptIn(114*GeV, 2200.0*GeV) && Cuts::absrap < 3.5);
      for (const Jet& j : jetsAK10) {
        _hist_sigmaAK10.fill(j.absrap(), j.pT(), weight);
        break;
      }


      const FastJets& fjAK11 = applyProjection<FastJets>(event, "JetsAK11");
      const Jets& jetsAK11 = fjAK11.jetsByPt(Cuts::ptIn(114*GeV, 2200.0*GeV) && Cuts::absrap < 3.5);
      for (const Jet& j : jetsAK11) {
        _hist_sigmaAK11.fill(j.absrap(), j.pT(), weight);
        break;
      }





    }


    // Finalize
    void finalize() {
      /// @todo What is the cross-section unit?
      _hist_sigmaAK2.scale(crossSection()/sumOfWeights()/2.0, this);
      _hist_sigmaAK3.scale(crossSection()/sumOfWeights()/2.0, this);
      _hist_sigmaAK4.scale(crossSection()/sumOfWeights()/2.0, this);
      _hist_sigmaAK5.scale(crossSection()/sumOfWeights()/2.0, this);
      _hist_sigmaAK6.scale(crossSection()/sumOfWeights()/2.0, this);
      _hist_sigmaAK7.scale(crossSection()/sumOfWeights()/2.0, this);
      _hist_sigmaAK8.scale(crossSection()/sumOfWeights()/2.0, this);
      _hist_sigmaAK9.scale(crossSection()/sumOfWeights()/2.0, this);
      _hist_sigmaAK10.scale(crossSection()/sumOfWeights()/2.0, this);
      _hist_sigmaAK11.scale(crossSection()/sumOfWeights()/2.0, this);
    }


    /// @name Histograms
    //@{
    BinnedHistogram<double> _hist_sigmaAK2;
    BinnedHistogram<double> _hist_sigmaAK3;
    BinnedHistogram<double> _hist_sigmaAK4;
    BinnedHistogram<double> _hist_sigmaAK5;
    BinnedHistogram<double> _hist_sigmaAK6;
    BinnedHistogram<double> _hist_sigmaAK7;
    BinnedHistogram<double> _hist_sigmaAK8;
    BinnedHistogram<double> _hist_sigmaAK9;
    BinnedHistogram<double> _hist_sigmaAK10;
    BinnedHistogram<double> _hist_sigmaAK11;

    //@}

  };


  // This global object acts as a hook for the plugin system.
  DECLARE_RIVET_PLUGIN(CMS_RAD_lead);

}
