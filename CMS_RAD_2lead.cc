// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  /// Inclusive jet pT at 13 TeV
  class CMS_RAD_2lead : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_RAD_2lead);


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
      // AK7
      _hist_sigmaAK7.addHistogram(0.0, 0.5, bookHisto1D("d01-x01-y01-AK7", refData(1,1,1)));
      _hist_sigmaAK7.addHistogram(0.5, 1.0, bookHisto1D("d02-x01-y01-AK7", refData(2,1,1)));
      _hist_sigmaAK7.addHistogram(1.0, 1.5, bookHisto1D("d03-x01-y01-AK7", refData(3,1,1)));
      _hist_sigmaAK7.addHistogram(1.5, 2.0, bookHisto1D("d04-x01-y01-AK7", refData(4,1,1)));
      _hist_sigmaAK7.addHistogram(2.0, 2.5, bookHisto1D("d05-x01-y01-AK7", refData(5,1,1)));
      _hist_sigmaAK7.addHistogram(2.5, 3.0, bookHisto1D("d06-x01-y01-AK7", refData(6,1,1)));
      _hist_sigmaAK7Forward =               bookHisto1D("d07-x01-y01-AK7", refData(7,1,1));
      // AK4
      _hist_sigmaAK2.addHistogram(0.0, 0.5, bookHisto1D("d01-x01-y01-AK2", refData(1,1,1)));
      _hist_sigmaAK2.addHistogram(0.5, 1.0, bookHisto1D("d02-x01-y01-AK2", refData(2,1,1)));
      _hist_sigmaAK2.addHistogram(1.0, 1.5, bookHisto1D("d03-x01-y01-AK2", refData(3,1,1)));
      _hist_sigmaAK2.addHistogram(1.5, 2.0, bookHisto1D("d04-x01-y01-AK2", refData(4,1,1)));
      _hist_sigmaAK2.addHistogram(2.0, 2.5, bookHisto1D("d05-x01-y01-AK2", refData(5,1,1)));
      _hist_sigmaAK2.addHistogram(2.5, 3.0, bookHisto1D("d06-x01-y01-AK2", refData(6,1,1)));
      _hist_sigmaAK2Forward =               bookHisto1D("d07-x01-y01-AK2", refData(7,1,1));

      _hist_sigmaAK3.addHistogram(0.0, 0.5, bookHisto1D("d01-x01-y01-AK3", refData(1,1,1)));
      _hist_sigmaAK3.addHistogram(0.5, 1.0, bookHisto1D("d02-x01-y01-AK3", refData(2,1,1)));
      _hist_sigmaAK3.addHistogram(1.0, 1.5, bookHisto1D("d03-x01-y01-AK3", refData(3,1,1)));
      _hist_sigmaAK3.addHistogram(1.5, 2.0, bookHisto1D("d04-x01-y01-AK3", refData(4,1,1)));
      _hist_sigmaAK3.addHistogram(2.0, 2.5, bookHisto1D("d05-x01-y01-AK3", refData(5,1,1)));
      _hist_sigmaAK3.addHistogram(2.5, 3.0, bookHisto1D("d06-x01-y01-AK3", refData(6,1,1)));
      _hist_sigmaAK3Forward =               bookHisto1D("d07-x01-y01-AK3", refData(7,1,1));

      _hist_sigmaAK4.addHistogram(0.0, 0.5, bookHisto1D("d01-x01-y01-AK4", refData(1,1,1)));
      _hist_sigmaAK4.addHistogram(0.5, 1.0, bookHisto1D("d02-x01-y01-AK4", refData(2,1,1)));
      _hist_sigmaAK4.addHistogram(1.0, 1.5, bookHisto1D("d03-x01-y01-AK4", refData(3,1,1)));
      _hist_sigmaAK4.addHistogram(1.5, 2.0, bookHisto1D("d04-x01-y01-AK4", refData(4,1,1)));
      _hist_sigmaAK4.addHistogram(2.0, 2.5, bookHisto1D("d05-x01-y01-AK4", refData(5,1,1)));
      _hist_sigmaAK4.addHistogram(2.5, 3.0, bookHisto1D("d06-x01-y01-AK4", refData(6,1,1)));
      _hist_sigmaAK4Forward =               bookHisto1D("d07-x01-y01-AK4", refData(7,1,1));

      _hist_sigmaAK5.addHistogram(0.0, 0.5, bookHisto1D("d01-x01-y01-AK5", refData(1,1,1)));
      _hist_sigmaAK5.addHistogram(0.5, 1.0, bookHisto1D("d02-x01-y01-AK5", refData(2,1,1)));
      _hist_sigmaAK5.addHistogram(1.0, 1.5, bookHisto1D("d03-x01-y01-AK5", refData(3,1,1)));
      _hist_sigmaAK5.addHistogram(1.5, 2.0, bookHisto1D("d04-x01-y01-AK5", refData(4,1,1)));
      _hist_sigmaAK5.addHistogram(2.0, 2.5, bookHisto1D("d05-x01-y01-AK5", refData(5,1,1)));
      _hist_sigmaAK5.addHistogram(2.5, 3.0, bookHisto1D("d06-x01-y01-AK5", refData(6,1,1)));
      _hist_sigmaAK5Forward =               bookHisto1D("d07-x01-y01-AK5", refData(7,1,1));

      _hist_sigmaAK6.addHistogram(0.0, 0.5, bookHisto1D("d01-x01-y01-AK6", refData(1,1,1)));
      _hist_sigmaAK6.addHistogram(0.5, 1.0, bookHisto1D("d02-x01-y01-AK6", refData(2,1,1)));
      _hist_sigmaAK6.addHistogram(1.0, 1.5, bookHisto1D("d03-x01-y01-AK6", refData(3,1,1)));
      _hist_sigmaAK6.addHistogram(1.5, 2.0, bookHisto1D("d04-x01-y01-AK6", refData(4,1,1)));
      _hist_sigmaAK6.addHistogram(2.0, 2.5, bookHisto1D("d05-x01-y01-AK6", refData(5,1,1)));
      _hist_sigmaAK6.addHistogram(2.5, 3.0, bookHisto1D("d06-x01-y01-AK6", refData(6,1,1)));
      _hist_sigmaAK6Forward =               bookHisto1D("d07-x01-y01-AK6", refData(7,1,1));

      _hist_sigmaAK8.addHistogram(0.0, 0.5, bookHisto1D("d01-x01-y01-AK8", refData(1,1,1)));
      _hist_sigmaAK8.addHistogram(0.5, 1.0, bookHisto1D("d02-x01-y01-AK8", refData(2,1,1)));
      _hist_sigmaAK8.addHistogram(1.0, 1.5, bookHisto1D("d03-x01-y01-AK8", refData(3,1,1)));
      _hist_sigmaAK8.addHistogram(1.5, 2.0, bookHisto1D("d04-x01-y01-AK8", refData(4,1,1)));
      _hist_sigmaAK8.addHistogram(2.0, 2.5, bookHisto1D("d05-x01-y01-AK8", refData(5,1,1)));
      _hist_sigmaAK8.addHistogram(2.5, 3.0, bookHisto1D("d06-x01-y01-AK8", refData(6,1,1)));
      _hist_sigmaAK8Forward =               bookHisto1D("d07-x01-y01-AK8", refData(7,1,1));

      _hist_sigmaAK9.addHistogram(0.0, 0.5, bookHisto1D("d01-x01-y01-AK9", refData(1,1,1)));
      _hist_sigmaAK9.addHistogram(0.5, 1.0, bookHisto1D("d02-x01-y01-AK9", refData(2,1,1)));
      _hist_sigmaAK9.addHistogram(1.0, 1.5, bookHisto1D("d03-x01-y01-AK9", refData(3,1,1)));
      _hist_sigmaAK9.addHistogram(1.5, 2.0, bookHisto1D("d04-x01-y01-AK9", refData(4,1,1)));
      _hist_sigmaAK9.addHistogram(2.0, 2.5, bookHisto1D("d05-x01-y01-AK9", refData(5,1,1)));
      _hist_sigmaAK9.addHistogram(2.5, 3.0, bookHisto1D("d06-x01-y01-AK9", refData(6,1,1)));
      _hist_sigmaAK9Forward =               bookHisto1D("d07-x01-y01-AK9", refData(7,1,1));

      _hist_sigmaAK10.addHistogram(0.0, 0.5, bookHisto1D("d01-x01-y01-AK10", refData(1,1,1)));
      _hist_sigmaAK10.addHistogram(0.5, 1.0, bookHisto1D("d02-x01-y01-AK10", refData(2,1,1)));
      _hist_sigmaAK10.addHistogram(1.0, 1.5, bookHisto1D("d03-x01-y01-AK10", refData(3,1,1)));
      _hist_sigmaAK10.addHistogram(1.5, 2.0, bookHisto1D("d04-x01-y01-AK10", refData(4,1,1)));
      _hist_sigmaAK10.addHistogram(2.0, 2.5, bookHisto1D("d05-x01-y01-AK10", refData(5,1,1)));
      _hist_sigmaAK10.addHistogram(2.5, 3.0, bookHisto1D("d06-x01-y01-AK10", refData(6,1,1)));
      _hist_sigmaAK10Forward =               bookHisto1D("d07-x01-y01-AK10", refData(7,1,1));

      _hist_sigmaAK11.addHistogram(0.0, 0.5, bookHisto1D("d01-x01-y01-AK11", refData(1,1,1)));
      _hist_sigmaAK11.addHistogram(0.5, 1.0, bookHisto1D("d02-x01-y01-AK11", refData(2,1,1)));
      _hist_sigmaAK11.addHistogram(1.0, 1.5, bookHisto1D("d03-x01-y01-AK11", refData(3,1,1)));
      _hist_sigmaAK11.addHistogram(1.5, 2.0, bookHisto1D("d04-x01-y01-AK11", refData(4,1,1)));
      _hist_sigmaAK11.addHistogram(2.0, 2.5, bookHisto1D("d05-x01-y01-AK11", refData(5,1,1)));
      _hist_sigmaAK11.addHistogram(2.5, 3.0, bookHisto1D("d06-x01-y01-AK11", refData(6,1,1)));
      _hist_sigmaAK11Forward =               bookHisto1D("d07-x01-y01-AK11", refData(7,1,1));


    }


    /// Per-event analysis
    void analyze(const Event &event) {

      const double weight = event.weight();
      int  count=0;
      // AK4 jets
      const FastJets& fjAK4 = applyProjection<FastJets>(event, "JetsAK4");
      const Jets& jetsAK4 = fjAK4.jetsByPt(Cuts::ptIn(114*GeV, 2200.0*GeV) && Cuts::absrap < 4.7);
      count=0;for (const Jet& j : jetsAK4) {
        _hist_sigmaAK4.fill(j.absrap(), j.pT(), weight);
        if (inRange(j.absrap(), 3.2, 4.7)) _hist_sigmaAK4Forward->fill(j.pT(), weight);count++;if(count>1)break;
      }

      // AK7 jets
      const FastJets& fjAK7 = applyProjection<FastJets>(event, "JetsAK7");
      const Jets& jetsAK7 = fjAK7.jetsByPt(Cuts::ptIn(114*GeV, 2200.0*GeV) && Cuts::absrap < 4.7);
      count=0;for (const Jet& j : jetsAK7) {
        _hist_sigmaAK7.fill(j.absrap(), j.pT(), weight);
        if (inRange(j.absrap(), 3.2, 4.7)) _hist_sigmaAK7Forward->fill(j.pT(), weight);count++;if(count>1)break;
      }



      const FastJets& fjAK2 = applyProjection<FastJets>(event, "JetsAK2");
      const Jets& jetsAK2 = fjAK2.jetsByPt(Cuts::ptIn(114*GeV, 2200.0*GeV) && Cuts::absrap < 4.7);
      count=0;for (const Jet& j : jetsAK2) {
        _hist_sigmaAK2.fill(j.absrap(), j.pT(), weight);
        if (inRange(j.absrap(), 3.2, 4.7)) _hist_sigmaAK2Forward->fill(j.pT(), weight);count++;if(count>1)break;
      }


      const FastJets& fjAK3 = applyProjection<FastJets>(event, "JetsAK3");
      const Jets& jetsAK3 = fjAK3.jetsByPt(Cuts::ptIn(114*GeV, 2200.0*GeV) && Cuts::absrap < 4.7);
      count=0;for (const Jet& j : jetsAK3) {
        _hist_sigmaAK3.fill(j.absrap(), j.pT(), weight);
        if (inRange(j.absrap(), 3.2, 4.7)) _hist_sigmaAK3Forward->fill(j.pT(), weight);count++;if(count>1)break;
      }


      const FastJets& fjAK5 = applyProjection<FastJets>(event, "JetsAK5");
      const Jets& jetsAK5 = fjAK5.jetsByPt(Cuts::ptIn(114*GeV, 2200.0*GeV) && Cuts::absrap < 4.7);
      count=0;for (const Jet& j : jetsAK5) {
        _hist_sigmaAK5.fill(j.absrap(), j.pT(), weight);
        if (inRange(j.absrap(), 3.2, 4.7)) _hist_sigmaAK5Forward->fill(j.pT(), weight);count++;if(count>1)break;
      }


      const FastJets& fjAK6 = applyProjection<FastJets>(event, "JetsAK6");
      const Jets& jetsAK6 = fjAK6.jetsByPt(Cuts::ptIn(114*GeV, 2200.0*GeV) && Cuts::absrap < 4.7);
      count=0;for (const Jet& j : jetsAK6) {
        _hist_sigmaAK6.fill(j.absrap(), j.pT(), weight);
        if (inRange(j.absrap(), 3.2, 4.7)) _hist_sigmaAK6Forward->fill(j.pT(), weight);count++;if(count>1)break;
      }


      const FastJets& fjAK8 = applyProjection<FastJets>(event, "JetsAK8");
      const Jets& jetsAK8 = fjAK8.jetsByPt(Cuts::ptIn(114*GeV, 2200.0*GeV) && Cuts::absrap < 4.7);
      count=0;for (const Jet& j : jetsAK8) {
        _hist_sigmaAK8.fill(j.absrap(), j.pT(), weight);
        if (inRange(j.absrap(), 3.2, 4.7)) _hist_sigmaAK8Forward->fill(j.pT(), weight);count++;if(count>1)break;
      }


      const FastJets& fjAK9 = applyProjection<FastJets>(event, "JetsAK9");
      const Jets& jetsAK9 = fjAK9.jetsByPt(Cuts::ptIn(114*GeV, 2200.0*GeV) && Cuts::absrap < 4.7);
      count=0;for (const Jet& j : jetsAK9) {
        _hist_sigmaAK9.fill(j.absrap(), j.pT(), weight);
        if (inRange(j.absrap(), 3.2, 4.7)) _hist_sigmaAK9Forward->fill(j.pT(), weight);count++;if(count>1)break;
      }


      const FastJets& fjAK10 = applyProjection<FastJets>(event, "JetsAK10");
      const Jets& jetsAK10 = fjAK10.jetsByPt(Cuts::ptIn(114*GeV, 2200.0*GeV) && Cuts::absrap < 4.7);
      count=0;for (const Jet& j : jetsAK10) {
        _hist_sigmaAK10.fill(j.absrap(), j.pT(), weight);
        if (inRange(j.absrap(), 3.2, 4.7)) _hist_sigmaAK10Forward->fill(j.pT(), weight);count++;if(count>1)break;
      }


      const FastJets& fjAK11 = applyProjection<FastJets>(event, "JetsAK11");
      const Jets& jetsAK11 = fjAK11.jetsByPt(Cuts::ptIn(114*GeV, 2200.0*GeV) && Cuts::absrap < 4.7);
      count=0;for (const Jet& j : jetsAK11) {
        _hist_sigmaAK11.fill(j.absrap(), j.pT(), weight);
        if (inRange(j.absrap(), 3.2, 4.7)) _hist_sigmaAK11Forward->fill(j.pT(), weight);count++;if(count>1)break;
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


      scale(_hist_sigmaAK2Forward,crossSection()/sumOfWeights()/3.0);
      scale(_hist_sigmaAK3Forward,crossSection()/sumOfWeights()/3.0);
      scale(_hist_sigmaAK4Forward,crossSection()/sumOfWeights()/3.0);
      scale(_hist_sigmaAK5Forward,crossSection()/sumOfWeights()/3.0);
      scale(_hist_sigmaAK6Forward,crossSection()/sumOfWeights()/3.0);
      scale(_hist_sigmaAK7Forward,crossSection()/sumOfWeights()/3.0);
      scale(_hist_sigmaAK8Forward,crossSection()/sumOfWeights()/3.0);
      scale(_hist_sigmaAK9Forward,crossSection()/sumOfWeights()/3.0);
      scale(_hist_sigmaAK10Forward,crossSection()/sumOfWeights()/3.0);
      scale(_hist_sigmaAK11Forward,crossSection()/sumOfWeights()/3.0);

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
    Histo1DPtr _hist_sigmaAK2Forward;
    Histo1DPtr _hist_sigmaAK3Forward;
    Histo1DPtr _hist_sigmaAK4Forward;
    Histo1DPtr _hist_sigmaAK5Forward;
    Histo1DPtr _hist_sigmaAK6Forward;
    Histo1DPtr _hist_sigmaAK7Forward;
    Histo1DPtr _hist_sigmaAK8Forward;
    Histo1DPtr _hist_sigmaAK9Forward;
    Histo1DPtr _hist_sigmaAK10Forward;
    Histo1DPtr _hist_sigmaAK11Forward;




    //@}

  };


  // This global object acts as a hook for the plugin system.
  DECLARE_RIVET_PLUGIN(CMS_RAD_2lead);

}
