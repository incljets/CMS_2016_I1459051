// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include <tuple> //< remove from R2.6 onward

// NOTES FROM JOEY
// pTjet > 30 GeV/c; upper limit 500 GeV/c for Higgs and 1 TeV for Z
// Rjet = 0.3,0.4,0.5,0.6,0.7,1.0
// plot, pT_leadjet, pT_2nd_leadjet, pT_3rd_leadjet
//     -10 GeV/c bins, which can later be combined
// njet distribution; inclusive and exclusive
// |yjet|<4.5
// |yH(Z)|<2.4 (Higgs(Z) stable)
// plot yH(Z)
// plot yjet distribution for 1st,2nd,3rd jets
// plot pTH(Z)
//     -10 GeV/c bins, same upper limits as above
// --
// for y_leadjet=0-1,1-2,2-3,3-4
// plot pTleadjet
// --
// deltaR_jet-jet>Rjet; jet algorithm should do this automatically
// do we want to require deltaR_jet-H(Z)>Rjet? No

namespace Rivet {


  // Inclusive jet pT
  class LH2017_ZHJETS : public Analysis {
  public:

    // Jet radii
    // static constexpr array<double,6>
    double JET_RADII[6] = {0.3, 0.4, 0.5, 0.6, 0.7, 1.0};

    // Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(LH2017_ZHJETS);


    // Book histograms and initialize projections:
    void init() {
      const FinalState fs;
      declare(fs, "FS");

      for (double R : JET_RADII) {

        // Initialize the projections
        const size_t iR = size_t(10 * R);
        declare(FastJets(fs, FastJets::ANTIKT, R, JetAlg::DECAY_MUONS), "JetsAK" + toString(iR));

        // Book histograms using this suffix
        const string hsuff = "_R" + string(iR < 10 ? "0" : "") + toString(iR);

        // Multiplicity histograms
        _jhists[make_tuple(iR, "njet_excl")] = bookHisto1D("njet_excl" + hsuff, 10, -0.5, 9.5);
        /// @todo Prefer to compute the inclusive spectrum from the excl one in the finalize()?
        // _jhists[make_tuple(iR, "njet_incl")] = bookHisto1D("njet_incl" + hsuff, 10, -0.5, 9.5);

        // Jet pT and rapidity spectra
        for (size_t ijet = 1; ijet <= 3; ++ijet) {
          const string hpre = "J" + toString(ijet);
          _jhists[make_tuple(iR, hpre + "_pT")] = bookHisto1D(hpre + "_pT" + hsuff, 100, 0, 1000);
          _jhists[make_tuple(iR, hpre + "_y")] = bookHisto1D(hpre + "_y" + hsuff, 50, 0, 5);
        }

        if (true){  // for the scope
          const string hpre = "J_incl";
          _jhists[make_tuple(iR, hpre + "_pT")] = bookHisto1D(hpre + "_pT" + hsuff, 100, 0, 1000);
          _jhists[make_tuple(iR, hpre + "_y")] = bookHisto1D(hpre + "_y" + hsuff, 50, 0, 5);
        }

        if (true){  // for the scope
	  _jhists_av[make_tuple(iR, "av_NJet_vs_ptlead" )]= bookProfile1D("av_NJet_vs_ptlead"  + hsuff, 100, 0, 1000);
          _jhists_av[make_tuple(iR, "av_pt_vs_Njet" )]    = bookProfile1D("av_pt_vs_Njet"      + hsuff, 10, -0.5,9.5);
        }

        // Lead jet pT spectra in |y| bins 0-1-2-3-4
        for (size_t iy = 0; iy <= 3; ++iy) {
          const string hpre = "J1dy" + toString(iy);
          _jhists[make_tuple(iR, hpre + "_pT")] = bookHisto1D(hpre + "_pT" + hsuff, 100, 0, 1000);
        }

        // Angularities (treat multiplicity differently)
        _jhists[make_tuple(iR, "GA0000")] = bookHisto1D("GA0000"+hsuff, 151, -0.5, 150.5);
        for (const string& s : {"GA1020", "GA1010", "GA1005", "GA2000"}) { //< without GA0000 = multiplicity
          _jhists[make_tuple(iR, s)] = bookHisto1D(s+hsuff, 200, 0.0, 1.0); //< ranges not quite [0,1] since no WTA jet axis
        }

      }

      // Boson pT and rapidity spectra
      _xhists["XpT"] = bookHisto1D("XpT", 100, 0, 1000);
      _xhists["Xy"] = bookHisto1D("Xy", 25, 0, 2.5);

    }


    // Analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      // Get Higgs or Z
      // const Particles bosons = apply<FinalState>(event, "FS") //< assuming status == 1
      //   .particles(Cuts::pid == PID::ZBOSON || Cuts::pid == PID::HIGGS);
      const Particles zs = event.allParticles(lastParticleWith(Cuts::pid == PID::ZBOSON));
      const Particles hs = event.allParticles(lastParticleWith(Cuts::pid == PID::HIGGS));
      const Particles bosons = zs + hs;
      if (bosons.size() > 1) {
        MSG_WARNING("More than one stable Z/H found... skipping event");
        vetoEvent;
      }
      // Fill boson pT and |y| spectra
      if (!bosons.empty()) {
        const Particle& boson = bosons.front();
        if (boson.absrap()>2.4)vetoEvent;
        _xhists["XpT"]->fill(boson.pT()/GeV, weight);
        _xhists["Xy"]->fill(boson.absrap(), weight);
      }

      for (double R : JET_RADII) {
        const size_t iR = size_t(10 * R);

        // Get jets
        const FastJets& fj = apply<FastJets>(event, "JetsAK" + toString(iR));
        const Jets& jets = fj.jetsByPt(Cuts::pT > 30*GeV && Cuts::absrap < 4.5);

        // Jet multiplicity histograms
        _jhists[make_tuple(iR, "njet_excl")]->fill(jets.size(), weight);



        // for (size_t ijet = 1; ijet <= jets.size(); ++ijet) {
        //   _jhists[make_tuple(iR, "njet_incl")]->fill(ijet, weight);
        // }

        // Need some R-jets from here on
        if (jets.empty()) continue;

        // Lead pT spectra in |y| bins
        const Jet& j1 = jets.front();
        auto httmp=0.;
        for (const Jet& j : jets) httmp+=j.pT()/GeV;
 
        _jhists_av[make_tuple(iR, "av_NJet_vs_ptlead" )]->fill( j1.pT()/GeV ,jets.size(), weight  );
        _jhists_av[make_tuple(iR, "av_pt_vs_Njet" )]    ->fill( jets.size() ,httmp/jets.size() , weight );

        const double y1 = j1.absrap();
        const size_t iy = y1 < 1 ? 0 : (y1 < 2 ? 1 : (y1 < 3 ? 2 : (y1 < 4 ? 3 : 4)));
        if (iy < 4) _jhists[make_tuple(iR, "J1dy" + toString(iy) + "_pT")]->fill(j1.pT()/GeV, weight);
        // if (y1 < 1) {
        //   _jhists[make_tuple(iR, "J1dy0_pT")]->fill(j1.pT()/GeV, weight);
        // } else if (y1 < 2) {
        //   _jhists[make_tuple(iR, "J1dy1_pT")]->fill(j1.pT()/GeV, weight);
        // } else if (y1 < 3) {
        //   _jhists[make_tuple(iR, "J1dy2_pT")]->fill(j1.pT()/GeV, weight);
        // } else if (y1 < 4) {
        //   _jhists[make_tuple(iR, "J1dy3_pT")]->fill(j1.pT()/GeV, weight);
        // }

        size_t ijet = 0;
        for (const Jet& j : jets) {
          ijet += 1;

          // // Everything in |y| bins
          // // Safe index decoding from Rivet 2.6.0:
          // // const int iy = binIndex(j.absrap(), RAP_BINEDGES);
          // // if (iy < 0) continue;
          // // For now:
          // if (!inRange(j.absrap(), RAP_BINEDGES.front(), RAP_BINEDGES.back())) continue; // outside range
          // const size_t iy = std::distance(RAP_BINEDGES.begin(), std::upper_bound(RAP_BINEDGES.begin(), RAP_BINEDGES.end(), j.absrap()));

          // Jet pT and rapidity spectra
          if (ijet <= 3) {
            const string hpre = "J" + toString(ijet);
            _jhists[make_tuple(iR, hpre+"_pT")]->fill(j.pT()/GeV, weight);
            _jhists[make_tuple(iR, hpre+"_y")]->fill(j.absrap(), weight);
          }

          if ( true ){ // for the scope 
          const string hpre = "J_incl"; 
          _jhists[make_tuple(iR, hpre+"_pT")]->fill(j.pT()/GeV, weight);                                                                                                          
          _jhists[make_tuple(iR, hpre+"_y")]->fill(j.absrap(), weight);
          }

          
          // Angularities
          /// @todo The GAs are computed across all jets -- right?
          double scalar_pt = 0; //scalar_pt2 = 0;
          double sum1020 = 0, sum1010 = 0, sum1005 = 0, sum0000 = 0, sum2000 = 0;
          for (const Particle& p : j.particles()) {
            const double pt = p.pT();
            const double dr = deltaR(p, j);
            scalar_pt += pt;
            //scalar_pt2 += sqr(pt);
            sum1020 += pt      * sqr(dr);
            sum1010 += pt      * dr;
            sum1005 += pt      * sqrt(dr);
            sum0000 += 1       * 1;
            sum2000 += sqr(pt) * 1;
          }
          const double ga1020 = sum1020 / scalar_pt      / sqr(R);
          const double ga1010 = sum1010 / scalar_pt      / R;
          const double ga1005 = sum1005 / scalar_pt      / sqrt(R);
          const double ga0000 = sum0000 / 1.             / 1.;
          const double ga2000 = sum2000 / sqr(scalar_pt) / 1.;
          //
          if (ga2000 > 1) MSG_INFO("ga2000 > 1: " << sum2000 << " / " << sqr(scalar_pt) << " = " << ga2000);
          //
          _jhists[make_tuple(iR, "GA1020")]->fill(ga1020, weight);
          _jhists[make_tuple(iR, "GA1010")]->fill(ga1010, weight);
          _jhists[make_tuple(iR, "GA1005")]->fill(ga1005, weight);
          _jhists[make_tuple(iR, "GA0000")]->fill(ga0000, weight);
          _jhists[make_tuple(iR, "GA2000")]->fill(ga2000, weight);

        }
      }

    }


    // Finalize
    void finalize() {
      for (auto k_hptr : _jhists)
        scale(k_hptr.second, crossSection()/sumOfWeights());
      for (auto k_hptr : _xhists)
        scale(k_hptr.second, crossSection()/sumOfWeights());


      /// @todo Compute inclusive Njet spectrum here
    }


  private:

    map<tuple<size_t,string>, Histo1DPtr> _jhists;
    map<tuple<size_t,string>, Profile1DPtr> _jhists_av;

    map<string, Histo1DPtr> _xhists;

  };


  // Hook for the plugin system.
  DECLARE_RIVET_PLUGIN(LH2017_ZHJETS);


}
