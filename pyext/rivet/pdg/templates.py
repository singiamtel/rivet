from string import Template
info_T = Template("""\
Name: ${rivetName}
Year: <Insert year of publication>
Summary: $$ ${latexName} $$ branching ratios
Status: VALIDATED NOHEPDATA
Reentrant: true
Authors:
 - Your Name <your@email.address>
References:
 - 'Example: Eur.Phys.J. C76 (2016) no.7, 392'
 - 'Example: DOI:10.1140/epjc/s10052-016-4184-8'
 - 'Example: arXiv:1605.03814'
RunInfo: Any process producing ${name}
Description:
  'Branching ratios from the PDG for testing simulations of $$ ${latexName} $$ decays.'
ValidationInfo:
  'A description of the process used to validate the Rivet code against
  the original experimental analysis. Cut-flow tables and similar information
  are welcome'
Keywords: []
BibKey: Workman:2022ynf
BibTeX: '@article{Workman:2022ynf,
    author = "Workman, R. L. and Others",
    collaboration = "Particle Data Group",
    title = "{Review of Particle Physics}",
    journal = "PTEP",
    volume = "2022",
    pages = "083C01",
    year = "2022"
}'
""")

cc_T = Template("""\
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {

  /// @brief ${name} decays
  class ${rivetName} : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(${rivetName});

    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(Cuts::pid==${pid}), "UFS");

      vector<unsigned int> ids = { ${histIDs} };
      for (unsigned int iloc : ids) {
       	book(_h_mult[iloc], 1, 1, iloc);
      }
      book(_n, "/TMP/n_part");
      ${initialize}
    }

    void findDecayProducts(const Particle& mother, map<int,unsigned int>& cStable,
                           map<int,unsigned int>& cUnStable, unsigned int& nstable,
                           const set<int>& unstable) const {
      for (const Particle& p : mother.children()) {
        int id = p.pid();
        if (p.children().empty() || unstable.find(id) != unstable.end()) {
          nstable +=1;
          cStable[id] += 1;
        }
        else {
          // check particle is not a child or itself, eg copy or from photon radiation
          bool isChild(false);
          for (const Particle& p2 : p.children()) {
            if (p2.pid()==id) {
              isChild = true;
              break;
            }
          }
          if (!isChild) cUnStable[id] += 1;
          findDecayProducts(p,cStable,cUnStable,nstable,unstable);
        }
      }
    }

    unsigned int nCharged(const Particle& mother) const {
      unsigned int ncharged=0;
      for (const Particle & p : mother.children()) {
    	int id = p.pid();
        if (p.children().empty()) {
          if (isCharged(p)) ++ncharged;
        }
        // skip pi0 K0
        else if (id==111 || id==310 || id==311 || id==-311 || id==310) {
          continue;
        }
        else {
          ncharged+= nCharged(p);
        }
      }
      return ncharged;
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle& meson : apply<UnstableParticles>(event, "UFS").particles()) {
       	_n->fill();
        map<int,unsigned int> cStable, cUnStable;
        unsigned int nstable=0;
        findDecayProducts(meson,cStable,cUnStable,nstable,set<int>({${default_stable}}));
        bool noIntermediates = cUnStable.empty();
        ${charged}
        // inclusive modes
        ${inclusive}
        // two body
        if (meson.children().size()==2) {
          unsigned int nks0(0),nkl0(0);
          for (unsigned int ix=0; ix<2; ++ix) {
            if (meson.children()[ix].abspid()==311 ||
                meson.children()[ix].pid()==PID::K0S ||
                meson.children()[ix].pid()==PID::K0L) {
              Particle parent=meson.children()[ix];
              while (parent.children().size()==1) {
                if (parent.pid()==PID::K0S || parent.pid()==PID::K0L) break;
                parent=parent.children()[0];
              }
              if (parent.pid()==PID::K0S) ++nks0;
              else if (parent.pid()==PID::K0L) ++nkl0;
            }
          }
          unsigned int nk0=nks0+nkl0;
          ${twobody}
        }
        // exclusive modes
        ${exclusive}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      ${finalize}
      // brs
      for (unsigned int ix=1; ix<${nHist}; ++ix) {
      	if (_h_mult[ix]) scale(_h_mult[ix], 100./  *_n);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _n;
    CounterPtr _h_mult[${nHist}];
    ${histos}
    /// @}


  };


  RIVET_DECLARE_PLUGIN(${rivetName});
}""")
