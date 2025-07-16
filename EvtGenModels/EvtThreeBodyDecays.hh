#ifndef EVTTHREEBODYDECAYS_HH
#define EVTTHREEBODYDECAYS_HH

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtResonance2.hh"
#include "EvtGenBase/EvtVector4R.hh"

#include <nlohmann/json.hpp>

#include <map>
#include <memory>    // for std::unique_ptr
#include <string>
#include <tuple>
#include <vector>

using json = nlohmann::json;


class EvtThreeBodyDecays : public EvtDecayAmp {
  public:
    EvtThreeBodyDecays();
    virtual ~EvtThreeBodyDecays();

    std::string getName();
    EvtDecayBase* clone();

    void init();
    void initProbMax();
    void decay( EvtParticle* p );

  protected:
    


  private:
    nlohmann::json decayDescription;
    std::map<std::string, nlohmann::json> functions;
    nlohmann::json domains;
    nlohmann::json misc;
    nlohmann::json parameterpoints;

    std::array<double, 2> calculateFormFactors(nlohmann::json chain, std::map<std::string, nlohmann::json>& functions,
    std::array<double, 3> σs,
    int id1, int id2, int id3,
    double mParent, double m1, double m2, double m3);

};



#endif
