#include "EvtGenModels/EvtThreeBodyDecays.hh"

#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtParticleFactory.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtResonance2.hh"
#include "EvtGenBase/EvtSpinType.hh"
#include "EvtGenBase/EvtdFunction.hh"

#include <nlohmann/json.hpp>
#include "ThreeBodyDecays/ThreeBodyDecays.hh"
#include "ThreeBodyDecays/ThreeBodyAmplitudeModel.hh"
#include "ThreeBodyDecays/FormFactors.hh"

#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <unistd.h>    // Für getcwd()

using json = nlohmann::json;

bool deubevt = false;
bool mandeldebug = false;
bool fractionout = true; // Für Debugging-Zwecke, um Brüche auszugeben
bool ffdebug =true; // Für Debugging-Zwecke, um Formfaktoren auszugeben
bool compjulia = true; // debug - compare with julia results


// Konstruktor
EvtThreeBodyDecays::EvtThreeBodyDecays()
{
}

// Destruktor
EvtThreeBodyDecays::~EvtThreeBodyDecays()
{
}


std::string EvtThreeBodyDecays::getName()
{
    return "THREEBODYDECAYS";
}

EvtDecayBase* EvtThreeBodyDecays::clone()
{
    return new EvtThreeBodyDecays;
}

complex evtparseComplex(const std::string& str)
{
    std::istringstream iss(str);
    double real, imag;
    
    // Try to read real part
    if (!(iss >> real)) {
        throw std::invalid_argument("Cannot parse real part of complex number");
    }
    
    // Check for + or - sign
    char plus_minus;
    iss >> plus_minus;
    
    if (plus_minus != '+' && plus_minus != '-') {
        // If no imaginary part, return just real part
        return complex(real, 0.0);
    }
    
    // Try to read imaginary part
    if (!(iss >> imag)) {
        throw std::invalid_argument("Cannot parse imaginary part of complex number");
    }
    
    // Adjust sign if necessary
    if (plus_minus == '-') {
        imag = -imag;
    }
    
    // Skip 'i' character if present
    char i;
    iss >> i;
    
    return complex(real, imag);
}

double parseHelicity(const std::string& helicityStr) {
    // Check if the string contains a fraction (has a '/')
    size_t slashPos = helicityStr.find('/');
    if (slashPos != std::string::npos) {
        // Extract numerator and denominator
        std::string numeratorStr = helicityStr.substr(0, slashPos);
        std::string denominatorStr = helicityStr.substr(slashPos + 1);
        
        // Convert to integers
        double numerator = std::stod(numeratorStr);
        double denominator = std::stod(denominatorStr);
        
        // Return the fraction
        return numerator / denominator;
    }
    
    // If not a fraction, just convert directly to double
    return std::stod(helicityStr);
}

double parseFractionlocal(const std::string& str)
{
    size_t pos = str.find('/');
    if (pos != std::string::npos) {
        // Es ist ein Bruch "a/b"
        int numerator = stoi(str.substr(0, pos));      // Zähler
        int denominator = stoi(str.substr(pos + 1));  // Nenner
        return static_cast<double>(numerator) / denominator;
    }
    // Es ist eine ganze Zahl
    return std::stod(str);
}



double BlattWeisskopf(double q, int L, double d) {
    // Calculate z² and z0²
    double z2 = (q * d) * (q * d);
    //double z02 = (q0 * d) * (q0 * d);
    
    // Return appropriate factor based on angular momentum
    switch(L) {
        case 0:
            return 1.0;  // Always 1 for  S-wave
            
        case 1:
            // P-wave: sqrt(z0²/chi_1(z0²) / z²/chi_1(z²))
            return sqrt( ( z2) / (1.0 + z2));
            
        case 2:
            // D-wave: sqrt(z0⁴/chi_2(z0²) / z⁴/chi_2(z²))
            return sqrt( (z2*z2)/ 
                         (9.0 + 3.0*z2 + z2*z2));
            
        default:
            return 1.0;  // Return 1.0 for unsupported L values
    }
}

  double breakupMomentum(double m, double m1, double m2) {
    return sqrt((m - (m1 + m2)) * (m + (m1 + m2))) * 
           sqrt((m - (m1 - m2)) * (m + (m1 - m2))) / (2*m);
}


double breakup(double M, double m1, double m2) {
    if (M < m1 + m2) return 0.0;
    double lambda = (M*M - (m1 + m2)*(m1 + m2)) * (M*M - (m1 - m2)*(m1 - m2));
    return 0.5 * std::sqrt(lambda) / M;
}


double estimateMaxProb(const std::vector<double>& formfactors, double baseProb = 1.0) {
    double factor = 1.0;
    for (double ff : formfactors) {
        // Quadratische Wirkung des Formfaktors auf die Amplitude
        factor *= (ff * ff);
    }

    // Sicherheitsfaktor, um Ausreißer abzufangen
    return baseProb * factor * 2.0; // oder * 5.0 bei mehr Absicherung
}

// Initialisiere das Modell
void EvtThreeBodyDecays::initProbMax()
{
    setProbMax( 1500000 );
    setProbMax( 100000 );
    //setProbMax( 10000000 );    // Setze die maximale Wahrscheinlichkeit für den Zerfall
    //setProbMax(10000);    // Setze die maximale Wahrscheinlichkeit für den Zerfall

    
}


// Initialisiere das Modell
void EvtThreeBodyDecays::init()
{

    // Check that there are the correct number of daughters
    checkNDaug( 3 );

    // We expect specific spin types for the parent and daughters
    checkSpinParent(
        EvtSpinType::DIRAC );    // Ändere dies entsprechend dem erwarteten Spin-Typ

    EvtSpinType::spintype d0type = EvtPDL::getSpinType( getDaug( 0 ) );
    EvtSpinType::spintype d1type = EvtPDL::getSpinType( getDaug( 1 ) );
    EvtSpinType::spintype d2type = EvtPDL::getSpinType( getDaug( 2 ) );

    // JSON-Datei laden
    //std::ifstream file("lc2ppik-lhcb-2683025_copy.json");
    //std::ifstream file( "lc2ppik-lhcb-2683025.json" );

    std::string jsonFilePath = getArgStr( 0 );

    // JSON-Datei laden
    std::ifstream file( jsonFilePath );

    if ( !file ) {
        // Zeige den aktuellen Arbeitsordner an, um Probleme mit dem Dateizugriff zu diagnostizieren
        char buf[256];
        std::cout << "Fehler beim Öffnen der Datei! Aktueller Pfad: "
                  << getcwd( buf, sizeof( buf ) ) << std::endl;
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "Konnte die JSON-Datei nicht öffnen!" << std::endl;
        return;
    }

    json decayData;
    try {
        file >> decayData;
    } catch ( const std::exception& e ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Error parsing JSON file: " << e.what() << std::endl;
        return;
    }

    // Log the entire JSON data for deubevtging
    //EvtGenReport(EVTGEN_INFO, "EvtGen") << "Loaded JSON data: " << decayData.dump(4) << std::endl;

    // Überprüfe, ob die Zerfallsbeschreibung vorhanden ist
    for ( const auto& dist : decayData["distributions"] ) {
        if ( dist["name"] == "default_model" ) {
            decayDescription = dist["decay_description"];
            break;
        }
    }

    if ( decayDescription.empty() ) {
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "Keine Zerfallsbeschreibung in der JSON-Datei gefunden!"
            << std::endl;
        return;
    }

    auto final_state = decayDescription["kinematics"]["final_state"];
    auto topology = decayDescription["reference_topology"];

    // Lade alle Funktionen aus der JSON-Datei in eine Lookup-Tabelle
    for ( const auto& func : decayData["functions"] ) {
        functions[func["name"]] = func;
    }

    // Überprüfe, ob die Funktionen geladen wurden
    if ( functions.empty() ) {
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "Keine Funktionsdefinitionen in der JSON-Datei gefunden!"
            << std::endl;
        return;
    }

    // Log the loaded functions for deubevtging
    EvtGenReport( EVTGEN_INFO, "EvtGen" )
        << "Loaded functions from JSON file:" << std::endl;
    for ( const auto& func : functions ) {
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "Function name: " << func.first << std::endl;
    }

    double mDaug[3] = {EvtPDL::getMeanMass( EvtPDL::getId( "K-" ) ),
        EvtPDL::getMeanMass( EvtPDL::getId( "pi+" ) ),
        EvtPDL::getMeanMass( EvtPDL::getId( "p+" ) )};
    double mParent = EvtPDL::getMass( EvtPDL::getId( "Lambda_c+" ) );

    ThreeBodyMasses ms = {mDaug[0], mDaug[1], mDaug[2], mParent};
    ThreeBodyMasses mssquared = {mDaug[0] * mDaug[0], mDaug[1] * mDaug[1],
                                  mDaug[2] * mDaug[2], mParent * mParent};
    ThreeBodySpins spins = {1, 0, 0, 1};    // h0=1 bezieht sich auf den Spin des Elternteilchens
    ThreeBodySystem tbs  = {ms, spins};

    //ThreeBodyParities Conserving = {'-', '+', '-', '+'};
    //ThreeBodyParities Violating = {'-', '+', '-', '-'};

    // Erstellen eines BreitWigner mit Masse 1.0 und Breite 0.1
    //BreitWigner bw(1.0, 0.1);

    if(deubevt) EvtGenReport( EVTGEN_INFO, "EvtGen" ) << "INFO" << std::endl;





}

std::vector<std::pair<std::string, double>> weighttuple; 
std::array<double, 9> intensities;


double intensity = 0;
double newtotalintensity = 0;
//std::vector<double> compintens;
std::array<double, 20> compintens2;
int num = 0;

std::array<std::vector<double>,9> intensityfractions;
std::array<std::vector<double>,20> allintensityfractions;


void EvtThreeBodyDecays::decay( EvtParticle* p )
{
    if ( !p ) {
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "Null Pointer auf Partikel!" << std::endl;
        return;
    }

    // Generate a phase space decay to initialize daughter particles
    p->initializePhaseSpace(getNDaug(), getDaugs());
    
    // Get the daughter particles
    EvtParticle* proton = p->getDaug(0);
    EvtParticle* pion = p->getDaug(1);
    EvtParticle* kaon = p->getDaug(2);
    
    // Create random Mandelstam variables for the decay
    ThreeBodyDecays tbDecays;
    double mDaug[3] = {
        EvtPDL::getMeanMass( EvtPDL::getId( "p+" ) ),
        EvtPDL::getMeanMass( EvtPDL::getId( "pi+" ) ),
        EvtPDL::getMeanMass( EvtPDL::getId( "K-" ) ),
    };
    double mParent = EvtPDL::getMass( EvtPDL::getId( "Lambda_c+" ) );

    ThreeBodyMasses ms = {mDaug[0], mDaug[1], mDaug[2], mParent};
    ms =  {0.938272046,
          0.13957018,
          0.493677,
          2.28646}; // masses in GeV/c^2
    
   
    //ms = {0.938,0.494,0.140,2.286};
    ThreeBodyMasses mssquared = {mDaug[0] * mDaug[0], mDaug[1] * mDaug[1],
                                  mDaug[2] * mDaug[2], mParent * mParent};
    //ThreeBodySpins spins = {1, 0, 0, 1};    // h0=1 bezieht sich auf den Spin des Elternteilchens
    ThreeBodySpins spins = {1,0,0,1};    // h0=1 bezieht sich auf den Spin des Elternteilchens
    ThreeBodySystem tbs  = {ms, spins};
    MandelstamTuple σstest = tbDecays.x2σs({0.3, 0.3}, ms, 1);

    // Get 4-momenta of all particles
    EvtVector4R pLc = p->getP4();
    EvtVector4R pp = p->getDaug(0)->getP4();
    EvtVector4R ppi = p->getDaug(1)->getP4();
    EvtVector4R pk = p->getDaug(2)->getP4();
    
    // Calculate invariant masses
    double mLc = pLc.mass();
    double mp = pp.mass();
    double mpi = ppi.mass();
    double mk = pk.mass();

    double s12 = (pp + ppi).mass2();
    double s23 = (ppi + pk).mass2();
    double s31 = (pk + pp).mass2();

    std::vector<double> formfactors;

    MandelstamTuple σs = {s23, s31, s12};

    if(compjulia ) {
        σs = {1.5,3.2,1.6714505793792584};

    }
    if(mandeldebug) std::cout << "Mandelstam variables: s12 = " << s12 << ", s23 = " << s23 << ", s31 = " << s31 << std::endl;
    if(mandeldebug) std::cout << "Invariant masses: mLc = " << mLc << ", mp = " << mp << ", mpi = " << mpi << ", mk = " << mk << std::endl;
    
    //σs = {1,3,2};
    //std::vector<std::pair<std::string, double>> weighttuple; 

    auto chains = decayDescription["chains"];
    auto final_state = decayDescription["kinematics"]["final_state"];
    ThreeBodyAmplitudeModel model;
    std::vector<std::string> paritylist = {"+", "-", "+", "-"};

    for ( const auto& chain : chains ) {
        std::string resonanceName = chain["name"];
        std::string paramType = chain["propagators"][0]["parametrization"];
        if ( functions.find( paramType ) != functions.end() ) {
            const auto& func = functions[paramType];
            auto topology = chain["topology"];

            if ( paramType.find( "_BW" ) != std::string::npos ) {
                double mass = func["mass"];
                double width = func["width"];
                double l = func["l"];
                double mb = func["mb"];
                double ma = func["ma"];
                double d = func["d"];
                std::string spin = chain["propagators"][0]["spin"];
                if(deubevt) EvtGenReport( EVTGEN_INFO, "EvtGen" ) << "INFO" << std::endl;

                complex weight = evtparseComplex(
                    chain["weight"]);
                double weightSquared = abs( weight ) * abs( weight );

                int id1 = topology[0][0];
                int id2 = topology[0][1];
                int id3 = topology[1];
                id1 = id1 - 1;
                id2 = id2 - 1;
                id3 = id3 - 1;

                EvtVector4R moms1 = p->getDaug( id1 )->getP4();
                EvtVector4R moms2 = p->getDaug( id2 )->getP4();
                EvtVector4R moms3 = p->getDaug( id3 )->getP4();
                /*double m1 = moms1.mass();
                double m2 = moms2.mass();
                double m3 = moms3.mass();
                */
                double m1 = p->getDaug(id1)->mass();
                double m2 = p->getDaug(id2)->mass();
                double m3 = p->getDaug(id3)->mass();

                double daughtermass[3] = { p->getDaug( id1 )->mass(),
                                           p->getDaug( id2 )->mass(),
                                           p->getDaug( id3 )->mass() };
                
                

                // Create ThreeBodySystem
                auto masses_tbs = ThreeBodyMasses{ p->getDaug( 0 )->mass(),
                                                   p->getDaug( 1 )->mass(),
                                                   p->getDaug( 2 )->mass(),
                                                   p->mass() };

                auto spins_tbs = ThreeBodySpins{
                    2 * p->getDaug( 0 )->getSpinType(),
                    2 * p->getDaug( 1 )->getSpinType(),
                    2 * p->getDaug( 2 )->getSpinType(), 2 * p->getSpinType() };

                //auto tbs = ThreeBodySystem( masses_tbs, spins_tbs );

                // Get spin-parity and parities from JSON
                if(deubevt) EvtGenReport( EVTGEN_INFO, "EvtGen" ) << "INFO" << std::endl;
                std::string spin_str = chain["propagators"][0]["spin"];
                if(deubevt) std::cout << spin_str << std::endl;
                std::string parity_str = chain["vertices"][1]["parity_factor"];
                int vertixy_ind = 0;
                std::string jp = spin_str + parity_str;
                if(deubevt) std::cout << jp << std::endl;


        


                std::vector<std::array<double, 3>> helicities;
                
                
                for ( const auto& vertice : chain["vertices"] )
                {
                    double helfactor1 = parseHelicity( vertice["helicities"][0]);
                    double helfactor2 = parseHelicity( vertice["helicities"][1]);
                    double heltype;
                    if ( vertice["type"] == "parity" ) {
                        std::string formfactor = vertice["formfactor"];
                        heltype = 1;
                    }
                    if ( vertice["type"] == "helicity" ) {
                        std::string formfactor = vertice["formfactor"];
                        heltype = 0;
                    }
                    helicities.push_back({ helfactor1, helfactor2, heltype });
                }

                /*
                "vertices": [
              {
                "type": "helicity",
                "helicities": ["-1", "-1/2"],
                "node": [[2, 3], 1],
                "formfactor": ""
              },
              */
                
                int helfactor1 = parseHelicity( chain["vertices"][0]["helicities"][0])*2;
                int helfactor2 = parseHelicity( chain["vertices"][0]["helicities"][1])*2;
                std::array<int, 2> hel = { helfactor1, helfactor2 };
                if(deubevt) {
                    std::cout << "Helicities: " << helfactor1 << " " << helfactor2 << std::endl;
                }

                int parfactor1 = parseHelicity( chain["vertices"][1]["helicities"][0])*2;
                int parfactor2 = parseHelicity( chain["vertices"][1]["helicities"][1])*2;
                std::array<int, 2> par = { parfactor1, parfactor2 };
                std::string par_type = chain["vertices"][1]["parity_factor"];
                bool par_sign = false;
                if(par_type == "+") {
                    par_sign = true;
                } else if(par_type == "-") {
                    par_sign = false;
                } else {
                    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                        << "Unbekannter Paritätsfaktor: " << par_type << std::endl;
                }

                if(deubevt) {
                    std::cout << "Parity factors: " << parfactor1 << " " << parfactor2 << std::endl;
                    std::cout << "Parity sign: " << par_sign << std::endl;
                }
                //Ps = {0,0,0,0};


                ///////  Blatt Weisskopf Faktoren  ///////
                // Berechne Zerfallsimpulse für Blatt-Weisskopf-Formfaktoren
                double qold = 0.0, q0old = 0.0;
                double formFactor1 = 1.0;
                double formFactor2 = 1.0;

                // Bestimme die aktuelle invariante Masse basierend auf der Topologie
                int k = topology[1].get<int>();
                double invariantMassSq = 0.0;
                
                
                

                auto Xlineshapeold = BreitWigner(mass, width);
                // Log the calculated q and q0 for deubevtging

                // Überprüfe auf vorhandene Formfaktoren in den Vertices
                for (const auto& vertex : chain["vertices"]) {
                    if (vertex.contains("formfactor") && !vertex["formfactor"].get<std::string>().empty()) {
                        std::string formfactorName = vertex["formfactor"];
                        if (functions.find(formfactorName) != functions.end()) {
                            const auto& ff = functions[formfactorName];
                            double radius = ff["radius"];
                            int L = ff["l"];
                            
                            
                            auto node = vertex["node"];
                            if (node.is_array() && node[0].is_array()) {


               
                                double ssub = σs[id3];
                                double msub = sqrt(ssub);

                                double q = breakup(mLc, msub, m3);
                                //double q0 = breakup(m0, mass, m3);

                                // Log the calculated q and q0 for deubevtging
                                // Boost particle 3 into resonance rest frame (p1 + p2)
                                /*
                                EvtVector4R resP4 = moms1 + moms2;
                                EvtVector4R p3_boosted = moms3;
                                p3_boosted.applyBoostTo(resP4);

                                q = p3_boosted.d3mag();             // boosted momentum*/
                               

                                
                                formFactor1 *= BlattWeisskopf(q, L, radius);
                                formfactors.push_back(formFactor1);

                                if(ffdebug) {
                                std::cout << resonanceName << " Formfactor: " << formfactorName << std::endl;
                                std::cout << "q: " << q << ", q0: " << " l " << L << " r " << radius << std::endl;
                                std::cout << "Formfactor1: " << formFactor1 << std::endl;
                                std::cout << " msub: " << msub << " idx " << id1 << " " << id2 << " " << id3 << "m" << m1 << " " << m2 << " " << m3 << std::endl;
                                }
                        

                            } else {
                            
                                double ssub = σs[id3];
                                double msub = sqrt(ssub);

                                // Impuls bei gegebener (laufender) Masse
                                double q = breakup(msub, m1, m2);

                                // Impuls bei Resonanzmasse
                                //double q0 = breakup(mass, m1, m2);

                                formFactor2 *= BlattWeisskopf(q, L, radius);
                                formfactors.push_back(formFactor2);

                                if(ffdebug) {
                                std::cout << resonanceName << " Formfactor: " << formfactorName << std::endl;
                                std::cout << "q: " << q << ", q0: "  << " l " << L << " r " << radius << std::endl;
                                std::cout << "Formfactor2: " << formFactor2 << std::endl;
                                std::cout << " msub: " << msub << " idx " << id1 << " " << id2 << " " << id3 << "m" << m1 << " " << m2 << " " << m3 << std::endl;

                                }
                            }

                            
                           

                            //Xlineshape2 = Xlineshape2 * complex(formFactor,formFactor); // Modifiziere den Breit-Wigner mit dem Formfaktor
                        }
                    }
                }
                //std::cout << resonanceName << "Formfaktor1: " << formFactor1 << std::endl;
                //std::cout << resonanceName << "Formfaktor2: " << formFactor2 << std::endl;
                
                int kint = topology[1].get<int>(); 

                // Modifiziere das Breit-Wigner mit dem Formfaktor
                auto originalBreitWigner = BreitWigner(mass, width);
                auto Xlineshape = [originalBreitWigner, formFactor1,formFactor2](double s) -> complex {
                    return originalBreitWigner(s) * formFactor1 * formFactor2;
                };
                if(compjulia) std::cout << resonanceName << "XLineshape: " << Xlineshape(σs[kint-1]) << " with s=" << σs[kint-1]<< std::endl;
                if(compjulia) std::cout << resonanceName << "Xlineshapeold: " << Xlineshapeold(σs[kint-1]) << " " << mass << " " << width << std::endl;
                if(compjulia) std::cout << resonanceName << "FF " << formFactor1 << " " << formFactor2 << std::endl;
                // Get k from topology
                   

                        if(deubevt) EvtGenReport( EVTGEN_INFO, "EvtGen" ) << "INFO" << std::endl;
                
                if(deubevt) std::cout << resonanceName << std::endl;
                if(deubevt) std::cout << "k: " << kint << std::endl;
                if(deubevt) std::cout << mass << " " << width << std::endl;
                if(deubevt) std::cout << "masses: " << ms[0] << " " << ms[1] << " "
                          << ms[2] << " " << ms[3] << std::endl;
                if(deubevt) std::cout << "spins: " << spins[0] << " " << spins[1] << " "
                          << spins[2] << spins[3] << std::endl;
                
                ThreeBodyParities Ps = {'+', '-', '-', '+'};

                // Create DecayChain
                auto dc = createDecayChainCoupling( kint, Xlineshape, jp, tbs, RecouplingType::NoRecoupling, hel, false, RecouplingType::ParityRecoupling, par, par_sign );
                if(compjulia) std::cout << resonanceName << " k " << kint << " Lineshape " << Xlineshape(0.0) << jp << " helicity " << hel[0] << hel[1] << " parity " << par[0] << par[1] << par_sign << std::endl;
                model.add(dc, resonanceName, weight);

                if(deubevt or compjulia) std::cout << weight << std::endl;
                std::vector<double> two_λs = { 0.5, 0.5, 0.5 };

                int weightint = 1;

                Tensor4Dcomp A_chain_values = tbDecays.amplitude4dcomp( *dc, σs , 1);
                if(compjulia or deubevt){ std::cout << "Single amp tensor:" << std::endl;
                for ( int i = 0; i < A_chain_values.size(); ++i ) {
                    for ( int j = 0; j < A_chain_values[0].size(); ++j ) {
                        for ( int k = 0; k < A_chain_values[0][0].size(); ++k ) {
                            for ( int z = 0; z < A_chain_values[0][0][0].size(); ++z ) {
                                std::cout << A_chain_values[i][j][k][z] << "\t";    // Tab für schöne Ausrichtung
                            }
                        }
                    }
                    std::cout << "\n";
                }
                for ( int i = 0; i < A_chain_values.size(); ++i ) {
                        for ( int j = 0; j < A_chain_values[0].size(); ++j ) {
                            for ( int k = 0; k < A_chain_values[0][0].size(); ++k ) {
                                for ( int z = 0; z < A_chain_values[0][0][0].size(); ++z ) {
                                    std::cout << A_chain_values[i][j][k][z]*weight << "\t";    // Tab für schöne Ausrichtung
                                }
                            }
                        }
                        std::cout << "\n";
                    }}

                //weight = complex(1.,1.);
                //std::cout << resonanceName << weight << std::endl;

                double realintensity = 0;
                for ( int i = 0; i < A_chain_values.size(); ++i ) {
                    for ( int j = 0; j < A_chain_values[0].size(); ++j ) {
                        for ( int k = 0; k < A_chain_values[0][0].size(); ++k ) {
                            for ( int z = 0; z < A_chain_values[0][0][0].size(); ++z ) {
                                //realintensity += A_chain_values[i][j][k][z].real() * A_chain_values[i][j][k][z].real() * weight.real() * weight.real() +
                                //                + A_chain_values[i][j][k][z].imag() * A_chain_values[i][j][k][z].imag() * weight.imag() * weight.imag();
                                double amp_intensity = A_chain_values[i][j][k][z].real() * A_chain_values[i][j][k][z].real() * weight.real() * weight.real() +
                                                A_chain_values[i][j][k][z].imag() * A_chain_values[i][j][k][z].imag() * weight.imag() * weight.imag();
                                realintensity += amp_intensity; // Add intensity for this amplitude
                                //std::cout << resonanceName << " Amplitude[" << i << "][" << j << "][" << k << "][" << z << "] = " 
                                //          << A_chain_values[i][j][k][z] << weight << amp_intensity << std::endl;
                                if(deubevt) std::cout << "Amplitude[" << i << "][" << j << "][" << k << "][" << z << "] = " << realintensity << std::endl;
                            }
                        }
                    }
                }


                bool frwithweight = true;
                if (!frwithweight) {
                    realintensity = tbDecays.intensity( *dc, σs, weightint ); // Calculate intensity without weight
                }
                
                intensity += realintensity; // Add real intensity to total intensity
                //std::cout << resonanceName << "Intensity: " << tbDecays.intensity( *dc, σs, kint) << std::endl;
                //std::cout << resonanceName << "Weight Intensity: " << tbDecays.intensity( *dc, σs, kint , weight) << std::endl;

                

                // Find existing entry for this resonance or add a new one
                bool found = false;
                for (size_t i = 0; i < weighttuple.size(); i++) {
                    if (resonanceName == weighttuple[i].first) {
                        // Update existing entry
                    
                        //weighttuple[i].second += tbDecays.intensity(*dc, σs, kint, weight);
                        weighttuple[i].second += realintensity; // Add real intensity to existing entry
                    
                        //weighttuple[i].second += tbDecays.intensity(*dc, σs, 1);

                        found = true;
                        //std::cout << "Updated entry for " << resonanceName << std::endl;
                        break; // Found the matching entry, no need to continue the loop
                    }
                }

                

                // If no entry exists for this resonance, add a new one
                if (!found) {
                    //weighttuple.push_back(std::make_pair(resonanceName, tbDecays.intensity(*dc, σs, kint)));
                    weighttuple.push_back(std::make_pair(resonanceName, realintensity));

                    //std::cout << "Added new entry for " << resonanceName << std::endl;
                }

                
                
                //weighttuple.push_back(std::make_pair(resonanceName, tbDecays.intensity(*dc, σs, kint)));//amplitude = amplitude(dc, sigma)
                //amplitude = amplitude(dc, sigma, two_λs);

                //std::array = {
                //    std::pow( ( daughtermass[0] + daughtermass[1] ), 2 ),
                //    std::pow( ( daughtermass[1] + daughtermass[2] ), 2 ),
                //    std::pow( ( daughtermass[2] + daughtermass[0] ), 2 ) };
                //A_chain_values[i] = aligned_amplitude(chain, invariants, masses); // matrices in helicities [+1/2,-1/2] x [+1/2,-1/2]
            }
        }
    }
    num++;
    

    double maxProbEstimate = estimateMaxProb(formfactors);
    //std::cout << "Empfohlene maxProb für Kette " << ": " << maxProbEstimate << std::endl;


   


    Tensor4Dcomp amp = model.amplitude4d(σs, 1);
    std::vector<std::vector<EvtComplex>> amplitude = {
        {EvtComplex(amp[0][0][0][0].real(), amp[0][0][0][0].imag()),
         EvtComplex(amp[0][0][0][1].real(), amp[0][0][0][1].imag())},
        {EvtComplex(amp[1][0][0][0].real(), amp[1][0][0][0].imag()),
         EvtComplex(amp[1][0][0][1].real(), amp[1][0][0][1].imag())}
    };

    bool ampout = false;
    if(deubevt or ampout or compjulia) {
    // Print the amplitude tensor
    EvtGenReport( EVTGEN_INFO, "EvtGen" ) << "RESULT" << std::endl;
    std::cout << "Amplitude tensor:" << std::endl;
    for ( int i = 0; i < amp.size(); ++i ) {
        for ( int j = 0; j < amp[0].size(); ++j ) {
            for ( int k = 0; k < amp[0][0].size(); ++k ) {
                for ( int z = 0; z < amp[0][0][0].size(); ++z ) {
                    std::cout << i << j << k << z<< amp[i][j][k][z] << "\t";    // Tab für schöne Ausrichtung                }
            }
        }
        std::cout << "\n";
    }
    }}
    
    vertex( 0, 0, amplitude[0][0] );
    vertex( 0, 1, amplitude[0][1] );
    vertex( 1, 0, amplitude[1][0] );
    vertex( 1, 1, amplitude[1][1] );
    if(deubevt) EvtGenReport( EVTGEN_EMERGENCY, "EvtGen" ) << "break" << std::endl;
    double modelintensity = model.intensity(σs, 1);
    double modintens = 0;
    const auto compintens = model.component_intensities(σs, 1);
    int j = 0;
    for (int i = 0; i < 8; i++) {
        intensities[i] += compintens[j] + compintens[j+1]; // Initialize intensities to zero
        modintens += compintens[j] + compintens[j+1];
        //std::cout << "Component " << i << ": " << intensities[i] << " " << intensities[i]/newtotalintensity << std::endl;
        intensityfractions[i].push_back((compintens[j] + compintens[j+1])/modelintensity); // Store the component intensities in the vector
        j+= 2; // Increment j by 2 to access the next pair of components

    }
    intensities[8] += compintens[19] + compintens[18] + compintens[17] + compintens[16] ;
    intensityfractions[8].push_back((compintens[19] + compintens[18] + compintens[17] + compintens[16])/modelintensity);
    //newtotalintensity += compintens[19] + compintens[18] + compintens[17] + compintens[16] ; // Add the last component intensity
    modintens += compintens[19] + compintens[18] + compintens[17] + compintens[16] ;

    //newtotalintensity += modintens;
    newtotalintensity += modelintensity; // Add the component intensities to the new total intensity

    for (size_t i=0; i < compintens.size(); i++) {
        compintens2[i] += (compintens[i]);
    }

    for (size_t i=0; i < compintens.size(); i++) {
        allintensityfractions[i].push_back(compintens[i]/modintens); // Store the component intensities in the vector
    } 


    if(num%1000 == 0) {
        if(fractionout) std::cout << "Intensity: " << intensity << std::endl;
        if(fractionout) std::cout << "Fraction:" << std::endl;
        for (const auto& wtup : weighttuple) {
            if(fractionout) std::cout << wtup.first << " " << wtup.second/intensity << std::endl;
        }
    }
    

    if(num%1000 == 0) {
    if(fractionout) 
    {std::cout << "Component intensities:" << std::endl;
    for (size_t i = 0; i < intensities.size(); i++) {
        // Calculate running mean
        double running_mean = intensities[i] / newtotalintensity;
        
        // Get individual fraction values from intensityfractions vector
        if (i < intensityfractions.size() && !intensityfractions[i].empty()) {
            // Calculate mean of stored fractions
            double mean = std::reduce(intensityfractions[i].begin(), intensityfractions[i].end(), 0.0) / intensityfractions[i].size();
            
            // Calculate variance
            double variance = 0.0;
            for (const auto& value : intensityfractions[i]) {
                variance += (value - mean) * (value - mean);
            }
            variance /= intensityfractions[i].size();
            double stddev = std::sqrt(variance);
            
            std::cout << weighttuple[i].first << " Intensity " << i << ": " 
                      << running_mean << " (Mean: " << mean << " ± " << stddev 
                      << ", N=" << intensityfractions[i].size() << ")" << std::endl;
        } else {
            std::cout << weighttuple[i].first << " Intensity " << i << ": " 
                      << running_mean << " (No variance data)" << std::endl;
        }
    }}
    }

    //std::cout << num << allintensityfractions.size() << " " << intensityfractions.size() << std::endl;
    if(num%1000 == 0) {
    if(fractionout) {
        for (int i = 0; i < allintensityfractions.size()-3; i+=2){
            auto const count = static_cast<float>(allintensityfractions[i].size());
            double mean1 =  std::reduce(allintensityfractions[i].begin(), allintensityfractions[i].end(), 0.0) / allintensityfractions[i].size();
            double mean2 =  std::reduce(allintensityfractions[i+1].begin(), allintensityfractions[i+1].end(), 0.0) / allintensityfractions[i+1].size();
            double mean3 =  std::reduce(allintensityfractions[i+2].begin(), allintensityfractions[i+2].end(), 0.0) / allintensityfractions[i+2].size();
            double mean4 =  std::reduce(allintensityfractions[i+3].begin(), allintensityfractions[i+3].end(), 0.0) / allintensityfractions[i+3].size();
            // Berechnung der Standardabweichung
            double variance1 = 0.0;
            double variance2 = 0.0;
            for (const auto& value : allintensityfractions[i]) {
                variance1 += (value - mean1) * (value - mean1);
            }
            for (const auto& value : allintensityfractions[i+1]) {
                variance2 += (value - mean2) * (value - mean2);
            }
            
            variance1 /= allintensityfractions[i].size();
            variance2 /= allintensityfractions[i+1].size();
            double stddev1 = std::sqrt(variance1);
            double stddev2 = std::sqrt(variance2);
            
            if(i > allintensityfractions.size()-3) {
                std::cout << i << " " << mean1 + mean2 + mean3 + mean4 << " ± " << stddev1 + stddev2 << std::endl;
            } else {
                std::cout << i << " " << mean1 + mean2 << " ± " << stddev1 + stddev2 << std::endl;
            }

            //std::cout << i << " " << mean1 << " ± " << stddev1 << std::endl;
            
        }
        std::cout << "Total intensity: " << newtotalintensity/num << std::endl;

    }
    }

    if(compjulia){
        std::cout << "Total intensity " << modelintensity << std::endl;

        for (size_t i = 0; i < intensities.size(); i++) {
        // Calculate running mean
        double running_mean = intensities[i] / newtotalintensity;
        
        // Get individual fraction values from intensityfractions vector
        if (i < intensityfractions.size() && !intensityfractions[i].empty()) {
            // Calculate mean of stored fractions
            double mean = std::reduce(intensityfractions[i].begin(), intensityfractions[i].end(), 0.0) / intensityfractions[i].size();
            
            // Calculate variance
            double variance = 0.0;
            for (const auto& value : intensityfractions[i]) {
                variance += (value - mean) * (value - mean);
            }
            variance /= intensityfractions[i].size();
            double stddev = std::sqrt(variance);
            
            std::cout << weighttuple[i].first << " Intensity " << i << ": " 
                      << running_mean << " (Mean: " << mean << " ± " << stddev 
                      << ", N=" << intensityfractions[i].size() << ")" << std::endl;
        } else {
            std::cout << weighttuple[i].first << " Intensity " << i << ": " 
                      << running_mean << " (No variance data)" << std::endl;
        }
    }
    }
    /*
    if(fractionout) 
    {std::cout << "Component intensities2:" << std::endl;
    for (size_t i = 0; i < compintens2.size(); i++) {
        std::cout << "Intensity " << i << ": " << compintens2[i]/newtotalintensity << std::endl;
    }}*/
    


}



