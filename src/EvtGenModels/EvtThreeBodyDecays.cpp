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

#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <unistd.h>    // Für getcwd()

using json = nlohmann::json;

bool deubevt = false;
bool mandeldebug = false;

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

// Initialisiere das Modell
void EvtThreeBodyDecays::initProbMax()
{
    setProbMax( 1500000 );
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
    //ms = {0.938,0.494,0.140,2.286};
    ThreeBodyMasses mssquared = {mDaug[0] * mDaug[0], mDaug[1] * mDaug[1],
                                  mDaug[2] * mDaug[2], mParent * mParent};
    //ThreeBodySpins spins = {1, 0, 0, 1};    // h0=1 bezieht sich auf den Spin des Elternteilchens
    ThreeBodySpins spins = {1,0,0,1};    // h0=1 bezieht sich auf den Spin des Elternteilchens
    ThreeBodySystem tbs  = {ms, spins};
    MandelstamTuple σstest = tbDecays.x2σs({0.3, 0.3}, ms, 1);

    // Get 4-momenta of all particles
    EvtVector4R p0 = p->getP4();
    EvtVector4R p1 = p->getDaug(0)->getP4();
    EvtVector4R p2 = p->getDaug(1)->getP4();
    EvtVector4R p3 = p->getDaug(2)->getP4();
    
    // Calculate invariant masses
    double m0 = p0.mass();
    double m1 = p1.mass();
    double m2 = p2.mass();
    double m3 = p3.mass();
    
    // Calculate Mandelstam variables (invariants)
    double s12 = (p1 + p2).mass2();
    double s23 = (p2 + p3).mass2();
    double s31 = (p3 + p1).mass2();
    MandelstamTuple σs = {s12, s23, s31};
    σs = {s23, s31, s12};
    if(mandeldebug) std::cout << "Mandelstam variables: s12 = " << s12 << ", s23 = " << s23 << ", s31 = " << s31 << std::endl;
    if(mandeldebug) std::cout << "Invariant masses: m0 = " << m0 << ", m1 = " << m1 << ", m2 = " << m2 << ", m3 = " << m3 << std::endl;
    



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

                auto Xlineshape = BreitWigner( mass, width );

                // Get k from topology
                int kint = topology[1].get<int>();    // Convert from 1-based to 0-based index

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
                if(deubevt) std::cout << weight << std::endl;
                MandelstamTuple sigma = {
                    std::pow( ( daughtermass[0] + daughtermass[1] ), 2 ),
                    std::pow( ( daughtermass[1] + daughtermass[2] ), 2 ),
                    std::pow( ( daughtermass[2] + daughtermass[0] ), 2 ) };
                std::vector<double> two_λs = { 0.5, 0.5, 0.5 };

                Tensor4Dcomp A_chain_values = tbDecays.amplitude4dcomp( *dc, σs , kint);
                if(deubevt) std::cout << "Single amp tensor:" << std::endl;
                for ( int i = 0; i < A_chain_values.size(); ++i ) {
                    for ( int j = 0; j < A_chain_values[0].size(); ++j ) {
                        for ( int k = 0; k < A_chain_values[0][0].size(); ++k ) {
                            for ( int z = 0; z < A_chain_values[0][0][0].size(); ++z ) {
                                if(deubevt) std::cout << A_chain_values[i][j][k][z] << "\t";    // Tab für schöne Ausrichtung
                            }
                        }
                    }
                    if(deubevt) std::cout << "\n";
                }

                model.add(dc, resonanceName, weight);
                //amplitude = amplitude(dc, sigma)
                //amplitude = amplitude(dc, sigma, two_λs);

                //std::array = {
                //    std::pow( ( daughtermass[0] + daughtermass[1] ), 2 ),
                //    std::pow( ( daughtermass[1] + daughtermass[2] ), 2 ),
                //    std::pow( ( daughtermass[2] + daughtermass[0] ), 2 ) };
                //A_chain_values[i] = aligned_amplitude(chain, invariants, masses); // matrices in helicities [+1/2,-1/2] x [+1/2,-1/2]
            }
        }
    }


    Tensor4Dcomp amp = model.amplitude4d(σs, 1);
    std::vector<std::vector<EvtComplex>> amplitude = {
        {EvtComplex(amp[0][0][0][0].real(), amp[0][0][0][0].imag()),
         EvtComplex(amp[0][0][0][1].real(), amp[0][0][0][1].imag())},
        {EvtComplex(amp[1][0][0][0].real(), amp[1][0][0][0].imag()),
         EvtComplex(amp[1][0][0][1].real(), amp[1][0][0][1].imag())}
    };

    if(deubevt) {
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

}



