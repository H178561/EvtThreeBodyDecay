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
#include "ThreeBodyDecays/include/ThreeBodyDecays.hh"
#include "ThreeBodyDecays/include/ThreeBodyAmplitudeModel.hh"
#include "ThreeBodyDecays/include/Lineshapes.hh"
#include "ThreeBodyDecays/include/FormFactors.hh"
//#include <muParserX/mpParser.h>

#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <unistd.h>    // Für getcwd()
#include <array>

using json = nlohmann::json;
using namespace FormFactors;



bool deubevt = false;
bool mandeldebug = false;
bool fractionout = false; // Für Debugging-Zwecke, um Brüche auszugeben
bool ffdebug =false; // Für Debugging-Zwecke, um Formfaktoren auszugeben
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

//   Bugg BW function    1/(0.824^2 - σ - i * 0.824 * (σ - 0.23397706275638377) / (0.824^2 - 0.23397706275638377) * 0.478 * exp(-0.941060 * σ))


std::function<complex(double)> make_Bugg_BW_K700() {
    // Hardcoded parameters
    const double mass = 0.824;
    const double mass_sq = mass * mass;
    const double threshold = 0.23397706275638377;
    const double g = 0.478;
    const double b = 0.941060;
    
    return [=](double s) -> complex {
        double denominator_real = mass_sq - s;
        double width_term = mass * (s - threshold) / (mass_sq - threshold) * g * std::exp(-b * s);
        complex denominator(denominator_real, -width_term);
        return 1.0 / denominator;
    };
}

//"expression": "1/(0.845^2 - σ - i * 0.845 * (σ - 0.2631968154742324) / (0.845^2 - 0.2631968154742324) * 0.468 * exp(-0.772134 * σ))",


std::function<complex(double)> make_Bugg_BW_K700_xic() {
    // Hardcoded parameters
    const double mass = 0.845;
    const double mass_sq = mass * mass;
    const double threshold = 0.2631968154742324;
    const double g = 0.468;
    const double b = 0.772134;

    return [=](double s) -> complex {
        double denominator_real = mass_sq - s;
        double width_term = mass * (s - threshold) / (mass_sq - threshold) * g * std::exp(-b * s);
        complex denominator(denominator_real, -width_term);
        return 1.0 / denominator;
    };
}


// Bugg BW function für K1430
std::function<complex(double)> make_Bugg_BW_K1430() {
    // Hardcoded parameters
    const double mass = 1.375;
    const double mass_sq = mass * mass;
    const double threshold = 0.23397706275638377;
    const double g = 0.190;
    const double b = 0.020981;
    
    return [=](double s) -> complex {
        double denominator_real = mass_sq - s;
        double width_term = mass * (s - threshold) / (mass_sq - threshold) * g * std::exp(-b * s);
        complex denominator(denominator_real, -width_term);
        return 1.0 / denominator;
    };
}

//            "expression": "1/(1.425^2 - σ - i * 1.425 * (σ - 0.2631968154742324) / (1.425^2 - 0.2631968154742324) * 0.27 * exp(-0.051629 * σ))",

std::function<complex(double)> make_Bugg_BW_K1430_xic() {
    // Hardcoded parameters
    const double mass = 1.425;
    const double mass_sq = mass * mass;
    const double threshold = 0.2631968154742324;
    const double g = 0.27;
    const double b = 0.051629;

    return [=](double s) -> complex {
        double denominator_real = mass_sq - s;
        double width_term = mass * (s - threshold) / (mass_sq - threshold) * g * std::exp(-b * s);
        complex denominator(denominator_real, -width_term);
        return 1.0 / denominator;
    };
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
    auto initial_state = decayDescription["kinematics"]["initial_state"];
    auto topology = decayDescription["reference_topology"];
    misc = decayData["misc"];

    std::cout << misc << std::endl;
    domains = decayData["domains"];
    parameterpoints = decayData["parameter_points"];

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
    



    //ThreeBodyParities Conserving = {'-', '+', '-', '+'};
    //ThreeBodyParities Violating = {'-', '+', '-', '-'};

    // Erstellen eines BreitWigner mit Masse 1.0 und Breite 0.1
    //BreitWigner bw(1.0, 0.1);

    if(deubevt) EvtGenReport( EVTGEN_INFO, "EvtGen" ) << "INFO" << std::endl;

    
    createDecayModel( );


}



void EvtThreeBodyDecays::createDecayModel( ){
    

     // Get the daughter particles

    
    // Create random Mandelstam variables for the decay
    ThreeBodyDecays tbDecays;
  

    auto final_state = decayDescription["kinematics"]["final_state"];
    auto initial_state = decayDescription["kinematics"]["initial_state"];


    double mParent = initial_state["mass"];
    double mDaug1 = final_state[0]["mass"];
    double mDaug2 = final_state[1]["mass"];
    double mDaug3 = final_state[2]["mass"];

    ThreeBodyMasses ms = {mDaug1, mDaug2, mDaug3, mParent};
    
    // get spin from json spin entry
    int spin = parseFractionlocal( initial_state["spin"] ) * 2; // Convert to doubled spin value
    int spindaug1 = parseFractionlocal( final_state[0]["spin"] ) * 2; // Convert to doubled spin value
    int spindaug2 = parseFractionlocal( final_state[1]["spin"] ) * 2; // Convert to doubled spin value
    int spindaug3 = parseFractionlocal( final_state[2]["spin"] ) * 2; // Convert to doubled spin value


    ThreeBodySpins spins = {spindaug1, spindaug2, spindaug3, spin};    // h0=1 bezieht sich auf den Spin des Elternteilchens
    ThreeBodySystem tbs  = {ms, spins};
    
    auto chains = decayDescription["chains"];
    ThreeBodyAmplitudeModel model;

    for ( const auto& chain : chains ) {
        std::string resonanceName = chain["name"];
        std::string paramType = chain["propagators"][0]["parametrization"];

        const auto& func = functions[paramType];
            auto topology = chain["topology"];
            auto LineshapeName = func["name"];

            
            int kint = topology[1].get<int>(); 

            int id1 = topology[0][0];
            int id2 = topology[0][1];
            int id3 = topology[1];
            id1 = id1 - 1;
            id2 = id2 - 1;
            id3 = id3 - 1;

            double m1 = final_state[id1]["mass"];
            double m2 = final_state[id2]["mass"];
            double m3 = final_state[id3]["mass"];


            // Get spin-parity and parities from JSON
            std::string spin_str = chain["propagators"][0]["spin"];
            if(deubevt) std::cout << resonanceName << " Spin " << spin_str << std::endl;
            std::string parity_str = "+";
            int vertixy_ind = 0;
            std::string jp = spin_str + parity_str;
            if(deubevt) std::cout << resonanceName << " Parity " << jp << std::endl;
            complex weight = evtparseComplex(chain["weight"]);
            // Define RecouplingType enum if not already defined
            
            struct HelicityConfig {
                RecouplingType recouplingType;
                std::array<int, 2> helicityFactors;
                bool paritySign;
            };

            std::vector<HelicityConfig> helicities;

                
                
            for (const auto& vertice : chain["vertices"]) {
                HelicityConfig config;
                
                if (vertice["type"] == "parity") {
                    config.recouplingType = RecouplingType::ParityRecoupling;
                    config.helicityFactors = {parseHelicity(vertice["helicities"][0]) * 2, 
                                            parseHelicity(vertice["helicities"][1]) * 2};
                    config.paritySign = (chain["vertices"][1]["parity_factor"] == "+");
                }
                else if (vertice["type"] == "helicity") {
                    config.recouplingType = RecouplingType::NoRecoupling;
                    config.helicityFactors = {parseHelicity(vertice["helicities"][0]) * 2, 
                                            parseHelicity(vertice["helicities"][1]) * 2};
                    config.paritySign = false;
                }
                else if (vertice["type"] == "ls") {
                    config.recouplingType = RecouplingType::LSRecoupling;
                    config.helicityFactors = {parseHelicity(vertice["l"]) * 2, 
                                            parseHelicity(vertice["s"]) * 2};
                    config.paritySign = false;
                }
                
                helicities.push_back(config);
            }

            // Debug-Ausgabe:
            std::cout << "Helicities: " << std::endl;
            for (const auto& hel : helicities) {
                std::cout << "RecouplingType: ";
                switch (hel.recouplingType) {
                    case RecouplingType::ParityRecoupling:
                        std::cout << "ParityRecoupling" << std::endl;
                        break;
                    case RecouplingType::NoRecoupling:
                        std::cout << "NoRecoupling" << std::endl;
                        break;
                    case RecouplingType::LSRecoupling:
                        std::cout << "LSRecoupling" << std::endl;
                        break;
                }
                std::cout << "Helicity Factors: " << hel.helicityFactors[0] << ", " << hel.helicityFactors[1] << std::endl;
                std::cout << "Parity Sign: " << (hel.paritySign ? "+" : "-") << std::endl;
            }

            std::function<complex(double)> Xlineshape;
            bool lineshapeInitialized = false;

            

            if ( func["type"] == "BreitWigner") {

                double mass = func["mass"];
                double width = func["width"];
                double l = func["l"];
                double mb = func["mb"];
                double ma = func["ma"];
                double d = func["d"];
                

                auto formFactorFunc = createFormFactorFunction(chain, functions, mParent, m1, m2, m3);

                
                int kint = topology[1].get<int>(); 

                // Modifiziere das Breit-Wigner mit dem Formfaktor
                auto originalBreitWigner = make_multichannel_bw_single(mass, width, ma, mb, l, d);
                Xlineshape = [originalBreitWigner, formFactorFunc](double s) -> complex {
                    auto ff = formFactorFunc(s);
                    return originalBreitWigner(s) * ff[0] * ff[1];
                };

                lineshapeInitialized = true;


                /// Debug ///////////
                //if(compjulia) std::cout << resonanceName << "FF " << formFactor1 << " " << formFactor2 << std::endl;
                // Get k from topology
                    

                        if(deubevt) EvtGenReport( EVTGEN_INFO, "EvtGen" ) << "INFO" << std::endl;
                
                if(deubevt) std::cout << resonanceName << std::endl;
                if(deubevt) std::cout << "k: " << kint << std::endl;
                if(deubevt) std::cout << mass << " " << width << std::endl;
                if(deubevt) std::cout << "masses: " << ms[0] << " " << ms[1] << " "
                            << ms[2] << " " << ms[3] << std::endl;
                if(deubevt) std::cout << "spins: " << spins[0] << " " << spins[1] << " "
                            << spins[2] << " " << spins[3] << std::endl;


                // Validations point check
                auto Xcheck = func["x"];

                // check if LineshapeName is in "misc": {"amplitude_model_checksums": [
                for ( const auto& check : misc["amplitude_model_checksums"] ) {
                    // Check if the name matches the LineshapeName
                    if ( check["distribution"] == LineshapeName ) {
                        auto valpoint_str = check["point"];
                        auto valvalue_str = check["value"];
                        complex valvalue = evtparseComplex(valvalue_str);
                      
                        // search for valpoint_str in  "parameter_points": [

                        for ( const auto& point : parameterpoints ) {
                            if ( point["name"] == valpoint_str ) {
                                auto valpoint = point["parameters"][0]["value"];
                                complex calculatedValue = originalBreitWigner(valpoint);

                                // check if values match in 1e-6 precision
                                // Use std::abs to compare complex numbers
                                bool realmatches = std::abs(valvalue.real() - calculatedValue.real()) < 1e-6;
                                bool imagmatches = std::abs(valvalue.imag() - calculatedValue.imag()) < 1e-6;
                                if(!(realmatches && imagmatches)) {
                                    EvtGenReport( EVTGEN_ERROR, "EvtGen" ) << "Validation point does not match for: " << LineshapeName << " "
                                              << valpoint << " calculated: "
                                              << calculatedValue << " expected: " << valvalue.real() << " + i*" << valvalue.imag()
                                              << std::endl;
                                    std::cout << "FF" << formFactorFunc(valpoint)[0] << " " << formFactorFunc(valpoint)[1] << std::endl;

                                    
                                }
                                break;
                            }
                            
                        }
                        break;
                    }
                }

               
            }
            else if( func["type"] == "MultichannelBreitWigner" ) {
                double mass = func["mass"];
                auto channels = func["channels"];
         

                //std::cout << channels.dump(4) << std::endl;

                double gsq1 = func["channels"][0]["gsq"];
                double ma1 = func["channels"][0]["ma"];
                double mb1 = func["channels"][0]["mb"];
                int l1 = func["channels"][0]["l"];
                double d1 = func["channels"][0]["d"];

                double gsq2 = func["channels"][1]["gsq"];
                double ma2 = func["channels"][1]["ma"];
                double mb2 = func["channels"][1]["mb"];
                int l2 = func["channels"][1]["l"];
                double d2 = func["channels"][1]["d"];

              


                auto formFactorFunc = createFormFactorFunction(chain, functions, mParent, m1, m2, m3);
    
                

                
                // Modifiziere das Breit-Wigner mit dem Formfaktor
                // make flatte function
                auto originalLineshape = make_flatte(mass, 
                                                       gsq1, ma1, mb1, l1, d1,
                                                       gsq2, ma2, mb2, l2, d2);
                //auto Xlineshapeold = BreitWigner(m
 
                Xlineshape = [originalLineshape, formFactorFunc](double s) -> complex {
                    auto ff = formFactorFunc(s);
                    return originalLineshape(s) * ff[0] * ff[1];
                };
                lineshapeInitialized = true;


                // Debug ///////////

                //if(compjulia) std::cout << resonanceName << "XLineshape: " << Xlineshape(σs[kint-1]) << " with s=" << σs[kint-1]<< std::endl;
                //if(compjulia) std::cout << resonanceName << "FF " << formFactor1 << " " << formFactor2 << std::endl;
                // Get k from topology
                    
                std::cout << resonanceName << "OLineshape: " << Xlineshape(3.2) << " with s=" << 3.2<< std::endl;


                auto Xcheck = func["x"];

                bool docheckforFlatte = false;
                if(docheckforFlatte){
                // check if LineshapeName is in "misc": {"amplitude_model_checksums": [
                for ( const auto& check : misc["amplitude_model_checksums"] ) {
                    // Check if the name matches the LineshapeName
                    if ( check["distribution"] == LineshapeName ) {
                        auto valpoint_str = check["point"];
                        auto valvalue_str = check["value"];
                        complex valvalue = evtparseComplex(valvalue_str);
                      
                        // search for valpoint_str in  "parameter_points": [

                        for ( const auto& point : parameterpoints ) {
                            if ( point["name"] == valpoint_str ) {
                                auto valpoint = point["parameters"][0]["value"];
                                complex calculatedValue = originalLineshape(valpoint);

                                // check if values match in 1e-6 precision
                                // Use std::abs to compare complex numbers
                                bool realmatches = std::abs(valvalue.real() - calculatedValue.real()) < 1e-6;
                                bool imagmatches = std::abs(valvalue.imag() - calculatedValue.imag()) < 1e-6;
                                if(!(realmatches && imagmatches)) {
                                    EvtGenReport( EVTGEN_ERROR, "EvtGen" ) << "Validation point does not match for: " << LineshapeName << " "
                                              << valpoint << " calculated: "
                                              << calculatedValue << " expected: " << valvalue
                                              << std::endl;
                                    std::cout << "FF" << formFactorFunc(valpoint)[0] << " " << formFactorFunc(valpoint)[1] << std::endl;
                                    
                                }
                                break;
                            }
                            
                        }
                        break;
                    }
                }}

            }
            else if( paramType.find( "_BuggBW" ) != std::string::npos ) {

                if(paramType.find( "K700_BuggBW" ) != std::string::npos){
                    Xlineshape = make_Bugg_BW_K700();
                    lineshapeInitialized = true;
                }

                if(paramType.find( "K(700)_BuggBW" ) != std::string::npos){
                    Xlineshape = make_Bugg_BW_K700_xic();
                    lineshapeInitialized = true;
                }



                if(paramType.find( "K1430_BuggBW" ) != std::string::npos) {
                    Xlineshape = make_Bugg_BW_K1430();
                    lineshapeInitialized = true;
                }

                if(paramType.find( "K(1430)_BuggBW" ) != std::string::npos) {
                    Xlineshape = make_Bugg_BW_K1430_xic();
                    lineshapeInitialized = true;
                }
            }
            else if (func["type"] == "BreitWignerMinL") {
                double mass = func["mass"];
                double width = func["width"];
                int l = func["l"];
                std::cout << "l: " << l << std::endl;
                int minl = func["minL"];
                std::cout << "minl: " << minl << std::endl;
                double ma = func["m1"];
                double mb = func["m2"];
                double mm = func["m0"];
                double mk = func["mk"];
                std::cout << "BreitWignerMinL parameters: " << mass << " " << width << " " << l << " " << minl << " " << ma << " " << mb << " " << mm << " " << mk << std::endl;

                auto formFactorFunc = createFormFactorFunction(chain, functions, mParent, m1, m2, m3);

                std::cout << "Form factor function created." << std::endl;
                int kint = topology[1].get<int>(); 

                // Modifiziere das Breit-Wigner mit dem Formfaktor
                auto originalBreitWigner = make_breit_wigner_minl(mass, width, l, minl, ma, mb, mk, mm);
                Xlineshape = [originalBreitWigner, formFactorFunc](double s) -> complex {
                    auto ff = formFactorFunc(s);
                    return originalBreitWigner(s) * ff[0] * ff[1];
                };

                MandelstamTuple sigma = {1.4,3.2,2.634279091379258};
                std::cout << "BreitWignerMinL lineshape created." << originalBreitWigner(sigma[kint-1]) << std::endl;

                lineshapeInitialized = true;
                std::cout << "BreitWignerMinL lineshape initialized." << std::endl;

                /// Debug ///////////
                //if(compjulia) std::cout << resonanceName << "FF " << formFactor1 << " " << formFactor2 << std::endl;
                // Get k from topology
                    

                        if(deubevt) EvtGenReport( EVTGEN_INFO, "EvtGen" ) << "INFO" << std::endl;
                
                if(deubevt) std::cout << resonanceName << std::endl;
                if(deubevt) std::cout << "k: " << kint << std::endl;
                if(deubevt) std::cout << mass << " " << width << std::endl;
                if(deubevt) std::cout << "masses: " << ms[0] << " " << ms[1] << " "
                            << ms[2] << " " << ms[3] << std::endl;
                if(deubevt) std::cout << "spins: " << spins[0] << " " << spins[1] << " "
                            << spins[2] << " " << spins[3] << std::endl;


                // Validations point check
                auto Xcheck = func["x"];
                

                // check if LineshapeName is in "misc": {"amplitude_model_checksums": [
                for ( const auto& check : misc["amplitude_model_checksums"] ) {
                    // Check if the name matches the LineshapeName
                    if ( check["distribution"] == LineshapeName ) {
                        auto valpoint_str = check["point"];
                        auto valvalue_str = check["value"];
                        complex valvalue = evtparseComplex(valvalue_str);
                      
                        // search for valpoint_str in  "parameter_points": [

                        for ( const auto& point : parameterpoints ) {
                            if ( point["name"] == valpoint_str ) {
                                auto valpoint = point["parameters"][0]["value"];
                                complex calculatedValue = originalBreitWigner(valpoint);

                                // check if values match in 1e-6 precision
                                // Use std::abs to compare complex numbers
                                bool realmatches = std::abs(valvalue.real() - calculatedValue.real()) < 1e-6;
                                bool imagmatches = std::abs(valvalue.imag() - calculatedValue.imag()) < 1e-6;
                                if(!(realmatches && imagmatches)) {
                                    EvtGenReport( EVTGEN_ERROR, "EvtGen" ) << "Validation point does not match for: " << LineshapeName << " "
                                              << valpoint << " calculated: "
                                              << calculatedValue << " expected: " << valvalue.real() << " + i*" << valvalue.imag()
                                              << std::endl;
                                    std::cout << "FF" << formFactorFunc(valpoint)[0] << " " << formFactorFunc(valpoint)[1] << std::endl;

                                    
                                }
                                break;
                            }
                            
                        }
                        break;
                    }
                }
            }




            
            
            ThreeBodyParities Ps = {'+', '-', '-', '+'};

            if (!lineshapeInitialized) {
                std::cout << "Warnung: Unbekannter Lineshape-Typ für Resonanz " << resonanceName 
                        << ": " << paramType << std::endl;
                // Fallback-Initialisierung
                Xlineshape = [](double s) -> complex { return complex(1.0, 0.0); };
            }
            
            std::shared_ptr<DecayChain> dc;
            if (helicities.size() >= 2) {
                 dc = createDecayChainCoupling(kint, Xlineshape, jp, tbs, 
                                                helicities[0].recouplingType, 
                                                helicities[0].helicityFactors, 
                                                helicities[0].paritySign,
                                                helicities[1].recouplingType,
                                                helicities[1].helicityFactors, 
                                                helicities[1].paritySign);
            } else {
                std::cout << "ERROR: Need at least 2 helicity configurations!" << std::endl;
            }
            //if(compjulia) std::cout << resonanceName << " k " << kint << " Lineshape " << Xlineshape(σs[kint-1]) << jp << " helicity " << hel[0] << hel[1] << " parity " << par[0] << par[1] << par_sign << std::endl;
            model.add(dc, resonanceName, weight);




            if(compjulia ) {
                MandelstamTuple sigma = {1.5,3.2,1.6714505793792584};
                //sigma = {5.,10.,17.707293};
                sigma = {1.5,3.2,2.5342790913792586};
                sigma = {1.4,3.2,2.634279091379258};
                // print out amplitude for this resonance
                std::cout << resonanceName <<  "XLineshape: " << Xlineshape(sigma[kint-1]) << " with s=" << sigma[kint-1]<< std::endl;
                Tensor4Dcomp amp = tbDecays.amplitude4dcomp(*dc,sigma, 1);
                EvtGenReport( EVTGEN_INFO, "EvtGen" ) << "RESULT" << std::endl;
                std::cout << "Amplitude tensor with weight" << weight << std::endl;
                for ( int i = 0; i < amp.size(); ++i ) {
                    for ( int j = 0; j < amp[0].size(); ++j ) {
                        for ( int k = 0; k < amp[0][0].size(); ++k ) {
                            for ( int z = 0; z < amp[0][0][0].size(); ++z ) {
                                std::cout << i << j << k << z << " " << amp[i][j][k][z] << "\t";
                        }
                    }
                    std::cout << "\n";
                }
                }}


            






        }

    // Set the decay model
    decayModel = std::make_unique<ThreeBodyAmplitudeModel>(model);

    double m12_sq = 3.2;
    double m23_sq = 1.4;
    double m31_sq = 1.77145;


    MandelstamTuple valpoints = {m23_sq, m12_sq, m31_sq}; // Set the Mandelstam variables in the correct order
    //valpoints = {m12_sq, m23_sq, m31_sq}; // Set the Mandelstam variables in the correct order
    MandelstamTuple valpointsJulia = {1.5,3.2,1.6714505793792584};


    double intens = decayModel->intensity(valpoints, 1);



    std::cout << "Initial intensity: " << intens << std::endl;


    if(compjulia ) {
        MandelstamTuple σs = {1.5,3.2,1.6714505793792584};
        //σs = {5.,10.,17.707293};
                        σs = {1.5,3.2,2.5342790913792586};
                        σs = {1.4,3.2,2.634279091379258};


        Tensor4Dcomp amp = decayModel->amplitude4d(σs, 1);
    size_t dim1 = amp.size();
    size_t dim2 = dim1 > 0 ? amp[0].size() : 0;
    size_t dim3 = dim2 > 0 ? amp[0][0].size() : 0;
    size_t dim4 = dim3 > 0 ? amp[0][0][0].size() : 0;

    

    
    if (dim1 == 2 && dim2 == 1 && dim3 == 1 && dim4 == 2) {
        vertex(0, 0, EvtComplex(amp[0][0][0][0].real(), amp[0][0][0][0].imag()));
        vertex(0, 1, EvtComplex(amp[0][0][0][1].real(), amp[0][0][0][1].imag()));
        vertex(1, 0, EvtComplex(amp[1][0][0][0].real(), amp[1][0][0][0].imag()));
        vertex(1, 1, EvtComplex(amp[1][0][0][1].real(), amp[1][0][0][1].imag()));
    }







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
    }


    


    

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


    // Get 4-momenta of all particles
    EvtVector4R pMother = p->getP4();
    EvtVector4R pDaug1 = p->getDaug(0)->getP4();
    EvtVector4R pDaug2 = p->getDaug(1)->getP4();
    EvtVector4R pDaug3 = p->getDaug(2)->getP4();



    double s12 = (pDaug1 + pDaug2).mass2();
    double s23 = (pDaug2 + pDaug3).mass2();
    double s31 = (pDaug3 + pDaug1).mass2();


    double mMother = pMother.mass();
    double m1 = pDaug1.mass();
    double m2 = pDaug2.mass();
    double m3 = pDaug3.mass();

 
   
    // Calculate proper kinematic boundaries
    double s12_min = (m1 + m2) * (m1 + m2);
    double s12_max = (mMother - m3) * (mMother - m3);
    double s23_min = (m2 + m3) * (m2 + m3);
    double s23_max = (mMother - m1) * (mMother - m1);
    double s31_min = (m3 + m1) * (m3 + m1);
    double s31_max = (mMother - m2) * (mMother - m2);

    // Check Dalitz plot constraint: s12 + s23 + s31 = M² + m1² + m2² + m3²
    double sum_constraint = s12 + s23 + s31 - (mMother*mMother + m1*m1 + m2*m2 + m3*m3);
    if (std::abs(sum_constraint) > 1e-6) {
        std::cout << "WARNING: Dalitz constraint violated: " << sum_constraint << std::endl;
    }

    // More restrictive Dalitz plot boundaries
    bool kinematically_allowed = true;
    double tolerance = 1e-6;

    if (s12 < s12_min - tolerance || s12 > s12_max + tolerance ||
        s23 < s23_min - tolerance || s23 > s23_max + tolerance ||
        s31 < s31_min - tolerance || s31 > s31_max + tolerance) {
        kinematically_allowed = false;
    }

    if (!kinematically_allowed) {
        std::cout << "Kinematically forbidden decay!" << std::endl;
        std::cout << "s12: " << s12 << " [" << s12_min << ", " << s12_max << "]" << std::endl;
        std::cout << "s23: " << s23 << " [" << s23_min << ", " << s23_max << "]" << std::endl;
        std::cout << "s31: " << s31 << " [" << s31_min << ", " << s31_max << "]" << std::endl;
        std::cout << "Actual masses: M=" << mMother << ", m1=" << m1 << ", m2=" << m2 << ", m3=" << m3 << std::endl;
        return; // Actually exit here
    }


    MandelstamTuple σs = {s23, s31, s12};
    //σs = {s31, s23, s12}; // Permutation to match the expected order in the model

    /*
    if(compjulia ) {
        σs = {1.5,3.2,1.6714505793792584};
        //σs = {5.,10.,17.707293};
        σs = {1.5,3.2,2.5342790913792586};

    }*/
    // Test point(σ1 = 0.7980703453578917, σ2 = 3.6486261122281745, σ3 = 2.7875826337931926)
    σs = {1.4,3.2,2.634279091379258};



    
    num++;

   
    //std::cout << "Empfohlene maxProb für Kette " << ": " << maxProbEstimate << std::endl;


    


    Tensor4Dcomp amp = decayModel->amplitude4d(σs, 1);
    size_t dim1 = amp.size();
    size_t dim2 = dim1 > 0 ? amp[0].size() : 0;
    size_t dim3 = dim2 > 0 ? amp[0][0].size() : 0;
    size_t dim4 = dim3 > 0 ? amp[0][0][0].size() : 0;

    if(std::isnan(amp[0][0][0][0].real())){
        std::cout << "Amplitude tensor contains NaN values!" << std::endl;
        std::cout << s12 << " " << s23 << " " << s31 << std::endl;

        std::cout << "Kinematically forbidden decay!" << std::endl;
        std::cout << "s12: " << s12 << " [" << s12_min << ", " << s12_max << "]" << std::endl;
        std::cout << "s23: " << s23 << " [" << s23_min << ", " << s23_max << "]" << std::endl;
        std::cout << "s31: " << s31 << " [" << s31_min << ", " << s31_max << "]" << std::endl;

        std::cout << pMother.mass() << " -> " 
                  << pDaug1.mass() << " + " 
                  << pDaug2.mass() << " + " 
                  << pDaug3.mass() << std::endl;

        std::cout << "Mandelstam tuple: " << σs[0] << " " << σs[1] << " " << σs[2] << std::endl;
    }
    

    
    if (dim1 == 2 && dim2 == 1 && dim3 == 1 && dim4 == 2) {
        vertex(0, 0, EvtComplex(amp[0][0][0][0].real(), amp[0][0][0][0].imag()));
        vertex(0, 1, EvtComplex(amp[0][0][0][1].real(), amp[0][0][0][1].imag()));
        vertex(1, 0, EvtComplex(amp[1][0][0][0].real(), amp[1][0][0][0].imag()));
        vertex(1, 1, EvtComplex(amp[1][0][0][1].real(), amp[1][0][0][1].imag()));
    }







    bool ampout = true;
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
    
    

    if(fractionout) calculateIntensities(*decayModel, σs);
    
    createDecayModel( );
    


}







// Neue Funktion - erstellt eine s-abhängige Formfaktor-Funktion
std::function<std::array<double, 2>(double)> EvtThreeBodyDecays::createFormFactorFunction(
    const nlohmann::json& chain,
    const std::map<std::string, nlohmann::json>& functions,
    double mParent, double m1, double m2, double m3)
{
    // Sammle alle Formfaktor-Konfigurationen aus dem Chain
    struct FFConfig {
        double radius;
        int L;
        bool isFirstLevel;
        std::string type;  // ← Neues Feld für Formfaktor-Typ
    };
    
    std::vector<FFConfig> ffConfigs;
    
    for (const auto& vertex : chain["vertices"]) {
        if (vertex.contains("formfactor") && !vertex["formfactor"].get<std::string>().empty()) {
            std::string formfactorName = vertex["formfactor"];
            if (functions.find(formfactorName) != functions.end()) {
                const auto& ff = functions.at(formfactorName);
                
                FFConfig config;
                config.radius = ff.contains("radius") ? ff["radius"].get<double>() : 1.0;
                config.L = ff["l"];
                std::cout << "Formfactor for vertex with L=" << config.L << " and radius=" << config.radius << std::endl;
                config.type = ff.contains("type") ? ff["type"].get<std::string>() : "BlattWeisskopf";
                config.isFirstLevel = vertex["node"].is_array() && vertex["node"][0].is_array();
                
                ffConfigs.push_back(config);
            }
        }
    }
    
    // Lambda-Funktion mit beiden Formfaktor-Typen:
    return [ffConfigs, mParent, m1, m2, m3](double s) -> std::array<double, 2> {
        double formFactor1 = 1.0;
        double formFactor2 = 1.0;
        
        for (const auto& config : ffConfigs) {
            double msub = std::sqrt(s);
            
            if (config.isFirstLevel) {
                double q = breakup(mParent, msub, m3);
                
                if (config.type == "MomentumPower") {
                    formFactor1 *= MomentumPower(q, config.L);
                } else {
                    // Default: Blatt-Weisskopf
                    formFactor1 *= BlattWeisskopf(q, config.L, config.radius);
                }
            } else {
                double q = breakup(msub, m1, m2);
                
                if (config.type == "MomentumPower") {
                    formFactor2 *= MomentumPower(q, config.L);
                } else {
                    // Default: Blatt-Weisskopf
                    formFactor2 *= BlattWeisskopf(q, config.L, config.radius);
                }
            }
        }
        
        return {formFactor1, formFactor2};
    };
}


void EvtThreeBodyDecays::calculateIntensities(ThreeBodyAmplitudeModel model, std::array<double, 3> σs)
{

    /// Calculate the intensities for this decay ///
    const auto& resonance_names = model.names();
    
    const double modelintensity = model.intensity(σs, 1);
    const auto compintens = model.component_intensities(σs, 1);
    totalintensity += modelintensity; // Add model intensity to total intensity
    allmodelintensities.push_back(modelintensity); // Store model intensity for later analysis
    // check that modelintensity is a number
    if (std::isnan(modelintensity)) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Model intensity is NaN!" << std::endl;
        return; // Skip this decay if model intensity is NaN
    }


    int ind = 0;
    for (const auto& name : resonance_names) {
        bool found = false;
        for (size_t i = 0; i < weighttuple.size(); i++) {
            if (name == weighttuple[i].first) {

                weighttuple[i].second.push_back( compintens[ind]); // Add real intensity to existing entry
                found = true;
            }
        }

    
        // If no entry exists for this resonance, add a new one
        if (!found) {
            weighttuple.push_back(std::make_pair(name, std::vector<double>{compintens[ind]}));
        }
        ind++;
    }


    if(num%1000 == 0) {
    
        std::cout << "Component intensities: N = " << allmodelintensities.size() << std::endl;
        
        // Mittelwert der Modellintensitäten berechnen
        double model_mean = std::accumulate(allmodelintensities.begin(), allmodelintensities.end(), 0.0) 
                           / allmodelintensities.size();
        
        // Über weighttuple iterieren
        for (const auto& entry : weighttuple) {
            const std::string& name = entry.first;
            const std::vector<double>& intensities = entry.second;
            
            if (!intensities.empty()) {
                // Mittelwert der Intensitäten berechnen
                double mean = std::accumulate(intensities.begin(), intensities.end(), 0.0) 
                            / allmodelintensities.size();
                
                // Varianz berechnen
                double variance = 0.0;
                for (const auto& value : intensities) {
                    variance += (value - mean) * (value - mean);
                }
                variance /= intensities.size();
                
                // Standardfehler berechnen
                double stddev = std::sqrt(variance / allmodelintensities.size());
                
                // Ausgabe wie im ursprünglichen Code, aber mit Prozentangabe
                std::cout << name << " Mean: " << (mean/model_mean)*100 << " ± " 
                          << (stddev/model_mean)*100 << ", N=" << allmodelintensities.size() << std::endl;
            } else {
                // Fallback für leere Intensitätsvektoren
                std::cout << name << " Intensity: 0.0 (No variance data)" << std::endl;
            }
        }
    
    }
}