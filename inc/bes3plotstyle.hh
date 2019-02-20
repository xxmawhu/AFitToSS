#ifndef bes3plotstyle__H
#define bes3plotstyle__H

#include "TH1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TArrow.h"
namespace bes3plotstyle {
    // Format for data poInt_ts
    void FormatData(TH1 * datahist);
    void Format(TPad *);
    // Format for graph data poInt_ts
    void FormatData(TGraph * datagraph);
    void Format(TCanvas *c);
    void Format(TPaveText *pt);
    // Format Axis
    void FormatAxis(TAxis *axis);
    void SetInSPad(TAxis * axis);
    void FormatYAxis(TAxis * axis);
    void Format(TLegend * leg);
    void Format(TArrow * arr);
    void FormatMC1(TH1 * mc1hist);
    void FormatMC1(TGraph * mc1hist);
    void FormatMC2(TH1 * mc2hist);
    void FormatMC2(TGraph * mc2hist);
    void FormatMC3(TGraph * mc3hist);
    void FormatMC4(TGraph * mc4hist);
    /*
    // Graph format for main MC (red line)
    // Format for second MC or background
    // (Blue shaded area)
    // Graph Format for second MC or background
    // (Blue line)
    // Graph Format for third MC or background
    // (Blue line)
    void FormatMC3(TH1 * mc3hist);

    void FormatMC4(TH1 * mc4hist);
    */

    // Name histogram axes
    void NameAxes(TH1 * datahist, char * xname, char * yname);

    // Write "BESIII" in the upper right corner
    void WriteBes3();
    // Write "Preliminary" below BESIII -
    // to be used together with WriteBes3()
    void WritePreliminary();

    // Make a legend; 
    // position will have to change depending on the data shape
    void MakeLegend(TH1 * datahist,   // Histogram with data
            char * dataname,  // Description of data
            TH1 * mc1hist =0, // Histogram with first MC
            char * mc1name=0, // Description of first MC
            TH1 * mc2hist =0, // Histogram with 2nd MC/BG
            char * mc2name=0, // Description of second MC/BG
            Double_t xlow = 0.55,      // Left edge of legend 
            //(fraction of canavas width)
            Double_t ylow = 0.5,       // Bottom edge of legend
            //(fraction of canavas height)
            Double_t xhi = 0.94,       // Right edge of legend 
            //(fraction of canavas width)
            Double_t yhi = 0.7);       // Top edge of legend
    //(fraction of canavas height)

    // Make a legend; 
    // position will have to change depending on the data shape
    void MakeLegend(TGraph * datahist,   // Graph with data
            char * dataname,  // Description of data
            TGraph * mc1hist =0, // Graph with first MC
            char * mc1name=0, // Description of first MC
            TGraph * mc2hist =0, // Graph with 2nd MC/BG
            char * mc2name=0, // Description of second MC/BG
            TGraph * mc3hist =0, // Graph with 3rd MC/BG
            char * mc3name=0, // Description of third MC/BG
            Double_t xlow = 0.55,      // Left edge of legend 
            //(fraction of canavas width)
            Double_t ylow = 0.5,       // Bottom edge of legend
            //(fraction of canavas height)
            Double_t xhi = 0.94,       // Right edge of legend 
            //(fraction of canavas width)
            Double_t yhi = 0.7);       // Top edge of legend
    //(fraction of canavas height)


    // Make a legend (version for fit functions
    // position will have to change depending on the data shape
    void MakeLegend(TH1 * datahist,   // Histogram with data
            char * dataname,  // Description of data
            char ** functionnames, // list of function names
            Double_t xlow = 0.55,      // Left edge of legend 
            //(fraction of canavas width)
            Double_t ylow = 0.5,       // Bottom edge of legend
            //(fraction of canavas height)
            Double_t xhi = 0.94,       // Right edge of legend 
            //(fraction of canavas width)
            Double_t yhi = 0.7);       // Top edge of legend
    //(fraction of canavas height)

    // Set the general style options
    void SetStyle();

    // Style options for "final" plots
    // (no stat/fit box)
    void SetPrelimStyle();

    // Style options for Int_ternal meetings
    // (stat/fit box)
    void SetMeetingStyle();

    // Plot a data MC plot
    void PlotDataMC(char * filename,  // Name for the output files, 
            // without extension 
            TH1 * datahist,   // Histogram with data
            char * dataname,  // Description of data
            TH1 * mc1hist =0, // Histogram with first MC
            char * mc1name=0, // Description of first MC
            TH1 * mc2hist =0, // Histogram with 2nd MC/BG
            char * mc2name=0, // Description of second MC/BG
            Int_t prelim = 1,    // Use 1 for Preliminary plot
            // 2 for a publication plot
            // and 0 for a meeting plot with 
            // stat and fit box
            Double_t xlow = 0.55, // Left edge of legend 
            //(fraction of canavas width)
            Double_t ylow = 0.5,  // Bottom edge of legend
            //(fraction of canavas height)
            Double_t xhi = 0.94,  // Right edge of legend 
            //(fraction of canavas width)
            Double_t yhi = 0.7);  // Top edge of legend
    //(fraction of canavas height)

    // Plot data with one or more (fit) functions
    // Functions should be part of the data histograms list of functions
    // (i.e. perform fits with the "+" option or add other functions via
    // datahist->GetListOfFunctions->Add(TF1 * function))
    // functionnames should have at least as many elements as the function
    // list
    void PlotDataFit(char * filename,  // Name for the output files, 
            // without extension 
            TH1F * datahist,   // Histogram with data
            char * dataname,  // Description of data
            char ** functionnames,// Names of associated functions
            Int_t prelim = 1,    // Use 1 for Preliminary plot
            // 2 for a publication plot
            // and 0 for a meeting plot with 
            // stat and fit box
            Double_t xlow = 0.35, // Left edge of legend 
            //(fraction of canavas width)
            Double_t ylow = 0.5,  // Bottom edge of legend
            //(fraction of canavas height)
            Double_t xhi = 0.94,  // Right edge of legend 
            //(fraction of canavas width)
            Double_t yhi = 0.7);  // Top edge of legend
    //(fraction of canavas height)

    // Scatter plot
    void PlotScatter(char * filename,  // Name for the output files, 
            TH1 * datahist,   // Histogram with data
            Int_t prelim = 1);
};
#endif
