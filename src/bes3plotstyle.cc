#include <TH1.h>
#include <TH1F.h>
#include <TF1.h>
#include <TLatex.h>
#include <TList.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TGraph.h>
#include <math.h> 
#include <iostream>
#include "TPad.h"

#include "bes3plotstyle.hh"
// bes3plotstyle::Format for data poInt_ts
/*
   void bes3plotstyle::bes3plotstyle::Format(RooPlot *frame){
   bes3plotstyle::FormatAxis(frame->GetXaxis());
   bes3plotstyle::FormatAxis(frame->GetYaxis());
   }
   */
void bes3plotstyle::Format(TPad *pad){
    pad->SetBottomMargin(0.15);
    pad->SetLeftMargin(0.14);
    pad->SetRightMargin(0.15);
    pad->SetTopMargin(0.12);
    pad->SetFillColor(0);
}
void bes3plotstyle::FormatData(TH1 * datahist){
    datahist->SetMarkerStyle(20);
    datahist->SetMarkerSize(1);
    datahist->SetLineWidth(2);

    bes3plotstyle::FormatAxis(datahist->GetXaxis());
    bes3plotstyle::FormatAxis(datahist->GetYaxis());
}

// bes3plotstyle::Format for graph data poInt_ts
void bes3plotstyle::FormatData(TGraph * datahist){
    datahist->SetMarkerStyle(20);
    datahist->SetMarkerSize(1);
    datahist->SetLineWidth(2);
}
void bes3plotstyle::Format(TCanvas * c){
    c->SetRightMargin(0.15);
    c->SetLeftMargin(0.2);
    c->SetTopMargin(0.15);
    c->SetBottomMargin(0.25);
}
void bes3plotstyle::Format(TPaveText *pt){
    pt->SetTextFont(132);
    pt->SetTextAlign(12);  //zuo dui qi
    pt->SetBorderSize(0);
    pt->SetFillColor(0);
}

void bes3plotstyle::FormatAxis(TAxis * axis){
    axis->SetLabelFont(42);
    axis->SetLabelSize(0.05);
    axis->SetLabelOffset(0.02);
    axis->SetNdivisions(505);
    axis->SetTitleFont(42);
    axis->SetTitleColor(1);
    axis->SetTitleSize(0.05);
    axis->SetNoExponent(kFALSE);
    axis->CenterTitle();
}
void bes3plotstyle::SetInSPad(TAxis * axis){
    axis->SetLabelFont(42);
    axis->SetLabelSize(0.10);
    axis->SetLabelOffset(0.02);
    axis->SetNdivisions(505);
    axis->SetTitleFont(42);
    axis->SetTitleColor(1);
    axis->SetTitleSize(0.15);
    axis->SetNoExponent(kFALSE);
    axis->CenterTitle();
}
void bes3plotstyle::FormatYAxis(TAxis * axis){
    axis->SetLabelFont(42);
    axis->SetLabelSize(0.05);
    axis->SetLabelOffset(0.01);
    axis->SetNdivisions(505);
    axis->SetTitleFont(42);
    axis->SetTitleOffset(1.2);
    axis->SetTitleColor(1);
    axis->SetTitleSize(0.05);
    axis->SetNoExponent(kFALSE);
    axis->CenterTitle();
}

void bes3plotstyle::NameAxes(TH1 * datahist, char * xname, char * yname){
    if(xname)
        datahist->GetXaxis()->SetTitle(xname);
    if(yname)
        datahist->GetYaxis()->SetTitle(yname);
}

// bes3plotstyle::Format for main MC (red line)
void bes3plotstyle::FormatMC1(TH1 * mc1hist){
    mc1hist->SetLineColor(2);
    mc1hist->SetLineWidth(2);
}
//fitting result
// Graph bes3plotstyle::Format for main MC (red line)
void bes3plotstyle::FormatMC1(TGraph * mc1hist){
    mc1hist->SetLineColor(4);
    mc1hist->SetLineWidth(2);
}

//backgrond
// bes3plotstyle::Format for second MC or background
// (Blue shaded area)
void bes3plotstyle::FormatMC2(TH1 * mc2hist){
    mc2hist->SetLineColor(3);
    mc2hist->SetFillColor(3);
    mc2hist->SetLineWidth(2);
    mc2hist->SetLineStyle(9);
    mc2hist->SetFillStyle(3001);
}
void bes3plotstyle::Format(TLegend * leg){
    leg->SetTextFont(42);
    leg->SetFillColor(0);
}
void bes3plotstyle::Format(TArrow * arr){
    arr->SetLineWidth(2);
    arr->SetLineColor(kRed);
}
// Graph bes3plotstyle::Format for second MC or background
// (Blue line)
void bes3plotstyle::FormatMC2(TGraph * mc2hist){
    mc2hist->SetLineColor(3);
    mc2hist->SetLineWidth(2);
    mc2hist->SetLineStyle(9);
    mc2hist->SetFillColor(3);
    //  mc2hist->SetFillStyle(3001);
}
//signal function
// Graph bes3plotstyle::Format for third MC or background
// (Blue line)
void bes3plotstyle::FormatMC3(TGraph * mc3hist){
    mc3hist->SetLineColor(2);
    mc3hist->SetLineWidth(2);
    mc3hist->SetLineStyle(7);
}

void bes3plotstyle::FormatMC4(TGraph * mc4hist){
    mc4hist->SetLineColor(2);
    mc4hist->SetLineWidth(2);
    mc4hist->SetLineStyle(9);
}

// Write "BESIII" in the upper right corner
void bes3plotstyle::WriteBes3(){
    TLatex * bes3 = new TLatex(0.94,0.94, "BESIII");
    bes3->SetNDC();
    bes3->SetTextFont(72);
    bes3->SetTextSize(0.1);
    bes3->SetTextAlign(33);
    bes3->Draw();
}

// Write "Preliminary" below BESIII -
// to be used together with WriteBes3()
void bes3plotstyle::WritePreliminary(){
    TLatex * prelim = new TLatex(0.94,0.86, "Preliminary");
    prelim->SetNDC();
    prelim->SetTextFont(62);
    prelim->SetTextSize(0.055);
    prelim->SetTextAlign(33);
    prelim->Draw();
}

// Make a legend; 
// position will have to change depending on the data shape
void bes3plotstyle::MakeLegend(TH1 * datahist,   // Histogram with data
        char * dataname,  // Description of data
        TH1 * mc1hist, // Histogram with first MC
        char * mc1name, // Description of first MC
        TH1 * mc2hist, // Histogram with 2nd MC/BG
        char * mc2name, // Description of second MC/BG
        Double_t xlow,      // Left edge of legend 
        //(fraction of canavas width)
        Double_t ylow,       // Bottom edge of legend
        //(fraction of canavas height)
        Double_t xhi,       // Right edge of legend 
        //(fraction of canavas width)
        Double_t yhi){       // Top edge of legend
    //(fraction of canavas height)

    TLegend * leg = new TLegend(xlow, ylow, xhi, yhi);
    if(datahist && dataname)
        leg->AddEntry(datahist, dataname, "LEP");
    if(mc1hist && mc1name)
        leg->AddEntry(mc1hist, mc1name, "L");
    if(mc2hist && mc2name)
        leg->AddEntry(mc2hist, mc2name, "LF");

    leg->SetFillColor(0);
    leg->SetTextFont(42);
    leg->Draw();

}


// Make a legend; 
// position will have to change depending on the data shape
void bes3plotstyle::MakeLegend(TGraph * datahist,   // Graph with data
        char * dataname,  // Description of data
        TGraph * mc1hist, // Graph with first MC
        char * mc1name, // Description of first MC
        TGraph * mc2hist, // Graph with 2nd MC/BG
        char * mc2name, // Description of second MC/BG
        TGraph * mc3hist, // Graph with 3rd MC/BG
        char * mc3name, // Description of third MC/BG

        Double_t xlow,      // Left edge of legend 
        //(fraction of canavas width)
        Double_t ylow,       // Bottom edge of legend
        //(fraction of canavas height)
        Double_t xhi,       // Right edge of legend 
        //(fraction of canavas width)
        Double_t yhi){       // Top edge of legend
    //(fraction of canavas height)

    TLegend * leg = new TLegend(xlow, ylow, xhi, yhi);
    if(datahist && dataname)
        leg->AddEntry(datahist, dataname, "LEP");
    if(mc1hist && mc1name)
        leg->AddEntry(mc1hist, mc1name, "L");
    if(mc2hist && mc2name)
        leg->AddEntry(mc2hist, mc2name, "L");
    if(mc3hist && mc3name)
        leg->AddEntry(mc3hist, mc3name, "L");

    leg->SetFillColor(0);
    leg->SetTextFont(42);
    leg->Draw();

}


// Make a legend (version for fit functions
// position will have to change depending on the data shape
void bes3plotstyle::MakeLegend(TH1 * datahist,   // Histogram with data
        char * dataname,  // Description of data
        char ** functionnames, // list of function names
        Double_t xlow,      // Left edge of legend 
        //(fraction of canavas width)
        Double_t ylow,       // Bottom edge of legend
        //(fraction of canavas height)
        Double_t xhi,       // Right edge of legend 
        //(fraction of canavas width)
        Double_t yhi){       // Top edge of legend
    //(fraction of canavas height)

    TLegend * leg = new TLegend(xlow, ylow, xhi, yhi);
    if(datahist && dataname)
        leg->AddEntry(datahist, dataname, "LEP");

    TList* list = datahist->GetListOfFunctions();
    Int_t nfun = list->GetEntries();

    for(Int_t i =0;  i < nfun; i++){
        TF1* f1 = (TF1*)(list->At(i));
        leg->AddEntry(f1, functionnames[i], "L");
    }
    leg->SetFillColor(0);
    leg->SetTextFont(42);
    leg->Draw();

}





// Set the general style options
void bes3plotstyle::SetStyle(){
    // No Canvas Border
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetCanvasBorderSize(0);
    // White BG
    gStyle->SetCanvasColor(10);
    // bes3plotstyle::Format for axes
    gStyle->SetLabelFont(42,"xyz");
    gStyle->SetLabelSize(0.04,"xyz");
    gStyle->SetLabelOffset(0.02,"y");
    gStyle->SetLabelOffset(0.02,"x");
    gStyle->SetNdivisions(510,"xyz");
    gStyle->SetTitleFont(42,"xyz");
    gStyle->SetTitleColor(1,"xyz");
    gStyle->SetTitleSize(0.09,"xyz");
    gStyle->SetTitleOffset(1.35,"x");
    gStyle->SetTitleOffset(1.09,"y");
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetHistLineWidth(1.85);
    // No pad borders
    gStyle->SetPadBorderMode(0);
    gStyle->SetPadBorderSize(0);
    // White BG
    gStyle->SetPadColor(10);
    // Margins for labels etc.
    gStyle->SetPadLeftMargin(0.25);
    gStyle->SetPadBottomMargin(0.25);
    gStyle->SetPadRightMargin(0.12);
    gStyle->SetPadTopMargin(0.12);
    // No error bars in x direction
    gStyle->SetErrorX(0);
    // bes3plotstyle::Format legend
    gStyle->SetLegendBorderSize(0);
    //
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
}

// Style options for "final" plots
// (no stat/fit box)
void bes3plotstyle::SetPrelimStyle(){
    gStyle->SetOptDate(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetOptTitle(0);
}

// Style options for Int_ternal meetings
// (stat/fit box)
void bes3plotstyle::SetMeetingStyle(){
    gStyle->SetOptDate(0);
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(1111);
    gStyle->SetStatBorderSize(1);
    gStyle->SetStatColor(10);
    gStyle->SetStatFont(42);
    gStyle->SetStatFontSize(0.03);
    gStyle->SetOptFit(1111);
}


// Plot a data MC plot
void bes3plotstyle::PlotDataMC(char * filename,  // Name for the output files, 
        // without extension 
        TH1 * datahist,   // Histogram with data
        char * dataname,  // Description of data
        TH1 * mc1hist, // Histogram with first MC
        char * mc1name, // Description of first MC
        TH1 * mc2hist, // Histogram with 2nd MC/BG
        char * mc2name, // Description of second MC/BG
        Int_t prelim,    // Use 1 for Preliminary plot
        // 2 for a publication plot
        // and 0 for a meeting plot with 
        // stat and fit box
        Double_t xlow, // Left edge of legend 
        //(fraction of canavas width)
        Double_t ylow,  // Bottom edge of legend
        //(fraction of canavas height)
        Double_t xhi,  // Right edge of legend 
        //(fraction of canavas width)
        Double_t yhi){  // Top edge of legend
    //(fraction of canavas height)

    SetStyle();
    if(prelim)
        SetPrelimStyle();
    else
        SetMeetingStyle();

    TCanvas * c1 = new TCanvas("bes3plots","BESIII Plots", 800,600);

    bes3plotstyle::FormatData(datahist);
    if(mc1hist)
        bes3plotstyle::FormatMC1(mc1hist);
    if(mc2hist)
        bes3plotstyle::FormatMC2(mc2hist);


    datahist->Draw("axis");
    if(mc2hist)
        mc2hist->Draw("same");
    if(mc1hist)
        mc1hist->Draw("same");
    datahist->Draw("Esame");
    datahist->Draw("axissame");
    if(prelim){
        WriteBes3();
        if(prelim == 1)
            WritePreliminary();
    }
    MakeLegend(datahist, dataname,
            mc1hist,  mc1name,
            mc2hist,  mc2name);


    c1->SaveAs(filename + TString(".eps") );
    c1->SaveAs(filename + TString(".png") );
}


// Plot data with one or more (fit) functions
// Functions should be part of the data histograms list of functions
// (i.e. perform fits with the "+" option or add other functions via
// datahist->GetListOfFunctions->Add(TF1 * function))
// functionnames should have at least as many elements as the function
// list
void bes3plotstyle::PlotDataFit(char * filename,  // Name for the output files, 
        // without extension 
        TH1F * datahist,   // Histogram with data
        char * dataname,  // Description of data
        char ** functionnames,// Names of associated functions
        Int_t prelim,    // Use 1 for Preliminary plot
        // 2 for a publication plot
        // and 0 for a meeting plot with 
        // stat and fit box
        Double_t xlow, // Left edge of legend 
        //(fraction of canavas width)
        Double_t ylow,  // Bottom edge of legend
        //(fraction of canavas height)
        Double_t xhi,  // Right edge of legend 
        //(fraction of canavas width)
        Double_t yhi){  // Top edge of legend
    //(fraction of canavas height)


    SetStyle();
    if(prelim)
        SetPrelimStyle();
    else
        SetMeetingStyle();

    TCanvas * c1 = new TCanvas("bes3plots","BESIII Plots", 800,600);

    bes3plotstyle::FormatData(datahist);

    Int_t linestyles[] = {1,2,3,7,9,10};
    //Int_t linecolors[]   = {2,4,kGreen+2,kOrange+7,kMagenta,2};
    Int_t linecolors[]   = {2,4,6, 9,8,2};

    TList* list = datahist->GetListOfFunctions();
    TH1F * datacopy = new TH1F(*datahist); 
    datacopy->Draw("axis");


    Int_t nfun = list->GetEntries();

    if(nfun > 6){
        std::cout << "ERROR: More than six associated functions not forseen" << std::endl;
        return;
    }


    for(Int_t i =0;  i < nfun; i++){
        TF1* f1 = (TF1*)(list->At(i));
        f1->SetLineColor(linecolors[i]);
        f1->SetLineStyle(linestyles[i]);
        f1->Draw("same");
    }

    MakeLegend(datahist, dataname, functionnames,xlow, ylow, xhi, yhi);

    datacopy->Draw("Esame");
    datacopy->Draw("axissame");

    if(prelim){
        WriteBes3();
        if(prelim==1)
            WritePreliminary();
    }

    c1->SaveAs(filename + TString(".png") );;
    c1->SaveAs(filename + TString(".eps") );;
}


// Scatter plot
void bes3plotstyle::PlotScatter(char * filename,  // Name for the output files, 
        // without extension 
        TH1 * datahist,   // Histogram with data
        Int_t prelim       // preliminary plot
        ){

    SetStyle();
    if(prelim)
        SetPrelimStyle();
    else
        SetMeetingStyle();

    TCanvas * c1 = new TCanvas("bes3plots","BESIII Plots", 800,600);

    bes3plotstyle::FormatData(datahist);

    if(datahist->Integral() > 5000)
        datahist->SetMarkerStyle(1);
    else if(datahist->Integral() > 500)
        datahist->SetMarkerSize(0.5);


    datahist->Draw("");

    if(prelim){
        WriteBes3();
        if(prelim==1)
            WritePreliminary();
    }


    c1->SaveAs(filename + TString(".eps") );
    c1->SaveAs(filename + TString(".png") );
}
/*
   void bes3plotstyle::DrawSigArea(Double_t xl,Double_t xh,Double_t yl,Double_t yh){
   TLine *l1 = new TLine(xl,yl,xl,yh);
   TLine *l2 = new TLine(xh,yl,xh,yh);
   TLine *l3 = new TLine(xl,yl,xh,yl);
   TLine *l4 = new TLine(xl,yh,xh,yh);
   l1->SetLineColor(3); 
   l2->SetLineColor(3); 
   l3->SetLineColor(3); 
   l4->SetLineColor(3); 
   l1->SetLineWidth(2); 
   l2->SetLineWidth(2); 
   l3->SetLineWidth(2); 
   l4->SetLineWidth(2); 
   l1->Draw();
   l2->Draw("same");
   l3->Draw("same");
   l4->Draw("same");
   }
   void bes3plotstyle::DrawSide(Double_t xl,Double_t xh,Double_t yl,Double_t yh){
   TLine *l1 = new TLine(xl,yl,xl,yh);
   TLine *l2 = new TLine(xh,yl,xh,yh);
   TLine *l3 = new TLine(xl,yl,xh,yl);
   TLine *l4 = new TLine(xl,yh,xh,yh);
   l1->SetLineColor(2); 
   l2->SetLineColor(2); 
   l3->SetLineColor(2); 
   l4->SetLineColor(2); 
   l1->SetLineWidth(2); 
   l2->SetLineWidth(2); 
   l3->SetLineWidth(2); 
   l4->SetLineWidth(2); 
   l1->SetLineStyle(2); 
   l2->SetLineStyle(2); 
   l3->SetLineStyle(2); 
   l4->SetLineStyle(2); 
   l1->Draw();
   l2->Draw("same");
   l3->Draw("same");
   l4->Draw("same");
   }
   Double_t sigma(RooRealVar sigma1, RooRealVar sigma2, RooRealVar mean1, RooRealVar mean2 , RooRealVar frac){
   Double_t s1 = sigma1.getVal();   
   Double_t s2 = sigma2.getVal();   
   Double_t m1 = mean1.getVal();   
   Double_t m1 = mean2.getVal();   
   Double_t f = frac.getVal();   
   Double_t e1 = f*(s1*s1+m1*m1);
   Double_t e2 = (1-f)*(s2*s2+m2*m2);
   Double_t e3 = f*m1+(1-f)*m2;
   Double_t error = e1 + e2 -e3*e3;
   return sqrt(error);
   }
   Double_t mean(RooRealVar sigma1, RooRealVar sigma2, RooRealVar mean1, RooRealVar mean2 , RooRealVar frac){
   Double_t s1 = sigma1.getVal();   
   Double_t s2 = sigma2.getVal();   
   Double_t m1 = mean1.getVal();   
   Double_t m1 = mean2.getVal();   
   Double_t f = frac.getVal();   
   return f*m1+f*m2;
   }
   void calculateChisq(RooPlot* frame, Double_t& m_chisq, Int_t& m_ndf) {
   RooCurve *curve = (RooCurve*)frame->getObject(1);
   RooHist *histo = (RooHist*)frame->getObject(0);
   Double_t x[1000], y[1000];
   Double_t xa[1000], ya[1000];
   Double_t y_fit[1000],x_fit[1000];
   Double_t tmp_y = 0;
   Double_t eyl[1000], eyh[1000], exl[1000], exh[1000];
   Int_t ndf = 0;
   for (Int_t i = 0; i<histo->GetN(); i++) {
histo->GetPoInt_t(i,x[i],y[i]);
eyl[i] = histo->GetEYlow()[i];
eyh[i] = histo->GetEYhigh()[i];
exl[i] = histo->GetEXlow()[i];
exh[i] = histo->GetEXhigh()[i];
}
for (Int_t i = 0; i<histo->GetN(); i++) {
    tmp_y += y[i];
    if ( tmp_y < 10 ) continue;
    xa[ndf] = x[i];
    ya[ndf] = tmp_y;
    tmp_y = 0;
    ndf++;
}
if ( tmp_y != 0 ) {
    xa[ndf] = x[histo->GetN()-1];
    ya[ndf] = tmp_y;

}
Double_t xstart, xstop, ystart, ystop;
curve->GetPoInt_t(0,xstart,ystart);
curve->GetPoInt_t(curve->GetN()-1,xstop,ystop);
for (Int_t i = 0; i<histo->GetN(); i++) {
    if (x[i] > xstop || x[i] < xstart ) continue;
    x_fit[i] = x[i];
    y_fit[i] = curve->average(x[i]-exl[i],x[i]+exh[i]);
}

Double_t y_fb[1000];
Int_t num = 0;
for (Int_t i=0; i<histo->GetN(); i++) {
    if ( y[i] < 1 ) continue;

    if ( x[i] <= xa[num] ) {
        y_fb[num] += y_fit[i];
    }
    else {
        num++;
        y_fb[num] += y_fit[i];
    }
}
Double_t tx[1000]={0};
Double_t tx_err[1000]={0};
Double_t ty[1000]={0};
Double_t chi_err[1000]={0};
Double_t chisq = 0;

for ( Int_t tk = 0; tk != histo->GetN(); tk++){
    tx[tk]=x[tk];
    if (y[tk] == 0) continue;
    if (y[tk] > y_fit[tk] ){
        ty[tk]=(y[tk] - y_fit[tk])/eyl[tk];
    }
    else {
        ty[tk]=(y[tk] - y_fit[tk])/eyh[tk];
    }
    if (fabs(ty[tk])<10.) chisq += ty[tk] * ty[tk];
}
m_chisq = chisq;
m_ndf = ndf;
}
void CalculateResidual( RooPlot *frame, Int_t &bins, Double_t xx[], Double_t *yy, 
        Double_t *yErrL, Double_t *yErrH){
    // the data should be  ploted first
    RooHist *histo = (RooHist*)frame->getObject(0);
    // then the model plot
    RooCurve *curve = (RooCurve*)frame->getObject(1);
    Int_t _npionts = histo->GetN();
    cout<<"the pionts in dataplot \t"<<_npionts<<endl;
    Double_t _xx [1000];
    Double_t _yy [1000];
    Double_t _xErr  [1000];
    Double_t _yErrL [1000];
    Double_t _yErrH [1000];

    // get the poInt_t in the frame
    for (Int_t i = 0; i < _npionts; i++) {
        histo->GetPoInt_t(i, _xx[i], _yy[i]);
        _yErrL[i] = histo->GetEYlow()[i];
        _yErrH[i] = histo->GetEYhigh()[i];
        //exl[i] = histo->GetEXlow()[i];
        //exh[i] = histo->GetEXhigh()[i];
    }
    // add the bin which contens less than 10 togethe
    Double_t _x_data [1000];
    Double_t _y_data [1000];
    Double_t _yErrL_data [1000];
    Double_t _yErrH_data [1000];

    Int_t kk =0; // index
    for(Int_t i =0 ; i< _npionts; i++){
        _x_data[i] = 0.0; 
        _y_data[i] = 0.0; 
        _yErrL_data[i] = 0.0; 
        _yErrH_data[i] = 0.0; 
    }
    for(Int_t i =0 ; i< _npionts; i++){
        // take the average of value of x as the x_value in data
        _x_data[kk] = ( _y_data[kk] * _x_data[kk] +  _xx[i] * _yy[i]) /
            ( _yy[i] + _y_data[kk]);
        _y_data[kk] += _yy[i];
        _yErrL_data[kk] = sqrt( _yErrL_data[kk]**2 + _yErrL[i]**2 );
        _yErrH_data[kk] = sqrt( _yErrH_data[kk]**2 + _yErrH[i]**2 );
        //cout<<"x,y="<<_x_data[kk]<<"\t"<<_y_data[kk]<<endl;
        // if( _y_data[kk] < 10){
        //     continue;
        // }
        // else{
        kk++;
        continue;
        //}
    }
    _npionts = kk;
    cout<<"new poInt_ts=\t"<<_npionts<<endl;
    // get the nearest poInt_t in the curve
    cout<<"total\t"<<curve->GetN()<<endl;
    Double_t _y_fit[1000];
    for(Int_t i=0;i<curve->GetN();i++){
        Double_t x1, x2, y1, y2;
        curve->GetPoInt_t(i, x1, y1);
        cout<<"index\t"<<i<<endl;
        cout<<"x,y=\t"<<x1<<"\t"<<y1<<endl;
    }
    for(Int_t i = 0; i < _npionts; i++ ){
        Double_t x1, x2, y1, y2;
        cout<<"i="<<i<<"\t"<<endl;
        Int_t idx =   curve->findPoInt_t(_x_data[i], 0.1);
        cout<<"index\t"<<idx<<endl;
        curve->GetPoInt_t(idx, x1, y1);
        cout<<"data\t"<<_x_data[i]<<" "<<_y_data[i]<<" fit\t"<<x1<<" "<<y1<<endl;
        curve->GetPoInt_t(idx, x1, y1);
        if(idx ==0){
            curve->GetPoInt_t(idx+1, x2, y2);
        }
        else if( idx >= curve->GetN()){
            curve->GetPoInt_t(idx-1, x2, y2);
        }
        else if(x1 < _xx[i]){
            curve->GetPoInt_t(idx+1, x2, y2);
        }
        else{
            curve->GetPoInt_t(idx-1, x2, y2);
        }
        Double_t x_fit = _x_data[i];
        _y_fit[i] = (y2 *( x_fit - x1 ) + y1* ( x2 - x_fit) ) / (x2 - x1);
        if(x_fit<1.9){
            x_fit = _x_data[i-1];
            _y_fit[i] = 0;
        }
    }


    // save value
    cout<<"----------------------- Save Value----------"<<endl;
    bins = _npionts;
    for(Int_t i = 0 ; i < bins ; i++ ){
        cout<<"sigma Y\t"<<sqrt(_y_data[i])<<endl;
        cout<<"i=\t"<<i<< _x_data[i]<<endl;
        xx[i] = _x_data[i];
        if( _y_data[i] ==0 ){ 
            xx[i] = xx[i-1];
            yErrH[i] = yErrH[i-1];
            yErrL[i] = yErrH[i-1];
            yy[i] = yy[i-1];
            continue;
        }

        yy[i] = ( _y_data[i] - _y_fit[i]) / sqrt( _y_data[i] );
        yErrL[i] = ( _yErrH_data[i]) / sqrt( _y_data[i] );
        yErrH[i] = ( _yErrH_data[i]) / sqrt( _y_data[i] );
    }
}


TString cut_On_Topo(Int_t* topos,Int_t size){
    TString m_cut = "";
    char cut[15]="";
    sprInt_tf(cut,"n_itopo==%d",topos[0]);
    m_cut +=cut;

    for(Int_t i=1;i<size;i++){
        char cut[15]="";
        sprInt_tf(cut,"||n_itopo==%d",topos[i]);
        m_cut +=cut;
    }
    return m_cut;
}
void getSigmaAndMean(RooAbsPdf *sigPdf, TString xname ="deltaE", Double_t &mean, Double_t &sigma){
    //sigPdf->PrInt_t("v");
    RooArgSet *params = sigPdf->getVariables();
    RooRealVar* x = (RooRealVar*)params->find(xname);

    RooDataSet* data_toy = sigPdf->generate(*x,1E6) ;

    mean =  data_toy->mean(*x);
    sigma = data_toy->sigma(*x);
}
*/
