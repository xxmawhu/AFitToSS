#ifndef RooSSCSPDF
#define RooSSCSPDF

#include <math.h>
#include <iomanip>
#include <iostream> 
#include <fstream>
#include "TComplex.h"
#include "TString.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TStopwatch.h"
#include "TH1F.h"
#if defined(USEROOT) || defined(__CINT__)
#include "RooStringVar.h"
#include "RooAbsPdf.h"
#include "RooAbsData.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "RooFitResult.h"
#include "RooArgSet.h"
#include "RooRealConstant.h"
//#include "TSpline.h"
#include "TTree.h"
//#include "TMap.h"
#else
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooAbsData.h"
#include "RooFitResult.h"
#include "RooRealVar.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "RooRandom.h"
#include "RooArgSet.h"
#include "RooRealConstant.h"
#endif
#include <string>
class RooSSCSPdf : public RooAbsPdf {
    public:
        RooSSCSPdf(const char *name, const char *title,
                RooAbsReal& _p11,
                RooAbsReal& _p12,
                RooAbsReal& _p13,
                RooAbsReal& _p14,
                RooAbsReal& _p21,
                RooAbsReal& _p22,
                RooAbsReal& _p23,
                RooAbsReal& _p24,
                RooAbsReal& _p31,
                RooAbsReal& _p32,
                RooAbsReal& _p33,
                RooAbsReal& _p34,
                RooAbsReal& _p41,
                RooAbsReal& _p42,
                RooAbsReal& _p43,
                RooAbsReal& _p44,
                RooAbsReal& _p51,
                RooAbsReal& _p52,
                RooAbsReal& _p53,
                RooAbsReal& _p54,
                RooArgList &params,
                const TString &PHSPDat);
        RooSSCSPdf(const RooSSCSPdf& other, const char* name=0);
        virtual TObject* clone(const char* newname) const { return new RooSSCSPdf(*this,newname);}
        inline virtual ~RooSSCSPdf();
        void project(const char* fname);
        void setPHSPDat(const TString &dat);


        void test();
        Double_t calEva(const Double_t *p1, const Double_t *p2, const
                Double_t *p3, const Double_t*p4, const Double_t *p5) const;
    protected:
        RooRealProxy p11, p12, p13, p14;
        RooRealProxy p21, p22, p23, p24;
        RooRealProxy p31, p32, p33, p34;
        RooRealProxy p41, p42, p43, p44;
        RooRealProxy p51, p52, p53, p54;
        Double_t evaluate() const;

        Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const;
        Double_t analyticalIntegral (Int_t code, const char* rangeName) const;
    private:
        //======= sad ... =====
        RooListProxy _ParameterCol;
        //RooArgList &parameters;
        TIterator *_parsItr;
        void initialize();
        TString  _PHSPDat;

        Double_t **_mcp1;
        Double_t **_mcp2;
        Double_t **_mcp3;
        Double_t **_mcp4;
        Double_t **_mcp5;
        Double_t *_weight;
        Int_t Nmc;
        Double_t bb(const Double_t& alphaLambda, const Double_t
                p4Proton[4], const Double_t p4Sigma[4]) const;
        Double_t bbBar(const Double_t& alphaLambdaBar, const Double_t
                p4Protonbar[4], const Double_t p4Sigmabar[4]) const;
        ClassDef(RooSSCSPdf,1)
};
#endif

