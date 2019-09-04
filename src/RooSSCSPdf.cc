// Copyright 2019 Ma Xinxin
#include "RooSSCSPdf.hh"
#include "TSystem.h"
#include <fstream>
#include "BaseFunc.hh"
#include "RooArgList.h"
#include "TTree.h"
#include "TRandom.h"
#include "TGenPhaseSpace.h"
const Double_t mLambda = 1.115683;
const Double_t mSigma = 1.192642;
const Double_t mProton = 0.938272081;
const Double_t mPion = 0.13957061;
using Func::dot;
using Func::det;
using std::cout;
using std::endl;


ClassImp(RooSSCSPdf)
    RooSSCSPdf::RooSSCSPdf(const char *name, const char *title,
            RooAbsReal& _p11, RooAbsReal& _p12, RooAbsReal& _p13,
            RooAbsReal& _p14, RooAbsReal& _p21, RooAbsReal& _p22,
            RooAbsReal& _p23, RooAbsReal& _p24, RooAbsReal& _p31,
            RooAbsReal& _p32, RooAbsReal& _p33, RooAbsReal& _p34,
            RooAbsReal& _p41, RooAbsReal& _p42, RooAbsReal& _p43,
            RooAbsReal& _p44, RooAbsReal& _p51, RooAbsReal& _p52,
            RooAbsReal& _p53, RooAbsReal& _p54,
            RooArgList & parameters,
            const TString& PHSPdat):
        RooAbsPdf(name, title),
        p11("p11", "p11", this, _p11),
        p12("p12", "p12", this, _p12),
        p13("p13", "p13", this, _p13),
        p14("p14", "p14", this, _p14),
        p21("p21", "p21", this, _p21),
        p22("p22", "p22", this, _p22),
        p23("p23", "p23", this, _p23),
        p24("p24", "p24", this, _p24),
        p31("p31", "p31", this, _p31),
        p32("p32", "p32", this, _p32),
        p33("p33", "p33", this, _p33),
        p34("p34", "p34", this, _p34),
        p41("p41", "p41", this, _p41),
        p42("p42", "p42", this, _p42),
        p43("p43", "p43", this, _p43),
        p44("p44", "p44", this, _p44),
        p51("p51", "p51", this, _p51),
        p52("p52", "p52", this, _p52),
        p53("p53", "p53", this, _p53),
        p54("p54", "p54", this, _p54),
        _ParameterCol("ParameterCol", "", this) {
    _ParameterCol.add(parameters);
    _parsItr = _ParameterCol.createIterator();
    // cout << __LINE__ << endl;
    _PHSPDat = PHSPdat;
    // cout << __LINE__ << endl;
    initialize();
    // cout << __LINE__ << endl;
}


RooSSCSPdf::RooSSCSPdf(const RooSSCSPdf& other, const char* name):
    RooAbsPdf(other, name),
    p11("p11", this, other.p11),
    p12("p12", this, other.p12),
    p13("p13", this, other.p13),
    p14("p14", this, other.p14),
    p21("p21", this, other.p21),
    p22("p22", this, other.p22),
    p23("p23", this, other.p23),
    p24("p24", this, other.p24),
    p31("p31", this, other.p31),
    p32("p32", this, other.p32),
    p33("p33", this, other.p33),
    p34("p34", this, other.p34),
    p41("p41", this, other.p41),
    p42("p42", this, other.p42),
    p43("p43", this, other.p43),
    p44("p44", this, other.p44),
    p51("p51", this, other.p51),
    p52("p52", this, other.p52),
    p53("p53", this, other.p53),
    p54("p54", this, other.p54),
    _ParameterCol("ParameterCol", this, other._ParameterCol) {
    //   cout << __LINE__ << endl;
    _parsItr = _ParameterCol.createIterator();
    Nmc = other.Nmc;
    //  cout << "Nmc:" << Nmc << endl;
    // cout << __LINE__ << endl;
    // _mcp1 = new Double_t*[Nmc];
    //   _mcp2 = new Double_t*[Nmc];
    //   _mcp3 = new Double_t*[Nmc];
    //   _mcp4 = new Double_t*[Nmc];
    //   _mcp4 = new Double_t*[Nmc];
    // _weight = new Double_t[Nmc];
    //  cout << __LINE__ << endl;
    //  for(Int_t i=0;i<Nmc;i++){
    //      cout << __LINE__ << endl;
    //      _mcp1[i] = new Double_t[4];
    //      _mcp2[i] = new Double_t[4];
    //      _mcp3[i] = new Double_t[4];
    //      _mcp4[i] = new Double_t[4];
    //      _mcp5[i] = new Double_t[4];
    //      cout << __LINE__ << endl;
    //  }
    //  cout << __LINE__ << endl;
    _mcp1 = other._mcp1;
    _mcp2 = other._mcp2;
    _mcp3 = other._mcp3;
    _mcp4 = other._mcp4;
    _mcp5 = other._mcp5;
    _weight = other._weight;
    //   cout << __LINE__ << endl;
    //   for(Int_t i=0;i<Nmc;i++){
    //       _mcp1[i][0] = other._mcp1[i][0];
    //       _mcp1[i][1] = other._mcp1[i][1];
    //       _mcp1[i][2] = other._mcp1[i][2];
    //       _mcp1[i][3] = other._mcp1[i][3];
    //       _mcp2[i][0] = other._mcp2[i][0];
    //       _mcp2[i][1] = other._mcp2[i][1];
    //       _mcp2[i][2] = other._mcp2[i][2];
    //       _mcp2[i][3] = other._mcp2[i][3];
    //       _mcp3[i][0] = other._mcp3[i][0];
    //       _mcp3[i][1] = other._mcp3[i][1];
    //       _mcp3[i][2] = other._mcp3[i][2];
    //       _mcp3[i][3] = other._mcp3[i][3];
    //       _weight[i]  = other._weight[i];
    //   }
}


void RooSSCSPdf::initialize() {
    Nmc = 1E6;
    //    cout << __LINE__ << endl;
    _mcp1   = new Double_t*[Nmc];
    _mcp2   = new Double_t*[Nmc];
    _mcp3   = new Double_t*[Nmc];
    _mcp4   = new Double_t*[Nmc];
    _mcp5   = new Double_t*[Nmc];
    _weight = new Double_t[Nmc];
    for (Int_t i = 0; i < Nmc; i++) {
        _mcp1[i] = new Double_t[4];
        _mcp2[i] = new Double_t[4];
        _mcp3[i] = new Double_t[4];
        _mcp4[i] = new Double_t[4];
        _mcp5[i] = new Double_t[4];
    }
    Double_t fx1, fy1, fz1, ft1,
             fx2, fy2, fz2, ft2,
             fx3, fy3, fz3, ft3,
             fx4, fy4, fz4, ft4,
             fx5, fy5, fz5, ft5,
             weight;
    FILE *fp;
    if ((fp = fopen(_PHSPDat, "r")) == NULL) {
        printf("can't open input file");
        exit(0);
        // cout << __LINE__ << endl;
    }
    Int_t i = 0;
    // cout << __LINE__ << endl;
    while (fscanf(fp, "%lf%lf%lf%lf\n%lf%lf%lf%lf\n%lf%lf%lf%lf\n%lf%lf%lf%lf\n%lf%lf%lf%lf%lf\n", // NOLINT
                &ft1, &fx1, &fy1, &fz1,
                &ft2, &fx2, &fy2, &fz2,
                &ft3, &fx3, &fy3, &fz3,
                &ft4, &fx4, &fy4, &fz4,
                &ft5, &fx5, &fy5, &fz5,
                &weight) != EOF) {
        if ( i >= Nmc) break;
        _mcp1[i][0] = ft1;
        _mcp1[i][1] = fx1;
        _mcp1[i][2] = fy1;
        _mcp1[i][3] = fz1;
        _mcp2[i][0] = ft2;
        _mcp2[i][1] = fx2;
        _mcp2[i][2] = fy2;
        _mcp2[i][3] = fz2;
        _mcp3[i][0] = ft3;
        _mcp3[i][1] = fx3;
        _mcp3[i][2] = fy3;
        _mcp3[i][3] = fz3;
        _mcp4[i][0] = ft4;
        _mcp4[i][1] = fx4;
        _mcp4[i][2] = fy4;
        _mcp4[i][3] = fz4;
        _mcp5[i][0] = ft5;
        _mcp5[i][1] = fx5;
        _mcp5[i][2] = fy5;
        _mcp5[i][3] = fz5;
        _weight[i] = weight;
        i++;
    }
    fclose(fp);
    Nmc = i;
    // cout << __LINE__ << endl;
}


Double_t RooSSCSPdf::evaluate() const {
    Double_t p4Sigma[4] = {p11, p12, p13, p14};
    Double_t p4Lambda[4] = {p21, p22, p23, p24};
    Double_t p4Proton[4] = {p31, p32, p33, p34};
    Double_t p4LambdaBar[4] = {p41, p42, p43, p44};
    Double_t p4Protonbar[4] = {p51, p52, p53, p54};
    Double_t pdf = calEva(p4Sigma, p4Lambda, p4Proton,
            p4LambdaBar, p4Protonbar);
    if (pdf < 0) {
        //       cout << "pdf is smaller than 0: " << pdf << endl;
        return 1e-30;
    }
    //   else{
    //
    //   }
    return pdf;
}


Int_t RooSSCSPdf::getAnalyticalIntegral(RooArgSet& allVars,
        RooArgSet& analVars, const char* rangeName) const {
    RooArgSet theSet1, theSet2, theSet3;
    theSet1.add(RooArgSet(p11.arg(), p12.arg(), p13.arg(), p14.arg(),
                p21.arg(), p22.arg(), p23.arg(), p24.arg() ));
    theSet2.add(RooArgSet(p31.arg(), p32.arg(), p33.arg(), p34.arg(),
                p41.arg(), p42.arg(), p43.arg(), p44.arg()));
    theSet3.add(RooArgSet(p51.arg(), p52.arg(), p53.arg(), p54.arg()));
    RooArgSet tmp(theSet1, theSet2, " ");
    RooArgSet theSet(tmp, theSet3, " ");
    if (matchArgs(allVars, analVars, theSet)) {
        return 1;
    }
    return 1;
}


Double_t RooSSCSPdf::analyticalIntegral(Int_t code,
        const char* rangeName) const {
    assert(code == 1);
    Double_t sum = 0;
    // cout << "sum :"  << sum/Nmc << endl;
    for (Int_t i = 0; i < Nmc; i++) {
        Double_t eva = calEva(_mcp1[i], _mcp2[i], _mcp3[i],
                _mcp4[i], _mcp5[i]);
        sum = sum + eva*_weight[i];
    }
    // cout << "sum :"  << sum/Nmc << endl;

    return sum / Nmc;
}


void RooSSCSPdf::project(const char* fname) {
    initialize();
    TFile f(fname, "recreate");

    TTree t("project", "the project of fit result");

    Double_t eva, weight;
    Double_t m_p4_1[4], m_p4_2[4], m_p4_3[4], m_p4_4[4], m_p4_5[4];
    t.Branch("eva", &eva, "eva/D");
    t.Branch("weight", &weight, "weight/D");
    t.Branch("p4_1", m_p4_1, "p4_1[4]/D");
    t.Branch("p4_2", m_p4_2, "p4_2[4]/D");
    t.Branch("p4_3", m_p4_3, "p4_3[4]/D");
    t.Branch("p4_4", m_p4_4, "p4_4[4]/D");
    t.Branch("p4_5", m_p4_5, "p4_5[4]/D");

    for (Int_t i = 0; i < Nmc; i++) {
        // cout << "weight " << _weight[i] << " ";
        // cout << __LINE__ << endl;
        weight = _weight[i];
        //    cout << "weight " << _weight[i] << " ";
        eva = calEva(_mcp1[i], _mcp2[i], _mcp3[i], _mcp4[i],
                _mcp5[i]);
        // cout << _mcp1[i][0] << " " << _mcp1[i][1] << " "
        // << _mcp1[i][2] << " " << _mcp1[i][3] << " " << endl;
        // cout << _mcp2[i][0] << " " << _mcp2[i][1] << " "
        // << _mcp2[i][2] << " " << _mcp2[i][3] << " " << endl;
        // cout << _mcp3[i][0] << " " << _mcp3[i][1] << " "
        // << _mcp3[i][2] << " " << _mcp3[i][3] << " " << endl;
        // cout << _mcp4[i][0] << " " << _mcp4[i][1] << " "
        // << _mcp4[i][2] << " " << _mcp4[i][3] << " " << endl;
        // cout << _mcp5[i][0] << " " << _mcp5[i][1] << " "
        //  // << _mcp5[i][2] << " " << _mcp5[i][3] << " " << endl;
            m_p4_1[3] = _mcp1[i][0];
            m_p4_1[0] = _mcp1[i][1];
            m_p4_2[3] = _mcp2[i][0];
            m_p4_2[0] = _mcp2[i][1];
            m_p4_3[3] = _mcp3[i][0];
            m_p4_3[0] = _mcp3[i][1];
            m_p4_4[3] = _mcp4[i][0];
            m_p4_4[0] = _mcp4[i][1];
            m_p4_5[3] = _mcp5[i][0];
            m_p4_5[0] = _mcp5[i][1];
            m_p4_1[1] = _mcp1[i][2];
            m_p4_1[2] = _mcp1[i][3];
            m_p4_2[1] = _mcp2[i][2];
            m_p4_2[2] = _mcp2[i][3];
            m_p4_3[1] = _mcp3[i][2];
            m_p4_3[2] = _mcp3[i][3];
            m_p4_4[1] = _mcp4[i][2];
            m_p4_4[2] = _mcp4[i][3];
            m_p4_5[1] = _mcp5[i][2];
            m_p4_5[2] = _mcp5[i][3];
        t.Fill();
    }
    t.Write();
    f.Close();
}

void RooSSCSPdf::test() {
    for (Int_t i = 0; i < Nmc; i++) {
        cout << _mcp1[i][0] << " " << _mcp1[i][1] << " " << _mcp1[i][2]
            << " " << _mcp1[i][3] << " " << endl;
        cout << _mcp2[i][0] << " " << _mcp2[i][1] << " " << _mcp2[i][2]
            << " " << _mcp2[i][3] << " " << endl;
        cout << _mcp3[i][0] << " " << _mcp3[i][1] << " " << _mcp3[i][2]
            << " " << _mcp3[i][3] << " " << endl;
        cout << _mcp4[i][0] << " " << _mcp4[i][1] << " " << _mcp4[i][2]
            << " " << _mcp4[i][3] << " " << endl;
        cout << _mcp5[i][0] << " " << _mcp5[i][1] << " " << _mcp5[i][2]
            << " " << _mcp5[i][3] << " " << endl;
        Double_t pdf = calEva(_mcp1[i], _mcp2[i], _mcp3[i],
                _mcp4[i], _mcp5[i]);
        if (pdf < 0) {
            cout << "pdf is less than 0: " << pdf << endl;
            cout << "p4Sigma ==> "
                 << _mcp1[i][0] << ", "
                 << _mcp1[i][1] << ", "
                 << _mcp1[i][2] << ", "
                 << _mcp1[i][3] << " "
                 << sqrt(dot(_mcp1[i], _mcp1[i])) << " " << endl;

            cout << "p4Lambda ==> "
                 << _mcp2[i][0] << ", "
                 << _mcp2[i][1] << ", "
                 << _mcp2[i][2] << ", "
                 << _mcp2[i][3] << " "
                 << sqrt(dot(_mcp2[i], _mcp2[i]) ) << " " << endl;

            cout << "p4Proton ==> "
                 << _mcp3[i][0] << ", "
                 << _mcp3[i][1] << ", "
                 << _mcp3[i][2] << ", "
                 << _mcp3[i][3] << " "
                 << sqrt(dot(_mcp3[i], _mcp3[i]) ) << " " << endl;


            cout << "p4LambdBar ==> "
                 << _mcp4[i][0] << ", "
                 << _mcp4[i][1] << ", "
                 << _mcp4[i][2] << ", "
                 << _mcp4[i][3] << " "
                 << sqrt(dot(_mcp4[i], _mcp4[i]) ) << " " << endl;

            cout << "p4ProtonBar ==> "
                 << _mcp5[i][0] << ", "
                 << _mcp5[i][1] << ", "
                 << _mcp5[i][2] << ", "
                 << _mcp5[i][3] << " "
                 << sqrt(dot(_mcp5[i], _mcp5[i])) << " " << endl;
        }
    }
}


void RooSSCSPdf::setPHSPDat(const TString &dat) {
    _PHSPDat = dat;
}


inline RooSSCSPdf::~RooSSCSPdf() {
    for (Int_t i = 0; i < Nmc; i++) {
        delete[] _mcp1[i];
        delete[] _mcp2[i];
        delete[] _mcp3[i];
        delete[] _mcp4[i];
        delete[] _mcp5[i];
    }
    delete []_weight;
}


Double_t RooSSCSPdf::bb(const Double_t& alphaLambda,
        const Double_t p4Proton[4], const Double_t p4Sigma[4]) const {
    return (-2*alphaLambda*mSigma*((pow(mLambda, 2) - pow(mPion, 2)
        + pow(mProton, 2))*(pow(mLambda, 2) + pow(mSigma, 2))
        - 4*pow(mLambda, 2)*dot(p4Proton, p4Sigma)))/
        (mLambda*sqrt(pow(mLambda, 2) + pow(pow(mPion, 2)
        - pow(mProton, 2), 2)/pow(mLambda, 2) - 2*(pow(mPion, 2)
        + pow(mProton, 2)))*pow(pow(mLambda, 2) - pow(mSigma, 2), 2));
}


Double_t RooSSCSPdf::bbBar(const Double_t& alphaLambdaBar,
        const Double_t p4Protonbar[4], const Double_t p4Sigmabar[4]) const {
    return (-2*alphaLambdaBar*mSigma*((pow(mLambda, 2) - pow(mPion, 2)
        + pow(mProton, 2))*(pow(mLambda, 2) + pow(mSigma, 2))
        - 4*pow(mLambda, 2)*dot(p4Protonbar, p4Sigmabar))) /
        (mLambda*sqrt(pow(mLambda, 2) + pow(pow(mPion, 2)
        - pow(mProton, 2), 2) / pow(mLambda, 2) - 2*(pow(mPion, 2)
        + pow(mProton, 2)))*pow(pow(mLambda, 2) - pow(mSigma, 2), 2));
}


Double_t RooSSCSPdf::calEva(const Double_t p4Sigma[4],
        const Double_t p4Lambda[4],  const Double_t p4Proton[4],
        const Double_t p4Lambdabar[4], const Double_t p4Protonbar[4]) const {
    _parsItr->Reset();
    // R, phi, alphaL, alphaLbar
    RooRealVar *aPara(0);
    aPara = reinterpret_cast<RooRealVar*>(_parsItr->Next());
    Double_t R = aPara->getVal();

    aPara = reinterpret_cast<RooRealVar*>(_parsItr->Next());
    Double_t phi = aPara->getVal();

    aPara = reinterpret_cast<RooRealVar*>(_parsItr->Next());
    Double_t alphaL = aPara->getVal();

    aPara = reinterpret_cast<RooRealVar*>(_parsItr->Next());
    Double_t alphaLbar = aPara->getVal();

    const Double_t k2[4]= {3.097/2, 0, 0, 3.097/2};
    const Double_t k1[4]= {3.097/2, 0, 0, -3.097/2};
    const Double_t p4Q[4]= {3.097, 0, 0, 0};
    const Double_t tau = dot(p4Q, p4Q) / (4*mSigma*mSigma);
    const Double_t s = pow(3.097, 2);
    // const Double_t b = bb(alphaL, p4Proton, p4Sigma);
    // const Double_t bBar = bbBar(alphaLbar, p4Protonbar, p4Sigmabar);
    Double_t alphapsi = R;
    Double_t alphaA = alphaL;
    Double_t alphaB = alphaLbar;
    Double_t deltaPhi = -phi;

    TLorentzVector P4S0(p4Sigma[1], p4Sigma[2], p4Sigma[3], p4Sigma[0]);
    TLorentzVector P4S0Bar(-p4Sigma[1], -p4Sigma[2], -p4Sigma[3],
            3.097 -  p4Sigma[0]);

    TLorentzVector P4Lambda(p4Lambda[1], p4Lambda[2],
            p4Lambda[3], p4Lambda[0]);

    TLorentzVector P4LambdaBar(p4Lambdabar[1], p4Lambdabar[2],
            p4Lambdabar[3], p4Lambdabar[0]);

    TLorentzVector P4Proton(p4Proton[1], p4Proton[2], p4Proton[3],
            p4Proton[0]);
    TLorentzVector P4ProtonBar(p4Protonbar[1], p4Protonbar[2],
            p4Protonbar[3], p4Protonbar[0]);

    // boost Lambda and proton into the rest frame of Sigma0
    P4Lambda.Boost(-P4S0.BoostVector());
    P4Lambda.RotateZ(-P4S0.Phi());
    P4Lambda.RotateY(-P4S0.Theta());

    P4Proton.Boost(-P4S0.BoostVector());
    P4Proton.RotateZ(-P4S0.Phi());
    P4Proton.RotateY(-P4S0.Theta());

    // boost the proton into the frame of Lambda
    P4Proton.Boost(-P4Lambda.BoostVector());
    P4Proton.RotateZ(-P4Lambda.Phi());
    P4Proton.RotateY(-P4Lambda.Theta());

    // boost anti-Lambda and anti-proton into the rest frame of anti-Sigma0
    P4LambdaBar.Boost(-P4S0Bar.BoostVector());
    P4LambdaBar.RotateZ(-P4S0Bar.Phi());
    P4LambdaBar.RotateY(-P4S0Bar.Theta());

    P4ProtonBar.Boost(-P4S0Bar.BoostVector());
    P4ProtonBar.RotateZ(-P4S0Bar.Phi());
    P4ProtonBar.RotateY(-P4S0Bar.Theta());

    // boost the anti-proton into the frame of anti-Lambda
    P4ProtonBar.Boost(-P4LambdaBar.BoostVector());
    P4ProtonBar.RotateZ(-P4LambdaBar.Phi());
    P4ProtonBar.RotateY(-P4LambdaBar.Theta());

    // give the angles
    Double_t theta = P4S0.Theta();
    Double_t theta1 = P4Lambda.Theta();
    Double_t phi1 = P4Lambda.Phi();
    Double_t theta2 = P4LambdaBar.Theta();
    Double_t phi2 = P4LambdaBar.Phi();
    Double_t thetaA = P4Proton.Theta();
    Double_t thetaB = P4ProtonBar.Theta();

    // Double_t bb
    // cout << "alphaL = " << alphaL << endl;
    // cout << "alphaLbar = " << alphaLbar << endl;
    // cout << "b = " << b << endl;
    // cout << "bBar = " << bBar << endl;
    return 1 + alphapsi*pow(cos(theta), 2) + alphaA*alphaB*( - alphapsi
            - pow(cos(theta), 2)) * cos(theta1) * cos(theta2)
        * cos(thetaA)*cos(thetaB) + alphaA*alphaB*sqrt(1 - pow(alphapsi, 2))
        * cos(deltaPhi)*cos(phi1)*cos(theta)*cos(theta2)*cos(thetaA)
        * cos(thetaB)*sin(theta)*sin(theta1) + alphaA
        * sqrt(1 - pow(alphapsi, 2))*cos(theta)*cos(thetaA)*sin(deltaPhi)
        * sin(phi1)*sin(theta)*sin(theta1) - alphaA*alphaB
        * sqrt(1 - pow(alphapsi, 2))*cos(deltaPhi)*cos(phi2)*cos(theta)
        * cos(theta1)*cos(thetaA)*cos(thetaB)*sin(theta)
        * sin(theta2) - alphaB*sqrt(1 - pow(alphapsi, 2))*cos(theta)
        * cos(thetaB)*sin(deltaPhi)*sin(phi2)*sin(theta)*sin(theta2)
        + alphaA*alphaB*cos(phi1)*cos(phi2)*cos(thetaA)*cos(thetaB)
        * pow(sin(theta), 2)*sin(theta1)*sin(theta2) + alphaA*alphaB
        * alphapsi*cos(thetaA)*cos(thetaB)*sin(phi1)*sin(phi2)
        * pow(sin(theta), 2)*sin(theta1)*sin(theta2);
}
