/*!  \file proton .C
     \brief plot proton spectrum used for IRF generation

     usage: plot_proton_spectrum();

     \author Gernot Maier
*/
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TMath.h"
#include "TSystem.h"

/*

   fit functions

   Note: all functions are multiplied by E^3

*/

double powerlaw( double* x, double* p )
{
	double e = x[0];    // energy in [TeV]
	double flux = p[0];
	double gamma = p[1];
	
	double f = flux;
	f       *= TMath::Power( e / 1., -gamma );
	
	// all values are multiplied by E^2.7
	f *= TMath::Power( e, 2.7 );
	
	return f;
}


double powerlaw_exp( double* x, double* p )
{
	double e = x[0];    // energy in [TeV]
	double flux = p[0];
	double gamma = p[1];
	double cutoff = p[2];
	double ep = p[3];
	
	double f = flux;
	f       *= TMath::Power( e / ep, -gamma );
	f       *= TMath::Exp( -1. * e / cutoff );
	
	// all values are multiplied by E^2.7
	f *= TMath::Power( e, 2.7 );
	
	return f;
}

/*
   all energies in [TeV]
*/
double powerlaw_broken( double* x, double* p )
{
	double e = x[0];
	double flux = p[0];
	double e_p = p[1];
	double gamma_1 = p[2];
	double gamma_2 = p[3];
	double alpha = p[4];
	
	double f = flux;
	f       *= TMath::Power( e / e_p, -1.* gamma_1 );
	f       *= TMath::Power( 1. + TMath::Power( e / e_p, 1. / alpha ), -1.*( gamma_2 - gamma_1 ) * alpha );
	
	// all values are multiplied by E^2.7
	f *= TMath::Power( e, 2.7 );
	
	return f;
}

double powerlaw_plus_lognormal( double* x, double* p )
{
	double e = x[0];
	double flux = p[0];
	double e_p = p[1];
	double gamma = p[2];
	double e_peak = p[3];
	double amp_peak = p[4];
	double width_peak = p[5];
	
	double f = flux;
	f       *= TMath::Power( e / e_p, -gamma );
	f       += amp_peak / (e*width_peak*sqrt(2.*TMath::Pi() ) ) * TMath::Gaus( log( e ), log( e_peak ), width_peak, false );

	// all values are multiplied by E^2.7
	f *= TMath::Power( e, 2.7 );
	
	return f;
}


/*

    plot spectra

    double iEScale :   scale energies from LAT data by this factor

*/
void plot_proton_spectrum()
{
    // proton spectrum from DAMPE
    // https://arxiv.org/pdf/1909.12860.pdf
    // data from Table 1 (statistical errors only)
    TGraphErrors *fDampe = new TGraphErrors( "proton_data/DAMPE_2019.dat", "%lg %lg %lg" );
    fDampe->SetMarkerStyle( 20 );
    fDampe->SetMarkerColor( 801 );
    fDampe->SetLineColor( 801 );

    double x, y, yerr;
    for( int i = 0; i < fDampe->GetN(); i++ )
    {
        fDampe->GetPoint( i, x, y );
        yerr = fDampe->GetErrorY( i );

        // GeV --> TeV
        x /= 1.e3;

        // 1/GeV/s/sr --> 1/TeV/s/sr
        y *= 1.e3;
        yerr *= 1.e3;

        // all values are multiplied by E^2.7
        y *= TMath::Power( x, 2.7 );
        yerr *= TMath::Power( x, 2.7 );

        fDampe->SetPoint( i, x, y );
        fDampe->SetPointError( i, 0., yerr );
    }

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// plotting
	TCanvas* c = new TCanvas( "cP", "proton spectrum", 10, 10, 900, 600 );
	c->SetLogx( 1 );
	c->SetLeftMargin( 0.15 );
	c->Draw();
	
	TH1D* h = new TH1D( "hnull", "", 10, 5. / 1.e2, 2.e2 );
	h->SetMinimum( 0.05 );
	h->SetMaximum( 0.15 );
	h->SetXTitle( "energy [TeV]" );
	h->SetYTitle( "E^{2.7} J(E) [TeV^{1.7}m^{-2}s^{-1}sr^{-1}]" );
	h->GetYaxis()->SetTitleOffset( 1.6 );
	h->SetStats( 0 );
	h->Draw();
	
	// ATIC fit to proton spectrum (S.Swordy, private communications)
	TF1* fPowerLaw_IRF = new TF1( "PowerLaw_IRF", powerlaw, 0.05, 1.e6, 2 );
	fPowerLaw_IRF->SetLineColor( 1 );
	//   fPowerLaw_IRF->SetLineStyle( 2 );
    // 1/cm2 --> 1/m2
	fPowerLaw_IRF->SetParameter( 0, 0.098e-4 * 1.e4 );
	fPowerLaw_IRF->SetParameter( 1, 2.62 );
	
	fPowerLaw_IRF->Draw( "same" );
    
	// MARS fit to proton spectrum (BESS from  ApJ 545, 1135 (2000))
	TF1* fPowerLaw_MARS = new TF1( "PowerLaw_MARS", powerlaw, 0.05, 1.e6, 2 );
	fPowerLaw_MARS->SetLineColor( 4 );
	//   fPowerLaw_MARS->SetLineStyle( 2 );
    // 1/cm2 --> 1/m2
	fPowerLaw_MARS->SetParameter( 0, 0.0969e-4 * 1.e4 );
	fPowerLaw_MARS->SetParameter( 1, 2.7 );
	
	fPowerLaw_MARS->Draw( "same" );


    fDampe->Draw( "p" );
	
}

