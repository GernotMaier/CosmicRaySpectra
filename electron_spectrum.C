/*!  \file electron_spectrum.C
     \brief plot electron spectrum from LAT and HESS

     usage: plot_electron_spectrum();

     Note: an energy independent systematic error is added to the LAT data

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
	
	// all values are multiplied by E^3
	
	f *= TMath::Power( e, 3. );
	
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
	
	// all values are multiplied by E^3
	
	f *= TMath::Power( e, 3. );
	
	return f;
}

/*
   all energies in [TeV]
*/
double powerlaw_broken_simple( double* x, double* p )
{
	double e = x[0];
	double flux = p[0];
	double e_p = p[1];
	double gamma_1 = p[2];
	double gamma_2 = p[3];
	
	double f = flux;
	if( e < e_p )
	{
		f *= TMath::Power( e / e_p, -1.* gamma_1 );
	}
	else
	{
		f *= TMath::Power( e / e_p, -1.* gamma_2 );
	}
	
	// all values are multiplied by E^3
	
	f *= TMath::Power( e, 3. );
	
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
	
	// all values are multiplied by E^3
	
	f *= TMath::Power( e, 3. );
	
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
	f       += amp_peak / ( e * width_peak * sqrt( 2.*TMath::Pi() ) ) * TMath::Gaus( log( e ), log( e_peak ), width_peak, false );
	
	// all values are multiplied by E^3
	f *= TMath::Power( e, 3. );
	
	return f;
}


/*

    plot spectra

    double iEScale :   scale energies from LAT data by this factor

*/
void plot_electron_spectrum( double iEScale = 1., bool bPlotTotal = false, bool bPlotFits = true )
{

	// read source data
	vector< string > inputfile;
	vector< string > legendName;
	inputfile.push_back( "./electron_data/HESS_AA_508_261_2009.dat" );
	legendName.push_back( "HESS A&A 508, 262 (2009)" );
	inputfile.push_back( "./electron_data/FERMI_LAT_PRL_102_181101_2009.dat" );
	legendName.push_back( "LAT PRL 102 181101 (2009)" );
	inputfile.push_back( "./electron_data/FERMI_LAT_LE_PRD_82_092004_2010.dat" );
	legendName.push_back( "LAT (LE) PRD 82, 092004 (2010)" );
	inputfile.push_back( "./electron_data/FERMI_LAT_HE_PRD_82_092004_2010.dat" );
	legendName.push_back( "LAT (HE) PRD 82, 092004 (2010)" );
	//        inputfile.push_back( "./electron_data/AMS_PRD_113_121102_2014.dat" );
	//        legendName.push_back( "AMS PRD 113, 121102 (2014)" );
	inputfile.push_back( "./electron_data/AMS_PRD_113_221102_2014.dat" );
	legendName.push_back( "AMS PRD 113, 221102 (2014)" );
	inputfile.push_back( "./electron_data/VERITAS_ICRC2015.dat" );
	legendName.push_back( "VERITAS PRD D 98, 062004 (2018)" );
	inputfile.push_back( "./electron_data/VER-300000-sedCR-1.ecsv" );
	legendName.push_back( "VTSCat" );
	inputfile.push_back( "./electron_data/FERMI_LAT_PRD_95_082007_2017.dat" );
	legendName.push_back( "LAT PRD 95, 082007 (2017)" );
	inputfile.push_back( "./electron_data/HESS_ICRC2017.dat" );
	legendName.push_back( "HESS ICRC 2017" );
	
	vector< TGraphAsymmErrors* > fSpectrum;
	TGraphAsymmErrors* fSpectrum_total = new TGraphAsymmErrors( 1 );
	fSpectrum_total->SetMarkerStyle( 24 );
	fSpectrum_total->SetMarkerColor( 8 );
	fSpectrum_total->SetLineColor( 8 );
	int z = 0;
	
	TGraphAsymmErrors* fSpectrum_syst = new TGraphAsymmErrors( 1 );
	fSpectrum_syst->SetMarkerStyle( 3 );
	fSpectrum_syst->SetMarkerColor( 4 );
	fSpectrum_syst->SetLineColor( 4 );
	fSpectrum_syst->SetLineStyle( 2 );
	int z_sys = 0;
	
	for( unsigned int i = 0; i < inputfile.size(); i++ )
	{
		fSpectrum.push_back( new TGraphAsymmErrors( 1 ) );
		
		ifstream is( inputfile[i].c_str() );
		if( !is )
		{
			cout << "Error: data file " << inputfile[i] << " not found" << endl;
			continue;
		}
		cout << "reading " << inputfile[i] << endl;
		
		int n = 0;
		double xmean = 0.;
		double x_low = 0.;
		double x_up = 0.;
		double y = 0.;
		double y_E_up = 0.;
		double y_E_low = 0.;
		double y_E_stat_low = 0.;
		double y_E_stat_up = 0.;
		for( ;; )
		{
			// energies in GeV
			if( inputfile[i].find( "HESS_AA_508_261_2009.dat" ) < inputfile[i].size() )
			{
				is >> xmean;
				x_low = xmean;
				x_up = xmean;
			}
			// energies in TeV
			else if( inputfile[i].find( "HESS_ICRC2017.dat" ) < inputfile[i].size()
					 || inputfile[i].find( "VERITAS" ) < inputfile[i].size() )
			{
				is >> xmean;
				is >> x_low;
				is >> x_up;
				if( inputfile[i].find( "HESS_ICRC2017.dat" ) < inputfile[i].size() )
				{
					xmean *= 1.e3;
					x_low *= 1.e3;
					x_up *= 1.e3;
				}
			}
            else if( inputfile[i].find( "VER-" ) < inputfile[i].size() )
            {
				is >> xmean;
				is >> x_low;
				is >> x_up;
                xmean *= 1.e3;
                x_low *= 1.e3;
                x_up *= 1.e3;
            }
			else
			{
				is >> x_low;
				is >> x_up;
				if( inputfile[i].find( "AMS" ) < inputfile[i].size() )
				{
					is >> xmean;
				}
				else
				{
					xmean = 0.5 * ( x_low + x_up );
				}
			}
			
			if( is.eof() )
			{
				break;
			}
			
			// fluxes as GeV^2 m^-2 s^-1 sr^-1
			is >> y;
			
			if( inputfile[i].find( "HESS_AA_508_261_2009.dat" ) < inputfile[i].size() )
			{
				is >> y_E_low;
				is >> y_E_up;
			}
            else if( inputfile[i].find( "VER-" ) < inputfile[i].size() )
            {
                is >> y_E_low;
                y_E_up = y_E_low;
            }
			else
			{
				is >> y_E_low;
				// fractional error given in LAT 2017 paper
				if( inputfile[i].find( "PRD_95_082007_2017" ) < inputfile[i].size() )
				{
					y_E_low *= y;
				}
				y_E_up = y_E_low;
				if( inputfile[i].find( "AMS_PRD_113_121102_2014.dat" ) < inputfile[i].size() )
				{
					is >> y_E_up;
				}
			}
			// last column: order of magnitude column
			if( inputfile[i].find( "AMS_PRD_113_221102_2014.dat" ) < inputfile[i].size() )
			{
				double ymult;
				is >> ymult;
				is >> ymult;
				y *= ymult;
				y_E_low *= ymult;
				y_E_up *= ymult;
			}
			// read systematic fluxes
			y_E_stat_low = -1.;
			y_E_stat_up  = -1.;
			if( inputfile[i].find( "FERMI_LAT_HE_PRD_82_092004_2010" ) < inputfile[i].size() )
			{
				is >> y_E_stat_low;
				is >> y_E_stat_up;
			}
			// Fermi/LAT 2010 data has different units: Multiply by E^3
			if( inputfile[i].find( "PRD_82_092004_2010" ) < inputfile[i].size()
					|| inputfile[i].find( "PRD_95_082007_2017" ) < inputfile[i].size()
					|| inputfile[i].find( "AMS" ) < inputfile[i].size()
                    || inputfile[i].find( "VER-" ) < inputfile[i].size() )
			{
				y       *= TMath::Power( xmean, 3. );
				y_E_low *= TMath::Power( xmean, 3. );
				y_E_up  *= TMath::Power( xmean, 3. );
				y_E_stat_low *= TMath::Power( xmean, 3. );
				y_E_stat_up  *= TMath::Power( xmean, 3. );
			}
			// GeV to TeV
			xmean /= 1.e3;
			x_up  /= 1.e3;
			x_low /= 1.e3;
			// GeV to TeV
			y       /= TMath::Power( 1.e3, 2. );
			y_E_up  /= TMath::Power( 1.e3, 2. );
			y_E_low /= TMath::Power( 1.e3, 2. );
			y_E_stat_low /= TMath::Power( 1.e3, 2. );
			y_E_stat_up  /= TMath::Power( 1.e3, 2. );

            if( inputfile[i].find( "VERITAS" ) < inputfile[i].size() )
            {
                cout << xmean << "  " << x_low << "  " << x_up << " ";
                // energies: TeV --> GeV
                cout << y / TMath::Power( xmean*1.e3, 3. ) * TMath::Power( 1.e3, 2. );
                // fluxes from TeV^2 /m/s/sr to 
                cout << "  " << y_E_up / TMath::Power( xmean*1.e3, 3. ) * TMath::Power( 1.e3, 2. );
                cout << "  " << y_E_low / TMath::Power( xmean*1.e3, 3. ) * TMath::Power( 1.e3, 2. );
                cout << endl;
            }  
			
			fSpectrum[i]->SetPoint( n, xmean, y );
			fSpectrum[i]->SetPointEYhigh( n, y_E_up );
			fSpectrum[i]->SetPointEYlow( n, y_E_low );
			fSpectrum[i]->SetPointEXhigh( n, x_up - xmean );
			fSpectrum[i]->SetPointEXlow( n, xmean - x_low );
			
			if( y_E_stat_low > 0. )
			{
				fSpectrum_syst->SetPoint( z_sys, xmean, y );
				fSpectrum_syst->SetPointEXhigh( z_sys, x_up - xmean );
				fSpectrum_syst->SetPointEXlow( z_sys, xmean - x_low );
				fSpectrum_syst->SetPointEYlow( z_sys, y_E_stat_low );
				fSpectrum_syst->SetPointEYhigh( z_sys, y_E_stat_up );
				z_sys++;
			}
			
			if( i == 0 )
			{
				fSpectrum_total->SetPoint( z, xmean * iEScale, y * TMath::Power( iEScale, 3. ) );
			}
			else
			{
				fSpectrum_total->SetPoint( z, xmean, y );
			}
			if( inputfile[i].find( "HESS" ) < inputfile[i].size() )
			{
				fSpectrum_total->SetPointEYhigh( z, y_E_up );
				fSpectrum_total->SetPointEYlow( z, y_E_low );
			}
			// add an energy independent systematic error for Fermi LAT data (larger for LE selection)
			else if( inputfile[i].find( "LE" ) < inputfile[i].size() )
			{
				fSpectrum_total->SetPointEYhigh( z, y_E_up + 14. / TMath::Power( 1.e3, 2. ) );
				fSpectrum_total->SetPointEYlow( z, y_E_low + 14. / TMath::Power( 1.e3, 2. ) );
			}
			// add an energy independent systematic error for Fermi LAT data
			else
			{
				fSpectrum_total->SetPointEYhigh( z, y_E_up + 7. / TMath::Power( 1.e3, 2. ) );
				fSpectrum_total->SetPointEYlow( z, y_E_low + 7. / TMath::Power( 1.e3, 2. ) );
			}
			fSpectrum_total->SetPointEXhigh( z, x_up - xmean );
			fSpectrum_total->SetPointEXlow( z, xmean - x_low );
			
			z++;
			n++;
		}
		cout << "\t total number of spectral points: " << n << endl;
	}
	
	////////////////////////////////////////////////////////////////////////////////////////////////////
	// plotting
	TCanvas* c = new TCanvas( "cS", "electron spectrum", 10, 10, 900, 600 );
	c->SetLogx( 1 );
	c->SetLeftMargin( 0.15 );
	c->Draw();
	
	TH1D* h = new TH1D( "hnull", "", 10, 5. / 1.e3, 30000. / 1.e3 );
	h->SetMinimum( 0.5e-6 );
	h->SetMaximum( 2.5e-4 );
	h->SetXTitle( "energy [TeV]" );
	h->SetYTitle( "E^{3} J(E) [TeV^{2}m^{-2}s^{-1}sr^{-1}]" );
	h->GetYaxis()->SetTitleOffset( 1.6 );
	h->SetStats( 0 );
	h->Draw();
	
	TLegend* iL = new TLegend( 0.2, 0.13, 0.5, 0.3 );
	
	unsigned int cc = 1;
	for( unsigned int i = 0; i < inputfile.size(); i++ )
	{
		if( fSpectrum[i] )
		{
			fSpectrum[i]->SetMarkerColor( cc );
			fSpectrum[i]->SetLineColor( cc );
			fSpectrum[i]->SetMarkerStyle( 20 + i );
			if( legendName[i].find( "VERITAS" ) != string::npos )
			{
				fSpectrum[i]->SetMarkerStyle( 20 );
				fSpectrum[i]->SetMarkerColor( 800 );
				fSpectrum[i]->SetLineColor( 800 );
			}
            else if( legendName[i].find( "VTSCat" ) != string::npos )
            {
                cout << "Flux at 710 GeV: " << fSpectrum[i]->Eval( 0.710 ) / TMath::Power( 0.710, 3. ) << endl;
            }
			
			fSpectrum[i]->Draw( "p" );
			
			iL->AddEntry( fSpectrum[i], legendName[i].c_str(), "pl" );
			cc++;
			if( cc == 5 )
			{
				cc++;
			}
		}
	}
	if( bPlotTotal )
	{
		fSpectrum_total->Draw( "p" );
	}
	iL->Draw();
	
	TLegend* iLF = new TLegend( 0.5, 0.7, 0.87, 0.87 );
	
	// do not plot fits
	if( !bPlotFits )
	{
		return;
	}
	
	//////////////////////
	// fits
	// (do not plot)
	TF1* fPowerLaw_Exp = new TF1( "PowerLaw_Exp", powerlaw_exp, 5. / 1.e3, 20000. / 1.e3, 4 );
	fPowerLaw_Exp->SetParameter( 0, 1.2e-4 );
	fPowerLaw_Exp->SetParameter( 1, 3. );
	fPowerLaw_Exp->SetParameter( 2, 3. );
	fPowerLaw_Exp->FixParameter( 2, 3.1 );
	fPowerLaw_Exp->SetParameter( 3, 1. );
	fPowerLaw_Exp->FixParameter( 3, 1. );
	
	cout << "Fit: powerlaw with exponential cut off" << endl;
	fSpectrum_total->Fit( "PowerLaw_Exp", "EMR" );
	cout << "======================================" << endl;
	
	fPowerLaw_Exp->SetLineStyle( 2 );
	// fPowerLaw_Exp->Draw( "sames" );
	// iLF->AddEntry( fPowerLaw_Exp, "power law with exp. cut off (global fit)", "l" );
	
	TF1* fPowerLawBroken = new TF1( "PowerLawBroken", powerlaw_broken, 5. / 1.e3, 20000. / 1.e3, 5 );
	fPowerLawBroken->SetLineColor( 2 );
	fPowerLawBroken->SetParameter( 0, 1.2e-4 );
	fPowerLawBroken->SetParameter( 1, 1. );
	fPowerLawBroken->SetParameter( 2, 0.9 );
	fPowerLawBroken->SetParameter( 3, 4. );
	fPowerLawBroken->SetParLimits( 3, 2., 5. );
	fPowerLawBroken->SetParameter( 4, 0.3 );
	fPowerLawBroken->FixParameter( 4, 0.4 );
	
	cout << "Fit: broken power law" << endl;
	fSpectrum_total->Fit( "PowerLawBroken", "EMR" );
	cout << "======================================" << endl;
	
	fPowerLawBroken->SetLineStyle( 2 );
	// fPowerLawBroken->Draw( "sames" );
	// iLF->AddEntry( fPowerLawBroken, " broken power law (global fit)", "l" );
	
	//////////////////////////////////////////////////////////////////////////////////////////////
	// values used in CTA analysis
	
	// old DESY analysis
	TF1* fPowerLawBroken_DESYAnalysis = new TF1( "PowerLawBroken_DESYAnalysis", powerlaw_broken, 7.e-3, 2., 5 );
	fPowerLawBroken_DESYAnalysis->SetLineColor( 5 );
	fPowerLawBroken_DESYAnalysis->SetParameter( 0, 2.3e-5 );
	fPowerLawBroken_DESYAnalysis->SetParameter( 1, 1.7 );
	fPowerLawBroken_DESYAnalysis->SetParameter( 2, 3.07 );
	fPowerLawBroken_DESYAnalysis->SetParameter( 3, 5.0 );
	fPowerLawBroken_DESYAnalysis->SetParameter( 4, 0.4 );
	
	// fPowerLawBroken_DESYAnalysis->Draw( "same" );
	// iLF->AddEntry( fPowerLawBroken_DESYAnalysis, "DESY analysis", "l" );
	
	///////////////////////////////////////////////////////////////////////////////////////////////
	// HESS broken power law from HESS A&A 508, 561 (2009)
	TF1* fPowerLawBroken_HESS = new TF1( "PowerLawBroken_HESS", powerlaw_broken, 350. / 1.e3, 2000. / 1.e3, 5 );
	fPowerLawBroken_HESS->SetLineColor( 1 );
	fPowerLawBroken_HESS->SetParameter( 0, 1.54e-4 );
	fPowerLawBroken_HESS->SetParameter( 1, 1.00 );
	fPowerLawBroken_HESS->SetParameter( 2, 3. );
	fPowerLawBroken_HESS->SetParameter( 3, 4.1 );
	fPowerLawBroken_HESS->SetParameter( 4, 0.2 );
	fPowerLawBroken_HESS->SetLineWidth( 1 );
	fPowerLawBroken_HESS->SetLineStyle( 2 );
	fPowerLawBroken_HESS->Draw( "same" );
	iLF->AddEntry( fPowerLawBroken_HESS, "HESS A&A 508, 561 (2009)", "l" );
	
	// HESS power law from HESS A&A 508, 561 (2009) and PRL 101, 261104 (2008)
	TF1* fPowerLaw_HESS = new TF1( "PowerLaw_HESS", powerlaw, 1000. / 1.e3, 20000. / 1.e3, 2 );
	fPowerLaw_HESS->SetLineColor( 3 );
	//   fPowerLaw_HESS->SetLineStyle( 2 );
	fPowerLaw_HESS->SetParameter( 0, 1.17e-4 );
	fPowerLaw_HESS->SetParameter( 1, 3.9 );
	
	// fPowerLaw_HESS->Draw( "same" );
	// iLF->AddEntry( fPowerLaw_HESS, "HESS PRL 101, 261104 (2008)", "l" );
	
	// Fit B from HESS/ATTIC/Pamela data (PRL 101, 261104 (2008))
	TF1* fPowerLawExp_HESS = new TF1( "PowerLawExp_HESS", powerlaw_exp, 5. / 1.e3, 20000. / 1.e3, 4 );
	fPowerLawExp_HESS->SetLineColor( 6 );
	//   fPowerLawExp_HESS->SetLineStyle( 4 );
	fPowerLawExp_HESS->SetParameter( 0, 1.17e-4 );
	fPowerLawExp_HESS->SetParameter( 1, 3.05 );
	fPowerLawExp_HESS->SetParameter( 2, 2.1 );
	fPowerLawExp_HESS->SetParameter( 3, 1.11 );
	
	//   fPowerLawExp_HESS->Draw( "same" );
	//   iLF->AddEntry( fPowerLawExp_HESS, "Combined Fit B from PRL 101, 261104 (2008)", "l" );
	
	///////////////////////////////////////////////////////////////////////////////////////////////
	// HESS ICRC 2017 spectrum
	
	TF1* fICRC2017_HESS = new TF1( "fICRC2017_HESS",
								   "[0]*TMath::Power( x, 3.-[1] ) * TMath::Power( 1+TMath::Power( x/[3], 1./[4] ), -1.*([2]-[1])*[4] )",
								   0.25, 20. );
	fICRC2017_HESS->SetLineColor( 8 );
	fICRC2017_HESS->SetParameter( 0, 105 / 1.e3 / 1.e3 );
	fICRC2017_HESS->SetParameter( 1, 3.04 );
	fICRC2017_HESS->SetParameter( 2, 3.78 );
	fICRC2017_HESS->SetParameter( 3, 0.94 );
	fICRC2017_HESS->SetParameter( 4, 0.12 );
	
	fICRC2017_HESS->Draw( "same" );
	iLF->AddEntry( fICRC2017_HESS, "HESS 2017", "l" );
	
	
	///////////////////////////////////////////////////////////////////////////////////////////////
	// spectrum from Aharonian & Bogovalov (2003)
	TF1* fAB = new TF1( "fAP", "[0]*TMath::Power( x*1.e3, -1. ) / (1.+TMath::Power( x*1.e3/5., 2.26 ) ) * TMath::Power( x, 3. )", 1.e-2, 1. );
	fAB->SetParameter( 0, 1.36e-3 * 1.e3 * 1.e4 );
	fAB->SetLineColor( 9 );
	//   fAB->Draw( "same" );
	//   iLF->AddEntry( fAB, "Aharonian & Bogovalov (2003)", "l" );
	
	// IFAE spectrum (AM April 15th 2011)
	
	TF1* fIFAE = new TF1( "fIFAE", "6.85e-12*1.e4*1.e3*TMath::Power(x,-3.21) * (1.+0.36*(TMath::Exp( TMath::Gaus( TMath::Log10( x ),-0.5, 0.6 )) -1.)) * TMath::Power( x, 3. )", 1.e-3, 1.e2 );
	// fIFAE->Draw( "same" );
	// iLF->AddEntry( fIFAE, "CTA IFAE (2011)", "l" );
	
	// IFAE spectrum (AM July 25th 2013)
	TF1* fIFAE_2013 = new TF1( "e3elec", "1.e4*x*x*x*2.385e-9*pow(x,-3.43)*(1.+1.95*(exp(TMath::Gaus(log10(x),-0.101, 0.741))-1.))", 0.006, 150. );
	fIFAE_2013->SetLineColor( 6 );
	fIFAE_2013->Draw( "same" );
	
	iLF->AddEntry( fIFAE_2013, "CTA IFAE (2013)", "l" );
	
	// KB (Bernloehr et al 2013 CTA MC paper; Table 3)
	TF1* fHD = new TF1( "HD", powerlaw_plus_lognormal, 1.e-3, 1.e2, 6 );
	fHD->SetParameter( 0, 6.85e-5 );
	fHD->SetParameter( 1, 1. );
	fHD->SetParameter( 2, 3.21 );
	fHD->SetParameter( 3, 0.107 );
	fHD->SetParameter( 4, 3.186e-3 );
	fHD->SetParameter( 5, 0.776 );
	fHD->SetLineStyle( 2 );
	fHD->SetLineWidth( 1 );
	fHD->SetLineColor( 2 );
	fHD->Draw( "same" );
	
	// KB (Bernloher et al 2013 CTA MC paper; Table 3)
	// (cross check)
	/*        TF1* fHD2 = new TF1( "HD2", "x*x*x* ( [0] * TMath::Power( x, [1] ) + [2]/(x*[4]*sqrt(2.*TMath::Pi()) ) *exp( -1. * (log(x/[3])*log(x/[3]))/2./[4]/[4]) )", 1.e-3, 1.e2 );
	//        F_peak(E) = Amp/(x*width*sqrt(2pi))*exp(-0.5*((log(E)-log(Epeak))/width)**2)
		fHD2->SetParameter( 0, 6.85e-5 );
		fHD2->SetParameter( 1, -3.21 );
		fHD2->SetParameter( 2, 3.19e-3 );
		fHD2->SetParameter( 3, 0.107 );
		fHD2->SetParameter( 4, 0.776 );
		fHD2->SetLineColor( 2 );
	        fHD2->SetLineStyle( 2 );
		fHD2->SetLineWidth( 1 );
		fHD2->Draw( "same" ); */
	iLF->AddEntry( fHD, "CTA KB (2010)", "l" );
	
	// Antonino D'Ai: FIT to HESS ICRC 2017 and AMS 2014 data
	// f(x)=Phi0* (x/Eb)**(-Gamma1) * (1+(x/Eb)**(1/alpha))**(-(Gamma2-Gamma1)*alpha)
	TF1* fADA = new TF1( "fADA",
						 "[0] * TMath::Power( x, -[1] ) * TMath::Power( 1+TMath::Power( x/[3], 1./[4] ), -1.*([2]-[1])*[4] )",
						 0.01, 20. );
	fADA->SetParameter( 0, 9.77788e-5 );  // Phi0
	fADA->SetParameter( 1, 0.176473 ); // Gamma1
	fADA->SetParameter( 2, 0.714868 ); //
	fADA->SetParameter( 3, 0.960377 );  // Eb
	fADA->SetParameter( 4, 0.12 ); // alpha
	fADA->SetLineWidth( 1 );
	fADA->SetLineColor( 4 );
	fADA->Draw( "same" );
	iLF->AddEntry( fADA, "AMS (2014) and HESS (ICRC 2017); ADA", "l" );
	
	
	iLF->Draw();
	
	fSpectrum_syst->Draw( "p" );
	
	//////////////////////////////////////////////////////////////////////////////////
	// VERITAS results (2018PhRvD..98f2004A)
	TF1* fVTS_PowerLawBroken = new TF1( "VTS_PowerLawBroken", powerlaw_broken_simple, 316. / 1.e3, 5.623, 4 );
	fVTS_PowerLawBroken->SetLineColor( 800 );
	fVTS_PowerLawBroken->SetParameter( 0, 0.000353289 );
	fVTS_PowerLawBroken->SetParameter( 1, 0.710 );
	fVTS_PowerLawBroken->SetParameter( 2, 3.2 );
	fVTS_PowerLawBroken->SetParameter( 3, 4.1 );
	fVTS_PowerLawBroken->Draw( "same" );
	
	
	//////////////////////////////////////////////////////////////////////////////////
	// calculate values at 1 and 0.1 TeV
	
	cout << "A&B (1TeV) " << fAB->Eval( 1. ) << endl;
	cout << "IFAE (2011) (1TeV) " << fIFAE->Eval( 1. ) << endl;
	cout << "IFAE (2013) (1TeV) " << fIFAE_2013->Eval( 1. ) << endl;
	cout << "HD (1TeV) " << fHD->Eval( 1. ) << endl;
	cout << "DESY (1TeV) " << fPowerLawBroken_DESYAnalysis->Eval( 1. ) << endl;
	
	cout << "M (PBL) (1TeV)  " << fPowerLawBroken->Eval( 1. ) * TMath::Power( 1.e3, 2. )  << endl;
	cout << "M (PBL) (1TeV) " << fPowerLawBroken->Eval( 1. ) << endl;
	cout << "HESS (PBL) (1TeV) " << fPowerLawBroken_HESS->Eval( 1. ) * TMath::Power( 1.e3, 2. )  << endl;
	cout << "HESS (PBL) (1TeV) " << fPowerLawBroken_HESS->Eval( 1. ) << endl;
	cout << "HESS (PBL) (0.1TeV) " << fPowerLawBroken_HESS->Eval( 0.1 ) * TMath::Power( 1.e3, 2. )  << endl;
	cout << "HESS (PBL) (0.1TeV) " << fPowerLawBroken_HESS->Eval( 0.1 ) << endl;
	
	////////////////////////////////////////////////////////////////
	// ratio plot to CTA function
	//
	TCanvas* cR = new TCanvas( "cSR", "electron spectrum (ratio to CTA function)", 110, 10, 900, 600 );
	cR->SetLogx( 1 );
	cR->SetLeftMargin( 0.15 );
	cR->Draw();
	
	TH1D* hR = new TH1D( "hnullR", "", 10, 5. / 1.e3, 30000. / 1.e3 );
	hR->SetMinimum( 0.1 );
	hR->SetMaximum( 2. );
	hR->SetXTitle( "energy [TeV]" );
	hR->SetYTitle( "flux ratio to CTA fit)" );
	hR->GetYaxis()->SetTitleOffset( 1.3 );
	hR->SetStats( 0 );
	hR->Draw();
	
	TLine* iLine1 = new TLine( 5. / 1.e3, 1., 30000. / 1.e3, 1. );
	iLine1->SetLineWidth( 3 );
	iLine1->Draw();
	
	for( unsigned int i = 0; i < inputfile.size(); i++ )
	{
		if( fSpectrum[i] && fIFAE_2013 )
		{
			TGraphAsymmErrors* fT = new TGraphAsymmErrors( fSpectrum[i]->GetN() );
			fT->SetLineColor( fSpectrum[i]->GetLineColor() );
			fT->SetMarkerColor( fSpectrum[i]->GetMarkerColor() );
			fT->SetMarkerStyle( fSpectrum[i]->GetMarkerStyle() );
			double x, y;
			for( int p = 0; p < fSpectrum[i]->GetN(); p++ )
			{
				fSpectrum[i]->GetPoint( p, x, y );
				fT->SetPoint( p, x, y / fIFAE_2013->Eval( x ) );
			}
			fT->Draw( "p" );
		}
	}
	
	
	
	
}

