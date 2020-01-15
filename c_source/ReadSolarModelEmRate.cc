//Code to read the opacity files and extract the required tables from them to be accessed later
//Author: Arshia Ruina

// compile: g++ ReadSolarModel.cc -o ReadSolarModel.out2 -I /usr/include/python2.7/ -I /usr/include/x86_64-linux-gnu/python2.7 -fno-strict-aliasing -Wdate-time -D_FORTIFY_SOURCE=2 -g -fdebug-prefix-map=/build/python2.7-MW0004/python2.7-2.7.15=. -fstack-protector-strong -Wformat -Werror=format-security  -DNDEBUG -g -fwrapv -O2 -Wall -Wstrict-prototypes -L/usr/lib/python2.7/config-x86_64-linux-gnu -L/usr/lib -lpython2.7 -lpthread -ldl  -lutil -lm  -Xlinker -export-dynamic -Wl,-O1 -Wl,-Bsymbolic-functions
// ./ReadSolarModel.out2

// to compile on mac, use g++ -std=c++11 <filename>
// also include sstream and string libraries

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include "ReadSolarModel.h"
//#include "pybind11/pybind11.h"
//#include <python2.7/Python.h>
//#include <pybind11>
#include "pybind11/embed.h"
#include "pybind11/numpy.h"
#include "pybind11/stl.h"         // mandatory for myPyObject.cast<std::vector<T>>()
#include "pybind11/functional.h"
//#include <complex>

namespace py = pybind11;

int main(){	
	SolarModel f;
	//f.ReadAndStoreSolarModel(); 
	//f.MetalMassFraction();
	f.ReadEnergyValues();
	f.ReadOpacityFileName();
	f.AccessSolarModel();
	f.ComputeEmissionRates();
	f.GetEnergyFile();
	f.ComputeFluxFractions();
	
	return 0;
}

double SquaredDistance(double x1, double x2){
	return((x1-x2)*(x1-x2));
}

double ComputeDistance(double x1, double x2){
	return(std::fabs(x1-x2));
}

bool SortDistances(DISVAL& v1, DISVAL& v2){
	return (v1.dist < v2.dist);
}

double TempTokeV(double temp){
	/*
		E = k_B * T
		The conversion factor is nothing but the value
		of the Boltzmann constant in units of [keV/K].
		Value of k_B taken from PDG 2018.
	*/
	return temp*8.617e-8;
}

double DensityTokeV3(double dens){
	/* 
		h_bar*c = 197.326 [MeV*fm]
		which gives a relation between length and electronvolts
		1 [cm**-3] = 7.683e-24 [keV**3]
		Value of h_bar*c taken from PDG 2018.
	*/
	return dens*7.683e-24;
}

int SolarModel::ReadEnergyValues(){

	int i = 0;
        std::string energy_filename = "../resources/axion_gae_flux.dat";
        std::string line;
        std::ifstream energy_file;

	energy_file.open(energy_filename.c_str());

        if(!energy_file.good()){
                std::cout << "Problem with energy file!" << std::endl;
                return 1;
        }

        std::cout << "[INFO] Reading and storing energy values... " << std::endl;

        while(!energy_file.eof()) {
                std::getline(energy_file, line);
                std::istringstream iss_line(line);
                if(line.find("#")==0 || line.empty()) {
                        continue;
                }
                else {
                	energy_vec.push_back(0.0);
			iss_line >> energy_vec[i];
			++i;
		}
	}
	return 0;
}

//double SolarModel::AccessOpacityFile(std::string s, std::string H, std::string He, std::string R, std::string T);

/*--------------------------------------------------------------------------------------------------------
This function reads and stores the contents of the Solar Model in a vector called row
and simultaneously computes and stores the metal mass fractions for each row in another vector called X_Z_metal.
The elements of X_Z_metal will be compared to the opacity filenames to find the closest matching metal composition
for a particular row and later access that opacity table.
----------------------------------------------------------------------------------------------------------*/
int SolarModel::AccessSolarModel() {

	int lineNumber = 0; // counter for each line in the solar model file
	int i = 0; //counter for each row in the solar model data table
	double total_massFrac = 0.0;
	double non_metal_massFrac = 1.0;
	double metal_massFrac_actual = 0.0;
	double metal_massFrac = 0.0;
        //std::vector<ROW> row;
        std::string solarmodel_filename = "../resources/AGSS09_solar_model.dat";
        std::string line;
        std::ifstream solarmodel_file;

	solarmodel_file.open(solarmodel_filename.c_str());

	if(!solarmodel_file.good()){
        	std::cout << "Problem with solar model file!" << std::endl;
		return 1;
	}

	std::cout << "[INFO] Opening output file for storing the computed results..." << std::endl;
	std::stringstream ss1;
	std::ofstream outfile1;
	outfile1.open("../results/electron_densities.txt");
        ss1 << std::setw(15) << std::left << "r [R_sun]";
        ss1 << std::setw(15) << std::left << "Ross mean op";
	ss1 << std::setw(25) << std::left << "e #density [cm**-3]";
	ss1 << std::setw(25) << std::left << "e #density [keV**3]";
	//ss1 << std::setw(20) << std::left << "Compton_emrate []";
	ss1 << std::setw(15) << std::left << "Debye scale";
	ss1 << std::setw(25) << std::left << "Sum(Z(j)^2*n_j) [keV**3]";
	//ss1 << std::setw(20) << std::left << "Brems_emrate";
	ss1 << std::setw(20) << std::left << "S_e";
	ss1 << std::setw(20) << std::left << "S_z";
	//ss1 << std::setw(20) << std::left << "Total_emrate";
	ss1 << std::endl;

        std::cout << "[INFO] Reading solar model file... " << std::endl;
	//++lineNumber;

        while(!solarmodel_file.eof() && lineNumber < 418) { //10% at 217 line
		std::getline(solarmodel_file, line);
                std::istringstream iss_line(line);
                if(line.find("#")==0 || line.empty()) {
                        std::cout << "[INFO] Skipping line " << lineNumber << std::endl;
                        //++lineNumber;
                }
                else {
                        std::cout << "[INFO] Reading line " << lineNumber << " of the Solar Model file = row " << i+1 << " of the Solar Model table" << std::endl;
			row.push_back(ROW());
			//TODO--> Make this efficient!
			iss_line >> row[i].massFrac >> row[i].radius >> row[i].temp >> row[i].density >> row[i].pressure >> row[i].lumiFrac >> row[i].H_massFrac >> row[i].He4_massFrac >> row[i].He3_massFrac >> row[i].C12_massFrac >> row[i].C13_massFrac >> row[i].N14_massFrac >> row[i].N15_massFrac >> row[i].O16_massFrac >> row[i].O17_massFrac >> row[i].O18_massFrac >> row[i].Ne_massFrac >> row[i].Na_massFrac >> row[i].Mg_massFrac >> row[i].Al_massFrac >> row[i].Si_massFrac >> row[i].P_massFrac >> row[i].S_massFrac >> row[i].Cl_massFrac >> row[i].Ar_massFrac >> row[i].K_massFrac >> row[i].Ca_massFrac >> row[i].Sc_massFrac >> row[i].Ti_massFrac >> row[i].V_massFrac >> row[i].Cr_massFrac >> row[i].Mn_massFrac >> row[i].Fe_massFrac >> row[i].Co_massFrac >> row[i].Ni_massFrac;
	
			row[i].temp_keV	= TempTokeV(row[i].temp);
	
			ss1 << std::setw(15) << std::left << row[i].radius;

			/*------------------------------------------------------------------------------------------------------------------
			   The isotopes don't need to be summed up in order to make the comparison with the OP filenames and select an OP file
			   because from the atomic weights mentioned in the OP tables, one can infer that they use the most abundant isotope.
			   That is why the lines below are commented out now and will be removed in future versions.
			------------------------------------------------------------------------------------------------------------------*/
                        //std::cout << "[INFO] Computing and storing isotope mass fractions..." << std::endl;
			//row[i].He_massFrac = row[i].He4_massFrac + row[i].He3_massFrac;
			//row[i].C_massFrac = row[i].C12_massFrac + row[i].C13_massFrac;
			//row[i].N_massFrac = row[i].N14_massFrac + row[i].N15_massFrac;
			//row[i].O_massFrac = row[i].O16_massFrac + row[i].O17_massFrac + row[i].O18_massFrac;
                        
			//std::cout << "[DEBUG] line: " << line << std::endl;
                        //std::cout << "[DEBUG] (iss_line) row contents: " << row[i].massFrac << " " <<  row[i].radius << " " <<  row[i].temp << " " <<  row[i].density << " " <<  row[i].pressure << " " <<  row[i].lumiFrac << " " <<  row[i].H_massFrac << " " <<  row[i].He_massFrac << " " <<  row[i].C12_massFrac << " " <<  row[i].C13_massFrac << " " <<  row[i].N14_massFrac << " " <<  row[i].N15_massFrac << " " <<  row[i].O16_massFrac << " " <<  row[i].O17_massFrac << " " <<  row[i].O18_massFrac << " " <<  row[i].Ne_massFrac << " " <<  row[i].Na_massFrac << " " <<  row[i].Mg_massFrac << " " <<  row[i].Al_massFrac << " " <<  row[i].Si_massFrac << " " <<  row[i].S_massFrac << " " <<  row[i].Ar_massFrac << " " <<  row[i].Ca_massFrac << " " <<  row[i].Cr_massFrac << " " <<  row[i].Mn_massFrac << " " <<  row[i].Fe_massFrac << " " <<  row[i].Ni_massFrac;
        
                        std::cout << "[INFO] Computing and storing total metal fraction and the required metal mass fractions (which will be later used to choose the opacity file)..." << std::endl;
			total_massFrac = row[i].H_massFrac + row[i].He4_massFrac + row[i].He3_massFrac + row[i].C12_massFrac + row[i].C13_massFrac + row[i].N14_massFrac + row[i].N15_massFrac + row[i].O16_massFrac + row[i].O17_massFrac + row[i].O18_massFrac + row[i].Ne_massFrac + row[i].Na_massFrac + row[i].Mg_massFrac + row[i].Al_massFrac + row[i].Si_massFrac + row[i].P_massFrac + row[i].S_massFrac + row[i].Cl_massFrac + row[i].Ar_massFrac + row[i].K_massFrac + row[i].Ca_massFrac + row[i].Sc_massFrac + row[i].Ti_massFrac + row[i].V_massFrac + row[i].Cr_massFrac + row[i].Mn_massFrac + row[i].Fe_massFrac + row[i].Co_massFrac + row[i].Ni_massFrac;
			non_metal_massFrac = row[i].H_massFrac + row[i].He_massFrac;
			metal_massFrac_actual = 1.0 - non_metal_massFrac; 
			metal_massFrac = row[i].C_massFrac + row[i].N_massFrac + row[i].O_massFrac + row[i].Ne_massFrac + row[i].Na_massFrac + row[i].Mg_massFrac + row[i].Al_massFrac + row[i].Si_massFrac + row[i].S_massFrac + row[i].Ar_massFrac + row[i].Ca_massFrac + row[i].Cr_massFrac + row[i].Mn_massFrac + row[i].Fe_massFrac + row[i].Ni_massFrac;
			/* metal_massFrac_actual and metal_massFrac are not equal
			   the former is slightly larger than the latter (difference is of the order 10**-5)
			   reason: while summing for metal_massFrac, we are inlcuding only those elements 
			   which are available in the opacity tables
			   for our purposes, we will consider the metal_massFrac as the total metal mass frac
			   reason: that's what corresponds to the opacity tables 
			*/

			//std::cout << "[DEBUG] total massFrac: " << total_massFrac << std::endl;
			//std::cout << "[DEBUG] H_massFrac: " << row[i].H_massFrac << std::endl;
			//std::cout << "[DEBUG] He_massFrac: " << row[i].He_massFrac << std::endl;
			//std::cout << "[DEBUG] C_massFrac: " << row[i].C_massFrac << std::endl;
			//std::cout << "[DEBUG] N_massFrac: " << row[i].N_massFrac << std::endl;
			//std::cout << "[DEBUG] O_massFrac: " << row[i].O_massFrac << std::endl;	
			//std::cout << "[DEBUG] Non-metal massFrac: " << non_metal_massFrac << std::endl;
			//std::cout << "[DEBUG] Metal massFrac: " << metal_massFrac << std::endl;
			//std::cout << "[DEBUG] Metal massFrac summed: " << metal_massFrac_summed << std::endl;
		
			row[i].total_metal_massFrac = metal_massFrac;
			for(int j=0;j<15;j++)
				row[i].X_Z_metal.push_back(0.0);
			//TODO--> Make this efficient!
			row[i].X_Z_metal[0] = row[i].C_massFrac  / row[i].total_metal_massFrac;	//row[i].X_C  
			row[i].X_Z_metal[1] = row[i].N_massFrac  / row[i].total_metal_massFrac;	//row[i].X_N  
			row[i].X_Z_metal[2] = row[i].O_massFrac  / row[i].total_metal_massFrac;	//row[i].X_O  
			row[i].X_Z_metal[3] = row[i].Ne_massFrac / row[i].total_metal_massFrac;	//row[i].X_Ne 
			row[i].X_Z_metal[4] = row[i].Na_massFrac / row[i].total_metal_massFrac;	//row[i].X_Na 
			row[i].X_Z_metal[5] = row[i].Mg_massFrac / row[i].total_metal_massFrac;	//row[i].X_Mg 
			row[i].X_Z_metal[6] = row[i].Al_massFrac / row[i].total_metal_massFrac;	//row[i].X_Al 
			row[i].X_Z_metal[7] = row[i].Si_massFrac / row[i].total_metal_massFrac;	//row[i].X_Si 
			row[i].X_Z_metal[8] = row[i].S_massFrac  / row[i].total_metal_massFrac;	//row[i].X_S  
			row[i].X_Z_metal[9] = row[i].Ar_massFrac / row[i].total_metal_massFrac;	//row[i].X_Ar 
			row[i].X_Z_metal[10] = row[i].Ca_massFrac / row[i].total_metal_massFrac;	//row[i].X_Ca 
			row[i].X_Z_metal[11] = row[i].Cr_massFrac / row[i].total_metal_massFrac;	//row[i].X_Cr 
			row[i].X_Z_metal[12] = row[i].Mn_massFrac / row[i].total_metal_massFrac;	//row[i].X_Mn 
			row[i].X_Z_metal[13] = row[i].Fe_massFrac / row[i].total_metal_massFrac;	//row[i].X_Fe 
			row[i].X_Z_metal[14] = row[i].Ni_massFrac / row[i].total_metal_massFrac;	//row[i].X_Ni 
	
                        std::cout << "[INFO] Choosing the closest opacity file... " << std::endl;
			std::vector<DISVAL> OpFileDistVec;
			for(unsigned int j=0; j < OpacityFiles.size(); j++) { 
				OpFileDistVec.push_back(DISVAL());
				OpFileDistVec.back().index = j;
				double d = 0.0;
                 		for(int k=0; k < 15; k++) {
                         		d += SquaredDistance(row[i].X_Z_metal[k],std::stod(XZ_VecForOpFile[j][k]));
                 		}	 
				OpFileDistVec.back().dist = std::sqrt(d);
         		}
			//std::sort(OpFileDistVec.begin(),OpFileDistVec.end(),SortDistances);
			sort(OpFileDistVec.begin(),OpFileDistVec.end(),SortDistances);
			int SelectedOpacityFile = OpFileDistVec[0].index;			

			/*
			OpFileDistVec now contains the opacity file names arranged in ascending order
			of their distance to the metal mass fractions of this row.
			*/
			//std::vector<std::string> MatchedOpacityFiles;	
			
			//MatchedOpacityFile.push_back(std::string);
			//MatchedOpacityFile.back() = OpFileDistVec[0].filename;
			//MatchedOpacityFile.push_back(std::string);
			//MatchedOpacityFile.back() = OpFileDistVec[1].filename;
			//MatchedOpacityFile.push_back(std::string);
			//MatchedOpacityFile.back() = OpFileDistVec[2].filename; 

                        std::cout << "[INFO] H mass fraction of this row: " << row[i].H_massFrac << std::endl;
                        std::cout << "[INFO] He mass fraction of this row: " << row[i].He_massFrac << std::endl;
                        std::cout << "[INFO] Choosing the closest H and He mass fractions... " << std::endl;
                        std::vector<DISVAL> HandHe_massFrac_distVec;
                        for(int j=0; j < 126; j++) { //H_massFrac_inOPfiles.size() = 126
                                //double d = SquaredDistance(row[i].H_massFrac,H_massFrac_inOPfiles[j]);
                                double d_H = ComputeDistance(row[i].H_massFrac,H_massFrac_inOPfiles[j]);
                                double d_He = ComputeDistance(row[i].He_massFrac,He_massFrac_inOPfiles[j]);
				HandHe_massFrac_distVec.push_back(DISVAL());
                                //H_massFrac_distVec.back().dist = std::sqrt(d);        
                                HandHe_massFrac_distVec.back().dist = std::sqrt(d_H*d_H + d_He*d_He);
                                HandHe_massFrac_distVec.back().index = j;
                        }
                        std::sort(HandHe_massFrac_distVec.begin(),HandHe_massFrac_distVec.end(),SortDistances);
                        int SelectedHandHemassFrac = -1;
                        if(HandHe_massFrac_distVec.size() > 0){
                                SelectedHandHemassFrac = HandHe_massFrac_distVec[0].index;
                        }
                        std::cout << "[INFO] H mass fraction found to be closest: " << H_massFrac_inOPfiles[SelectedHandHemassFrac] << std::endl;
                        std::cout << "[INFO] He mass fraction found to be closest: " << He_massFrac_inOPfiles[SelectedHandHemassFrac] << std::endl;

 			////--------------------------------------------------------------------------------//
                        //std::cout << "[INFO] H mass fraction of this row: " << row[i].H_massFrac << std::endl;
                        //std::cout << "[INFO] Choosing the closest H mass fraction... " << std::endl;
			//std::vector<DISVAL> H_massFrac_distVec;
			//for(int j=0; j < H_massFrac_inOPfiles.size(); j++) { //H_massFrac_inOPfiles.size() = 126
			//	//double d = SquaredDistance(row[i].H_massFrac,H_massFrac_inOPfiles[j]);
			//	H_massFrac_distVec.push_back(DISVAL());
			//	//H_massFrac_distVec.back().dist = std::sqrt(d);	
			//	H_massFrac_distVec.back().dist = ComputeDistance(row[i].H_massFrac,H_massFrac_inOPfiles[j]);
			//	H_massFrac_distVec.back().value = std::to_string(H_massFrac_inOPfiles[j]);
			//}
			//std::sort(H_massFrac_distVec.begin(),H_massFrac_distVec.end(),SortDistances);
			//std::string SelectedHmassFrac = "";
			//if(H_massFrac_distVec.size() > 0){
			//	SelectedHmassFrac = H_massFrac_distVec[0].value;
			//}
                        //std::cout << "[INFO] H mass fraction found to be closest: " << SelectedHmassFrac << std::endl;

			////--------------------------------------------------------------------------------//
                        //std::cout << "[INFO] He mass fraction of this row: " << row[i].He_massFrac << std::endl;
                        //std::cout << "[INFO] Choosing the closest He mass fraction... " << std::endl;
			//std::vector<DISVAL> He_massFrac_distVec;
			//for(int j=0; j < He_massFrac_inOPfiles.size(); j++) { //He_massFrac_inOPfiles.size() = 126
			//	//double d = SquaredDistance(row[i].He_massFrac,He_massFrac_inOPfiles[j]);
			//	He_massFrac_distVec.push_back(DISVAL());
			//	//He_massFrac_distVec.back().dist = std::sqrt(d);	
			//	He_massFrac_distVec.back().dist = ComputeDistance(row[i].He_massFrac,He_massFrac_inOPfiles[j]);
			//	He_massFrac_distVec.back().value = std::to_string(He_massFrac_inOPfiles[j]);
			//}
			//std::sort(He_massFrac_distVec.begin(),He_massFrac_distVec.end(),SortDistances);
			//std::string SelectedHemassFrac 	= "";
			//if(He_massFrac_distVec.size() > 0) {
			//	SelectedHemassFrac = He_massFrac_distVec[0].value;
			//}
                        //std::cout << "[INFO] He mass fraction found to be closest: " << SelectedHemassFrac << std::endl;

			//--------------------------------------------------------------------------------//
                        std::cout << "[INFO] Density value of this row: " << row[i].density << std::endl;
			double var = std::log10( row[i].density / pow(row[i].temp * 1e-6,3) ); // density in the OP files is stored in this format  
                        std::cout << "[INFO] log R value of this row: " << var << std::endl;
                        std::cout << "[INFO] Choosing the closest log R value... " << std::endl;
			std::vector<DISVAL> logR_distVec;
			for(unsigned int j=0; j < logR.size(); j++) { //logR.size() = 18
				/* x stores the density in a form that makes is comparable to how it is defined
				in the opacity tables: x = log R = log( row[i].density / (row[i].temp * 10e-6) ) */ 
				//double d = SquaredDistance(x,logR[j]);	
				logR_distVec.push_back(DISVAL());	
				//logR_distVec.back().dist = std::sqrt(d);	
				logR_distVec.back().dist = ComputeDistance(var,logR[j]);	
				//logR_distVec.back().value = std::to_string(logR[j]);	
				logR_distVec.back().index = j;//index of the logR value	
			}
			std::sort(logR_distVec.begin(),logR_distVec.end(),SortDistances);
			int SelectedlogR = -1;
			if(logR_distVec.size() > 0) {
				SelectedlogR = logR_distVec[0].index;
                        }
			std::cout << "[INFO] log R found to be closest: " << logR[SelectedlogR] << std::endl;

			//--------------------------------------------------------------------------------//
                        std::cout << "[INFO] Temperature value of this row: " << row[i].temp << std::endl;
                        std::cout << "[INFO] log T value of this row: " << std::log10(row[i].temp) << std::endl;
                        std::cout << "[INFO] Choosing the closest log T value... " << std::endl;
			std::vector<DISVAL> logT_distVec;
			for(int unsigned j=0; j < logT.size(); j++) { //logT.size() = 70
				//double d = SquaredDistance(std::log(row[i].temp),logT[j]);
				logT_distVec.push_back(DISVAL());	
				//logT_distVec.back().dist = std::sqrt(d);	
				logT_distVec.back().dist = ComputeDistance(std::log10(row[i].temp),logT[j]);
				//logT_distVec.back().value = std::to_string(logT[j]);	
				logT_distVec.back().index = j;//index of the logT value	
			}
			std::sort(logT_distVec.begin(),logT_distVec.end(),SortDistances);
			//std::string SelectedlogT = logT_distVec[0].value;
			int SelectedlogT = -1;
			if(logT_distVec.size() > 0) {
				SelectedlogT = logT_distVec[0].index;
                        }
			std::cout << "[INFO] log T found to be closest: " << logT[SelectedlogT] << std::endl;



			// Add up all the distances to find the least distance!
			//std::vector<DISVAL> TotalDistVec;
			//for(j=0; j<3; j++) {
			//TotalDistVec.push_back(DISVAL);
			//TotalDistVec.back().dist = sqrt( std::pow(OpFileDistVec[0].dist,2) + std::pow(H_massFrac_distVec[0].dist,2) + std::pow(He_massFrac_distVec[0].dist,2) + std::pow(logR_distVec[0].dist,2) + std::pow(logT_distVec[0].dist,2) );
			//TotalDistVec.value = OpFileDistVec[j].value;
			//}
			//sort(TotalDistVec.begin(),TotalDistVec.end(),SortDistances);
			
                        std::cout << "[INFO] Calling the function to access the chosen opacity file and obtain the opacity value for the chosen H and He mass fractions, log R and log T values..." << std::endl;
			row[i].opacity_value = pow(10,AccessOpacityFile(SelectedOpacityFile, SelectedHandHemassFrac, SelectedlogR, SelectedlogT));
			std::cout << "Radiative opacity stored for this row: " << row[i].opacity_value << std::endl; 
			ss1 << std::setw(15) << std::left << row[i].opacity_value;
						
			/*-------------------------------------- Number density ------------------------------------*/
		
			std::cout << "[INFO] Number densities n_Z being computed..." << std::endl;
			for(int j=0; j<29; j++)
				row[i].n_Z.push_back(0.0);
			//TODO--> Make this efficient!
			row[i].n_Z[0] = (row[i].H_massFrac / atomic_mass[0]) * (row[i].density / amu); 
			row[i].n_Z[1] = (row[i].He4_massFrac / atomic_mass[1]) * (row[i].density / amu); 
			row[i].n_Z[2] = (row[i].He3_massFrac / atomic_mass[2]) * (row[i].density / amu);
			row[i].n_Z[3] = (row[i].C12_massFrac / atomic_mass[3]) * (row[i].density / amu);
			row[i].n_Z[4] = (row[i].C13_massFrac / atomic_mass[4]) * (row[i].density / amu);
			row[i].n_Z[5] = (row[i].N14_massFrac / atomic_mass[5]) * (row[i].density / amu);
			row[i].n_Z[6] = (row[i].N15_massFrac / atomic_mass[6]) * (row[i].density / amu);
			row[i].n_Z[7] = (row[i].O16_massFrac / atomic_mass[7]) * (row[i].density / amu);
			row[i].n_Z[8] = (row[i].O17_massFrac / atomic_mass[8]) * (row[i].density / amu);
			row[i].n_Z[9] = (row[i].O18_massFrac / atomic_mass[9]) * (row[i].density / amu);
			row[i].n_Z[10] = (row[i].Ne_massFrac / atomic_mass[10]) * (row[i].density / amu);
			row[i].n_Z[11] = (row[i].Na_massFrac / atomic_mass[11]) * (row[i].density / amu);
			row[i].n_Z[12] = (row[i].Mg_massFrac / atomic_mass[12]) * (row[i].density / amu);
			row[i].n_Z[13] = (row[i].Al_massFrac / atomic_mass[13]) * (row[i].density / amu);
			row[i].n_Z[14] = (row[i].Si_massFrac / atomic_mass[14]) * (row[i].density / amu);
			row[i].n_Z[15] = (row[i].P_massFrac / atomic_mass[15]) * (row[i].density / amu);
			row[i].n_Z[16] = (row[i].S_massFrac / atomic_mass[16]) * (row[i].density / amu);
			row[i].n_Z[17] = (row[i].Cl_massFrac / atomic_mass[17]) * (row[i].density / amu);
			row[i].n_Z[18] = (row[i].Ar_massFrac / atomic_mass[18]) * (row[i].density / amu);
			row[i].n_Z[19] = (row[i].K_massFrac / atomic_mass[19]) * (row[i].density / amu);
			row[i].n_Z[20] = (row[i].Ca_massFrac / atomic_mass[20]) * (row[i].density / amu);
			row[i].n_Z[21] = (row[i].Sc_massFrac / atomic_mass[21]) * (row[i].density / amu);
			row[i].n_Z[22] = (row[i].Ti_massFrac / atomic_mass[22]) * (row[i].density / amu);
			row[i].n_Z[23] = (row[i].V_massFrac / atomic_mass[23]) * (row[i].density / amu);
			row[i].n_Z[24] = (row[i].Cr_massFrac / atomic_mass[24]) * (row[i].density / amu);
			row[i].n_Z[25] = (row[i].Mn_massFrac / atomic_mass[25]) * (row[i].density / amu);
			row[i].n_Z[26] = (row[i].Fe_massFrac / atomic_mass[26]) * (row[i].density / amu);
			row[i].n_Z[27] = (row[i].Co_massFrac / atomic_mass[27]) * (row[i].density / amu);
			row[i].n_Z[28] = (row[i].Ni_massFrac / atomic_mass[28]) * (row[i].density / amu);

			std::cout << "[INFO] Number densities converted to units of keV**3" << std::endl;
			for(int j=0; j<29; j++)
				row[i].n_Z_keV3.push_back(DensityTokeV3(row[i].n_Z[j]));
	
			/*-------------------------------------- Electron number density ------------------------------------*/

			std::cout << "[INFO] Computing the electron number density..." << std::endl;
			row[i].n_e = (row[i].density/amu) * (1+row[i].H_massFrac)/2;
			//row[i].n_e_keV = row[i].n_e * 7.645e-24;
			row[i].n_e_keV = DensityTokeV3(row[i].n_e);
			std::cout << "Electron number density stored for this row: " << row[i].n_e << std::endl;
			std::cout << "Electron number density stored for this row [keV]: " << row[i].n_e_keV << std::endl;

			ss1 << std::setw(25) << std::left << row[i].n_e;
			ss1 << std::setw(25) << std::left << row[i].n_e_keV;

			/*--------------------------------- Debye screening scale ------------------------------------*/
	
			std::cout << "[INFO] Computing the Debye screening scale..." << std::endl; 
			//row[i].debye_scale = (4 * M_PI * alpha / row[i].temp) * (row[i].n_e + row[i].n_Z[0] + row[i].n_Z[1]);//making the same approximation as for n_e calculation 
			row[i].debye_scale = std::sqrt( (4 * M_PI * alpha / row[i].temp_keV) * (row[i].n_e + row[i].n_Z[0] + 4*row[i].n_Z[1]) * 7.645e-24 );//making the same approximation as for n_e calculation 
			row[i].y = row[i].debye_scale/(std::sqrt(2*m_e_keV*row[i].temp_keV));
			std::cout << "Debye screening scale for this row: " << row[i].debye_scale << std::endl;
			std::cout << "y for this row: " << row[i].y << std::endl;
			std::cout << "------------------------------------------" << std::endl;
			
			ss1 << std::setw(15) << std::left << row[i].debye_scale;
			ss1 << std::setw(25) << std::left << (row[i].n_Z[0] + 4*row[i].n_Z[1]) * 7.645e-24;		


                       /*----------------------- to compare values with those of M. Giannotti ------------------------*/

                        row[i].Se = (alpha * alpha * (4./3.) * std::sqrt(M_PI) * row[i].n_e_keV * row[i].n_e_keV) / (std::sqrt(row[i].temp_keV) * std::pow(m_e_keV,3.5));
                        //row[i].Se = 1.422e-10 * (pow(470.8,2)) / (std::sqrt(row[i].temp_keV));
                        ss1 << std::setw(20) << std::left << row[i].Se;

                        row[i].Sz = (alpha * alpha * (4./3.) * std::sqrt(2.*M_PI) * row[i].n_e_keV * (row[i].n_Z[0] + 4*row[i].n_Z[1]) * 7.683e-24) / (std::sqrt(row[i].temp_keV) * std::pow(m_e_keV,3.5));
                        //row[i].Sz = (alpha * alpha * (4./3.) * std::sqrt(2*3.142)) / (std::pow(510.998,3.5));
                        ss1 << std::setw(20) << std::left << row[i].Sz;
			
			ss1 << std::endl;
		
			/*--------------------------------------------------------------------------------------------*/

			++i;
		}
		++lineNumber;
	}
	std::cout << std::endl;
	std::cout << "-----------  Finished accessing the solar model ---------"  << std::endl;
	std::cout << row.size() << std::endl;
	outfile1 << ss1.str() << std::endl;
	outfile1.close();
	return 0;
}

void SolarModel::ComputeEmissionRates(){
	
	//TODO:computed values in this function may not need to be stored in the row vector, they can just be stored in output files

	std::stringstream ss_emrate;
	std::ofstream outfile_emrate;
	outfile_emrate.open("../results/emission_rates.txt");

	GaussLagQuad gauss;
    
	for(unsigned int i=0; i<397; i++ ){//
	
		////std::cout << "Computing values for row " << i+1 << " of the Solar Model table" << std::endl;

		for(unsigned int j=200; j<energy_vec.size()-101; j=j+100){ //energy_vec.size()-23574have to start at 200 because flux fraction doesnt work for energies below 0.3 keV 

			/*------------------------------- start of use of energy values ------------------------------*/	

			double energy = energy_vec[j];//keV
			row[i].w = energy/(row[i].temp_keV); //temp in keV;
			std::cout << "energy: " << energy_vec[j] << std::endl;
			std::cout << "w for this row: " << row[i].w << std::endl;

			/*---------------------------------- Absorption coefficient ----------------------------------*/

                        //std::cout << "[INFO] Computing the absorption coefficient..." << std::endl;
			//for(int j=0; j<29; j++)
			//	//row[i].abs_coeff += row[i].n_Z[j]i;
			//	row[i].abs_coeff += row[i].n_Z[j] * 7.645e-24;
			row[i].abs_coeff *= row[i].opacity_value * (1-exp(row[i].w));
			//row[i].abs_coeff *= row[i].opacity_value * (1-exp(row[i].w));
			//std::cout << "Absorption coefficient for this row: " << row[i].abs_coeff << std::endl;
                       
			/* Radiative opacity is defined as the absorption coefficient per unit mass of target
			 * as is expressed in section 2.2 of Redondo's paper after equation 2.6
			 * So the absorption coefficient k(w) used in equation 2.21 will be calculated
			 * as the product of the radiative opacity and the solar density
			 */

			row[i].abs_coeff = row[i].opacity_value * row[i].density;
			

			/*---------------------------------- Compton emission rate -----------------------------------*/

			std::cout << "[INFO] Computing the Compton Emission Rate..." << std::endl;
			//row[i].compton_emrate = (alpha * g_ae * g_ae * energy * energy * row[i].n_e * row[i].n_e)/(3 * m_e * m_e * (exp(row[i].w)-1));
			row[i].compton_emrate = (alpha * g_ae * g_ae * energy * energy * row[i].n_e_keV)/(3. * m_e_keV * m_e_keV * (exp(row[i].w)-1.));
			std::cout << "Compton emission rate for this row: " << row[i].compton_emrate << std::endl;

			//ss1 << std::setw(20) << std::left << row[i].compton_emrate;

			/*------------------------------- Bremsstrahlung emission rate -------------------------------*/
	
			std::cout << "[INFO] Computing the Bremsstrahlung emission rate..." << std::endl; 
			row[i].brems_emrate = (alpha * alpha * g_ae * g_ae * (4./3.) * std::sqrt(M_PI) * row[i].n_e_keV * row[i].n_e_keV * exp(-row[i].w) * gauss.F(row[i].w,std::sqrt(2)*row[i].y)) / (std::sqrt(row[i].temp_keV) * std::pow(m_e_keV,3.5) * energy);
			//ss1 << std::setw(20) << std::left << row[i].brems_emrate;
			std::cout << "Bremsstrahlung emission rate for this row: " << row[i].brems_emrate << std::endl;
			

			/*----------------------------------- Total emission rate ----------------------------------*/
            
			double term1 = (g_ae * g_ae * energy * energy * row[i].abs_coeff) / (2. * e_charge * e_charge * m_e_keV * m_e_keV * (exp(row[i].w)-1.)) ; // includes contribution from ff, fb and bb processes and a part of the Comption contribution
			double term2 = ((exp(row[i].w)-2.) * row[i].compton_emrate) / (2.* (exp(row[i].w)-1.)) ; // completes the Compton contribution
			double term3 = row[i].brems_emrate; // contribution from ee-bremsstahlung
			row[i].total_emrate = (term1 + term2 + term3) / (6.58e-19); // in 1/sec
			std::cout << "Term 1: " << term1 << std::endl;
			std::cout << "Term 2: " << term2 << std::endl;
			std::cout << "Term 3: " << term3 << std::endl;
			std::cout << "Total emission rate for this row: " << row[i].total_emrate << std::endl;
			//total_em[i][j] = row[i].total_emrate;
			

            if (j<energy_vec.size()-201) { //also change if change energy steps ( 101+stepsize= 201) 
			ss_emrate << row[i].total_emrate << "\n";
			}
			else {
			//std::cout << "idiot"  << std::endl;
			ss_emrate << row[i].total_emrate;
		    }

			/*------------------------------ Differential flux fraction ------------------------------- */
			//std::cout << "start" << row[i].diff_fluxfrac << std::endl;	
            /*row[i].diff_fluxfrac_keV = 0.0;
			for(unsigned int l=i; l<row.size(); l++ ){ //row.size()
			
				//std::cout << "START" << std::endl;
                //std::cout << row[i].total_emrate * 1e29 << std::endl;
						/*------------------------------ Differential flux fraction ------------------------------- */	
            	//row[i].w = energy/(row[i].temp_keV);
				/*row[i].k = std::sqrt((energy * energy)- ((4* M_PI* alpha * row[i].n_e_keV) / m_e_keV)); //However, at energies near and below a typical solar plasma frequency, i.e., for energies near or below 0.3 keV,this calculation is not appropriate because the charged particles were treated as staticsources of electric fields, neglecting both recoil effects and collective motions. (those are the first 180-200 entries)
				//if(i != l){
				row[i].diff_fluxfrac_keV +=  (( energy / ((exp(energy/(row[l].temp_keV))-1.0) )) * (std::sqrt((row[l].radius + 0.00049) * (row[l].radius + 0.00049) - row[i].radius * row[i].radius) - std::sqrt((row[l].radius) * (row[l].radius) - row[i].radius * row[i].radius)));  //(1/(4 * M_PI * r_earth * r_earth)) * (2 * energy * energy) / (3 * M_PI)* ((row[i].radius * row[i].radius * row[i].radius) - (row[i-1].radius * row[i-1].radius * row[i-1].radius)) * row[i].total_emrate ;
				//}* row[l].total_emrate)row[i].k
				//else {
					//row[i].diff_fluxfrac_keV = row[i].diff_fluxfrac_keV + (r_sun * r_sun * r_sun * energy * row[i].k * row[l].total_emrate) / ((exp(energy/(row[l].temp_keV))-1.0) * 2 * M_PI * M_PI * M_PI * r_sunearth * r_sunearth); //* row[l].radius / (std::sqrt((row[l].radius) * (row[l].radius) - row[i].radius * row[i].radius))//row[i].diff_fluxfrac = row[i].diff_fluxfrac + ((r_sun * r_sun * r_sun * energy * row[l].k) / (2 * M_PI * M_PI * M_PI * r_sunearth * r_sunearth) * (((row[i].radius * row[i].radius) / (std::sqrt((row[l].radius + 0.00049) * (row[l].radius + 0.00049) - (row[i].radius * row[i].radius)) * (std::sqrt((row[l].radius + 0.00049) * (row[l].radius + 0.00049) - (row[i].radius * row[i].radius))) * (std::sqrt((row[l].radius + 0.00049) * (row[l].radius + 0.00049) - (row[i].radius * row[i].radius))))) -  ((row[i].radius * row[i].radius) / (std::sqrt((row[l].radius +0.00001) * (row[l].radius + 0.00001) - (row[i].radius * row[i].radius)) * (std::sqrt((row[l].radius + 0.00001) * (row[l].radius + 0.00001) - (row[i].radius * row[i].radius))) * std::sqrt((row[l].radius + 0.00001) * (row[l].radius + 0.00001) - (row[i].radius * row[i].radius))))) * row[l].total_emrate / (exp(row[l].w)-1.));  //row[i].diff_fluxfrac = (1/(4 * M_PI * r_earth * r_earth)) * (2 * energy * energy) / (3 * M_PI)* ((row[i].radius * row[i].radius * row[i].radius) ) * row[i].total_emrate   ;
				//}
            	//std::cout << (r_sun * r_sun * r_sun * energy * row[i].k * row[l].total_emrate) / ((exp(energy/(row[l].temp_keV))-1.0) * 2 * M_PI * M_PI * M_PI * r_sunearth * r_sunearth) * 1e29 * (std::sqrt((row[l].radius + 0.00049) * (row[l].radius + 0.00049) - row[i].radius * row[i].radius) - std::sqrt((row[l].radius) * (row[l].radius) - row[i].radius * row[i].radius)) << std::endl;
				//std::cout << row[i].total_emrate * 1e29 << std::endl;
				//std::cout << "ENDE" << std::endl;
			}
			row[i].diff_fluxfrac_keV = row[i].diff_fluxfrac_keV * row[i].total_emrate * (r_sun * r_sun * r_sun * energy  ) / (2 * M_PI * M_PI * M_PI * r_sunearth * r_sunearth);
			row[i].diff_fluxfrac = row[i].diff_fluxfrac_keV * 8.065548 * 8.065548 * 2.41799 * 1e29; // in 1/(cm²s keV) should be 1e29
            ////std::cout << "The differential flux fraction for this row: " << (row[i].diff_fluxfrac ) << std::endl;
			//std::cout << row[i].total_emrate << std::endl;
			std::cout << row[i].diff_fluxfrac  << std::endl;

			////std::cout << "--------------------------------------------------------------------------" << std::endl;*/

		}	
	       
		ss_emrate << std::endl;

	}
	
	outfile_emrate << ss_emrate.str() << 0.0 << std::endl;
	outfile_emrate.close();
}


void SolarModel::GetEnergyFile(){

    std::stringstream ss_energies;
	std::ofstream outfile_energies;
	outfile_energies.open("../results/energies.txt");
    
	for(unsigned int j=200; j<energy_vec.size()-101; j=j+100){ //energy_vec.size()-23574have to start at 200 because flux fraction doesnt work for energies below 0.3 keV
	    double energy = energy_vec[j];//keV
	    if (j<energy_vec.size()-201) { //also change if change energy steps ( 101+stepsize= 201)
			
		    //std::cout << energy  << std::endl;
		    ss_energies << energy << "\n";
		}
		else {
			//std::cout << "idiot"  << std::endl;
			ss_energies << energy;
		}

	}

	outfile_energies <<  ss_energies.str() << std::endl;
	outfile_energies.close();


}

void SolarModel::ComputeFluxFractions(){
	
	
	/*std::stringstream ss_fluxfrac;
	std::ofstream outfile_fluxfrac;
	outfile_fluxfrac.open("../results/flux_fractions.txt");*/

	for(unsigned int i=0; i<row.size(); i++ ){ //101
		
	    // not completely correct jet since the emission rate also depends on the energy and not only on the radius
		//std::cout << "Computing values for row " << i+1 << " of the Solar Model table" << std::endl;

		for(unsigned int j=20000; j<20001 ; j++){ //energy_vec.size() - 3576; j=j+100){ //-23520 star at 200 at least
			double energy = energy_vec[j];//keV
			//std::cout << "LOOP j" << std::endl;
			row[i].diff_fluxfrac_keV = 0.0;
			for(unsigned int l=i; l<row.size(); l++ ){ //row.size()
			
				//std::cout << "START" << std::endl;
                //std::cout << row[i].diff_fluxfrac_keV * 1e29 << std::endl;
						/*------------------------------ Differential flux fraction ------------------------------- */	
            	//row[i].w = energy/(row[i].temp_keV);
				row[i].k = std::sqrt((energy * energy)- ((4* M_PI* alpha * row[i].n_e_keV) / m_e_keV)); //However, at energies near and below a typical solar plasma frequency, i.e., for energies near or below 0.3 keV,this calculation is not appropriate because the charged particles were treated as staticsources of electric fields, neglecting both recoil effects and collective motions. (those are the first 180-200 entries)
				//if(i != l){
				row[i].diff_fluxfrac_keV +=  (r_sun * r_sun * r_sun * energy * row[i].k * row[l].total_emrate) / ((exp(energy/(row[l].temp_keV))-1.0) * 2 * M_PI * M_PI * M_PI * r_sunearth * r_sunearth) * (std::sqrt((row[l].radius + 0.00049) * (row[l].radius + 0.00049) - row[i].radius * row[i].radius) - std::sqrt((row[l].radius) * (row[l].radius) - row[i].radius * row[i].radius));  //(1/(4 * M_PI * r_earth * r_earth)) * (2 * energy * energy) / (3 * M_PI)* ((row[i].radius * row[i].radius * row[i].radius) - (row[i-1].radius * row[i-1].radius * row[i-1].radius)) * row[i].total_emrate ;
				//}
				//else {
					//row[i].diff_fluxfrac_keV = row[i].diff_fluxfrac_keV + (r_sun * r_sun * r_sun * energy * row[i].k * row[l].total_emrate) / ((exp(energy/(row[l].temp_keV))-1.0) * 2 * M_PI * M_PI * M_PI * r_sunearth * r_sunearth); //* row[l].radius / (std::sqrt((row[l].radius) * (row[l].radius) - row[i].radius * row[i].radius))//row[i].diff_fluxfrac = row[i].diff_fluxfrac + ((r_sun * r_sun * r_sun * energy * row[l].k) / (2 * M_PI * M_PI * M_PI * r_sunearth * r_sunearth) * (((row[i].radius * row[i].radius) / (std::sqrt((row[l].radius + 0.00049) * (row[l].radius + 0.00049) - (row[i].radius * row[i].radius)) * (std::sqrt((row[l].radius + 0.00049) * (row[l].radius + 0.00049) - (row[i].radius * row[i].radius))) * (std::sqrt((row[l].radius + 0.00049) * (row[l].radius + 0.00049) - (row[i].radius * row[i].radius))))) -  ((row[i].radius * row[i].radius) / (std::sqrt((row[l].radius +0.00001) * (row[l].radius + 0.00001) - (row[i].radius * row[i].radius)) * (std::sqrt((row[l].radius + 0.00001) * (row[l].radius + 0.00001) - (row[i].radius * row[i].radius))) * std::sqrt((row[l].radius + 0.00001) * (row[l].radius + 0.00001) - (row[i].radius * row[i].radius))))) * row[l].total_emrate / (exp(row[l].w)-1.));  //row[i].diff_fluxfrac = (1/(4 * M_PI * r_earth * r_earth)) * (2 * energy * energy) / (3 * M_PI)* ((row[i].radius * row[i].radius * row[i].radius) ) * row[i].total_emrate   ;
				//}
            	//std::cout << (r_sun * r_sun * r_sun * energy * row[i].k * row[l].total_emrate) / ((exp(energy/(row[l].temp_keV))-1.0) * 2 * M_PI * M_PI * M_PI * r_sunearth * r_sunearth) * 1e29 * (std::sqrt((row[l].radius + 0.00049) * (row[l].radius + 0.00049) - row[i].radius * row[i].radius) - std::sqrt((row[l].radius) * (row[l].radius) - row[i].radius * row[i].radius)) << std::endl;
				//std::cout << row[i].diff_fluxfrac_keV * 1e29 << std::endl;
				//std::cout << "ENDE" << std::endl;
			}
			row[i].diff_fluxfrac = row[i].diff_fluxfrac_keV * 8.065548 * 8.065548 * 2.41799 * 1e32 ; // in 1/(cm²s keV) should be 1e29
			//std::cout << "The differential flux fraction for this row: " << (row[i].diff_fluxfrac ) << std::endl;
            //std::cout << row[i].diff_fluxfrac << std::endl;
			//std::cout << row.size() << std::endl;	//1 
			/*ss_fluxfrac << row[i].diff_fluxfrac << "\t";*/
		}
	       
		/*ss_fluxfrac << std::endl;*/

	}
	
	/*outfile_fluxfrac << ss_fluxfrac.str() << std::endl;
	outfile_fluxfrac.close();*/
}


///*----------------------------------------------------------------------------------------
//This function computes the integral as given in equation 2.16 in l. Redondo's paper.
//It is required in the computation of the Bremsstrahlung emission rate. It is an improper 
//integral in that the bounds of x are 0 and inf. It has been attempted here to solve this 
//https://rosettacode.org/wiki/Numerical_integration/Gauss-Legendre_Quadrature
//----------------------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------------------------------*/

std::vector<double> GaussLagQuad::lagCoef() {

        //double lroots[N];
        //double weight[N];
	std::cout << "N: " << N << std::endl;
        std::vector<double> res_coef;
        //double lcoef[N + 1][N + 1] = {{0}}; // --> uncomment when on linux
        double lcoef[N + 1][N + 1];           // --> uncomment when on mac
        std::memset(lcoef,0,sizeof lcoef); // --> uncomment when on mac

        int n, i;
        lcoef[0][0] = lcoef[1][0] = 1.; lcoef[1][1]  = -1.;//coeffs of the first two polynomials
        for (n = 2; n <= N; n++) { //n-th polynomial
                lcoef[n][0] = 1.; //constants of all Laguerre polynomials are = 1
                for (i = 1; i <= n; i++) //i-th power of x in the n-th polynomial
                        lcoef[n][i] = ( (2*n-1)*lcoef[n-1][i] - lcoef[n-1][i-1] + (1-n)*lcoef[n-2][i] ) / n;
        }

        std::cout << "-------- Coefficients for i-th power (i-th column) of x in the n-th polynomial (n-th row) --------" << std::endl;
        for(int b=0; b<=N; b++){
                std::cout << "\t" << b ;
        }
        std::cout << std::endl;
        for(int a=0; a<=N; a++){
                std::cout << a << "\t";
                for(int b=0; b<=a; b++){
                        std::cout << std::setprecision(4) << std::fixed << lcoef[a][b] << "\t";
                }
                std::cout << std::endl;
        }

        //storing the coefficients of the N-th order polynomial only whose roots will be computed
        for(int i=N; i>=0; i--){
                res_coef.push_back(lcoef[N][i]);
        }

        return res_coef;

}

std::vector<double> GaussLagQuad::lagRootsPy(std::vector<double> poly){

        static py::scoped_interpreter guard{};
        std::cout << "import numpy" << std::endl;
        py::module np = py::module::import("numpy");
        std::cout << "cast poly" << std::endl;
        py::array_t<double> polyNumpy = py::cast(poly);
        std::cout << "access roots" << std::endl;
        py::object roots = np.attr("roots");
        std::cout << "call roots" << std::endl;
        py::object retVal = roots(polyNumpy);
        std::cout << "echo result" << std::endl;
        std::cout << retVal << std::endl;
	std::cout << std::endl << "-------- Finished computing Laguerre polynomials --------" << std::endl << std::endl;
        //return retVal.cast<std::vector<std::complex<double>>>();
        return retVal.cast<std::vector<double>>();

}

double GaussLagQuad::lagEval(int n, double x){
        if(n>0)
                return (2*n-1-x)*GaussLagQuad::lagEval(n-1,x) - (1.-(1./n))*GaussLagQuad::lagEval(n-2,x);
        if(n==1)
                return 1.-x;
        if(n==0)
                return 1.;
	return 0;
}

double GaussLagQuad::lagDeriv(int m, double x){
        if(m>0)
                return GaussLagQuad::lagDeriv(m-1,x) - GaussLagQuad::lagEval(m-1,x);
        if(m==0)
                return 0.; 
	return 0;
}

double GaussLagQuad::quadWeight(double x){
        return 1./(x * pow(GaussLagQuad::lagDeriv(N,x),2) );
}

double GaussLagQuad::inner_integral(double t){
        return (1./2.) * ( ((y*y) / (t*t + y*y)) + log( t*t + y*y ) );
} 

double GaussLagQuad::quadFunc(double x){
	double up_lim = std::sqrt(x+w) + std::sqrt(x);
	double lo_lim = std::sqrt(x+w) - std::sqrt(x);
	return inner_integral(up_lim) - inner_integral(lo_lim);
}

double GaussLagQuad::F(double _w, double _y){

	w = _w;
	y = _y; 	
        double result = 0.0;

	for(int i=0; i<N; i++){
		weights_vec.push_back(quadWeight(roots_vec[i]));
	}

	for(int i=0; i<N; i++){
		result += weights_vec[i] * quadFunc(roots_vec[i]);
	}

	result *= (1./2.);
	////std::cout << "[INFO] F for this row: " << result << std::endl;
	return result;
}

GaussLagQuad::GaussLagQuad() {

        coeff_vec = lagCoef();
	roots_vec = lagRootsPy(coeff_vec);

}

GaussLagQuad::~GaussLagQuad() {}


/*----------------------------------------------------------------------------------------
This function reads all the opacity filenames, extracts the metal mass fractions from them
and stores them in a vector of string vectors called XZ_VecForOpFile where each element is 
a vector of the mass fractions corresponding to the opacity file denoted by its position.
The opacity filenames are stored in a global vector called OpacityFiles. 
----------------------------------------------------------------------------------------*/
void SolarModel::ReadOpacityFileName(){


	for(unsigned int i=0; i < OpacityFiles.size(); i++) {
	//for(int i=0; i < 20; i++)  {
		//std::vector<std::string> X_Z;
		std::cout << "Reading filename: " << OpacityFiles[i] << std::endl;
		//std::string name = OpacityFiles[i];	
		//http://www.cplusplus.com/forum/beginner/68340/d::string name = OpacityFiles[i];	

		//XZ_VecForOpFile.push_back(std::vector<std::string>);
		std::vector<std::string> X_Z_metal;
		XZ_VecForOpFile.push_back(X_Z_metal);
		
		//std::cout << "DEBUG 1" << std::endl;

		std::vector<double> X_Z_metal_numeric;
		XZ_VecForOpFile_numeric.push_back(X_Z_metal_numeric); 

		//std::cout << "DEBUG 2" << std::endl;
		
		std::string name = OpacityFiles[i];
		std::size_t start_pos = name.find("OP17") + 5;
		std::size_t end_pos = name.find("-",start_pos+1); // no need to check for E- for the first end_pos because
							 	  // the carbon mass fraction is given in full decimals
								  // and not in E notation for the files I have 
		
		//std::cout << "DEBUG 3" << std::endl;
	
		std::size_t possible_end_pos = 0;
		//std::cout << "[DEBUG] " << "end_pos: " << end_pos << std::endl;
		//std::cout << "[DEBUG] " << "end: " << name.find("stored") << std::endl;
		//while(end_pos != std::string::npos){
		
		//std::cout << "[DEBUG] " << "name.find(\"stored\"): " << name.find("stored") << std::endl;

		//std::cout << "DEBUG 4" << std::endl;

		while(end_pos < name.find("stored")) {

			//std::cout << "[DEBUG] " << start_pos << std::endl;
			//std::cout << "[DEBUG] " << end_pos << std::endl;

			//std::cout << "DEBUG 5" << std::endl;
			
			std::string cooz = name.substr(start_pos,end_pos-start_pos);
			//std::cout << "[DEBUG] " << cooz << std::endl;

			X_Z_metal.push_back(cooz);
			X_Z_metal_numeric.push_back(std::stod(cooz));

		 	start_pos = end_pos + 1;	
			possible_end_pos = name.find("-",start_pos+1);
			//std::cout << "[DEBUG] " << "possible_end_pos:  " << possible_end_pos << std::endl;

			if(possible_end_pos!=std::string::npos && name.substr(possible_end_pos-1,1) == "E") { 
				//std::cout << "[DEBUG] " << "Found an " << name.substr(possible_end_pos-1,1) << std::endl;
				end_pos = name.find("-",possible_end_pos+1);
			}
			else
				end_pos = possible_end_pos;
			//std::cout << "[DEBUG] " << "start_pos: " << start_pos << std::endl;
			//std::cout << "[DEBUG] " << "end_pos: " << end_pos << std::endl;

		}

		//std::cout << "[DEBUG] " << "Hi!" << start_pos << "\t" << end_pos << std::endl;
		//std::cout << "DEBUG 6" << name.length() << std::endl;
		//std::cout << "DEBUG 7" << std::endl;
		std::string lastcooz = name.substr(start_pos,name.length()-7-start_pos);
		//std::cout << "[DEBUG] " << lastcooz << std::endl;	
		X_Z_metal.push_back(lastcooz);
		X_Z_metal_numeric.push_back(std::stod(lastcooz));
		//std::cout << "[DEBUG] " << "Length of string vector: " << X_Z_metal.size() << std::endl;	
		XZ_VecForOpFile.back() = X_Z_metal;
		XZ_VecForOpFile_numeric.back() = X_Z_metal_numeric;
		//std::cout << "[DEBUG] " << "XZ_VecForOpFile stored" << std::endl;
		//std::cout << "[DEBUG] " << X_Z_metal[0] << std::endl;
		//std::cout << "[DEBUG] " << XZ_VecForOpFile[0][0] << std::endl;
		//std::cout << "DEBUG 8" << std::endl;
	}	

	std::cout << std::endl;
	std::cout << "OpacityFiles.size() = " << OpacityFiles.size() << std::endl; 
	std::cout << "XZ_VecForOpFile.size() = " << XZ_VecForOpFile.size() << std::endl; 

	for(unsigned int i=0; i < OpacityFiles.size(); i++) {
	//for(int i=0; i < 20; i++) {
		std::cout << "File: " << OpacityFiles[i] << std::endl;
		for(int j=0; j < 15; j++) {
			std::cout << "XZ_VecForOpFile: " << XZ_VecForOpFile[i][j] << std::endl;
			std::cout << "XZ_VecForOpFile_numeric: " << XZ_VecForOpFile_numeric[i][j] << std::endl;
		}
		std::cout << "End of file" << std::endl;
	} 
}

/*
compare XZ vector for ith row with XZ vector for jth file
find distance between row[i].X_Z_metal[k] and XZ_VecForOpFile[j][k]
*/

/*----------------------------------------------------------------------------------------
This function compares the metal mass fraction vectors for a row from the Solar Model 
to those extracted from the opacity filenames. The latter have to be first coverted
to numeric type with the correct precision before making the comparison.

For every row, a distance is computed and stored corresponding to every opacity filename.
The distances are then sorted such that we can choose the filename with the least distance
for that row. 
-----------------------------------------------------------------------------------------*/	
/*int SolarModel::CompareMetalMassFractions() { 
	
	for(int i=0; i < OpacityFiles.size(); i++) {
		for(int j=0; j < 15; j++) {
			SquareDistance(row[i],XZ_VecForOpFile[i][j]
		}	
	}

}*/


/*---------------------------------------------------------------------------------------
This function takes as arguments the indices of the selected opacity file, H and He mass 
fractions (thereby the table index inside an opacity file), the R value and the T value 
to access the opacity file, the table inside it and the position inside the table so that 
it can return the value of the Rosseland mean cross-section for each row in the solar model.
---------------------------------------------------------------------------------------*/

//double AccessOpacityFile(std::string s, std::string H, std::string He, std::string R, std::string T) {
double SolarModel::AccessOpacityFile(int s, int HandHe, int R, int T) {

        int lineNumber = 0;
	int startTable_lineNumber = 0;                                
        //int posX = 0;
        //int posY = 0;
	std::string tableIndex;
	std::string selected_tableIndex;
	//std::vector<std::string> tableIndex;
			

	std::ifstream opacity_file;
	opacity_file.open(OpacityFiles[s].c_str());

	std::cout << "[INFO] Accessing opacity file: " << OpacityFiles[s].c_str() << std::endl;
	std::cout << "[INFO] Looking for H mass fraction: " << H_massFrac_inOPfiles[HandHe] << std::endl;
	std::cout << "[INFO] Looking for He mass fraction: " << He_massFrac_inOPfiles[HandHe] << std::endl;
	std::cout << "[INFO] Looking for log R value: " << logR[R] << std::endl;
	std::cout << "[INFO] Looking for log T value: " << logT[T] << std::endl;

	if(!opacity_file.good()){
		std::cout << "Problem with file!" << std::endl;
		return 1;
	}

	while(!opacity_file.eof()){
		std::string op_line = "";
		std::getline(opacity_file, op_line);
                std::istringstream iss_op_line(op_line);

		++lineNumber;	
		
		//std::cout << "[INFO] Finding table index by comparing the H and He mass fractions..." << std::endl;
		//if(lineNumber > 63 && lineNumber < 189){

		//	posX = op_line.find("X");
		//	posY = op_line.find("Y");
		//	std::string _H_massFrac = op_line.substr(posX+2,6);
		//	std::string _He_massFrac = op_line.substr(posY+2,6);
		//	//std::cout << "[DEBUG] _H_massFrac: " << _H_massFrac << std::endl;
		//	//std::cout << "[DEBUG] _He_massFrac: " << _He_massFrac << std::endl;
		//	
		//	//if(std::stod(_H_massFrac) == row[i].H_massFrac && std::stod(_He_massFrac) == row[i].He_massFrac){
		//	//if(std::stod(_H_massFrac) == H.c_str() && std::stod(_He_massFrac) == He.c_str()){
		//	//if(_H_massFrac == H.c_str() && _He_massFrac == He.c_str()){
		//	if(std::stod(_H_massFrac) == H_massFrac_inOPfiles[HandHe] && std::stod(_He_massFrac) == He_massFrac_inOPfiles[HandHe]){
        	//		selected_tableIndex = op_line.substr(8,3);
		//		//tableIndex.push_back(line.substr(8,3));
        	//		//select[i].tableIndex = tableIndex;
        	//		//std::cout << "[DEBUG] Selected mass fraction of H: " << select[i].H_massFrac << std::endl;
        	//		//std::cout << "[DEBUG] Selected mass fraction of He: " << select[i].He_massFrac << std::endl;
       		//		std::cout << "[INFO] Selected table: " << selected_tableIndex << std::endl;
       		//		//std::cout << "[DEBUG] Aur kuch? " << std::endl;
       		//	}		
		//}
		//--------------------------------------------------------------//
		// above lines replaced by -->
		//selected_tableIndex = std::to_string(HandHe+1); //--> not required, just use HandHe+1 in the if-statement for selection
		//--------------------------------------------------------------//
	
		//std::cout << "Selected table number: " << selected_tableIndex << std::endl; 
		//std::cout << "Locating the selected table in the opacity file..." << selected_tableIndex << std::endl; 

		if(lineNumber > 240) {
		//if(lineNumber > 3850 && lineNumber < 3950) { //testing for table 48

			//std::string tableIndex;
			std::string log_rad_opa = "";
                        std::size_t posTableIndex = op_line.find("TABLE");

                        if (posTableIndex!=std::string::npos){ // 1. if at the first line of a table
                        	//std::cout << "DEBUG: at first line of a table " << line << std::endl;
                                tableIndex = op_line.substr(posTableIndex+7,3);
                                //if(tableIndex == selected_tableIndex){ // 1.a) if at the start of a table we want to access
                                if(std::stoi(tableIndex) == HandHe+1){ // 1.a) if at the start of a table we want to access
                                	//std::cout << "DEBUG: at first line of a table we want to access" << std::endl;
                                        //std::cout << "DEBUG: " << line << std::endl;
                                        //s << line << std::endl;
                                        //select[i].lineNumber = lineNumber; // storing the location of start of that table
					//std::cout << "[DEBUG 1] tableIndex: " << tableIndex << std::endl;
					//std::cout << "[DEBUG 2] selected_tableIndex: " << selected_tableIndex << std::endl;
					//std::cout << "[DEBUG 3] lineNumber " << lineNumber << std::endl;
                                        startTable_lineNumber = lineNumber; // storing the location of start of that table
                                        //std::cout << "[DEBUG 4] Selected table starts on line: " << startTable_lineNumber << std::endl;
                                	//std::cout << "[DEBUG 5] startTable_lineNumber+5 " << startTable_lineNumber+5 << std::endl;
                                	//std::cout << "[DEBUG 6] startTable_lineNumber+76 " << startTable_lineNumber+76 << std::endl;
                                                 //foundTable = true;
                                }
                                else { // 1.b) if at the start of a table we don't want to access
                                       //std::cout << "DEBUG: at the start of a table we don't want to access" << std::endl;
                                	continue;
                                }
			}
                        else { // 2. if inside a table
                                //std::cout << "[DEBUG] inside a table" << std::endl;
                                //std::cout << "[DEBUG] line number " << lineNumber << std::endl;
                                //std::cout << "[DEBUG] startTable_lineNumber " << startTable_lineNumber << std::endl;
                                //std::cout << "[DEBUG] line number " << lineNumber << std::endl;
				//std::cout << "[DEBUG 7] startTable_lineNumber+5 " << startTable_lineNumber+5 << std::endl;
                                //std::cout << "[DEBUG 8] startTable_lineNumber+76 " << startTable_lineNumber+76 << std::endl;
                        	//if(lineNumber > startTable_lineNumber+5 ) {std::cout << "[DEBUG 9] line number: " << lineNumber << " and lineNumber > startTable_lineNumber+5" << std::endl;}
                        	//if(lineNumber < startTable_lineNumber+76 ) {std::cout << "[DEBUG 10] lineNumber: " << lineNumber << " and lineNumber < startTable_lineNumber+76" << std::endl;}
                        	if(lineNumber > startTable_lineNumber+5 && lineNumber < startTable_lineNumber+76 ) { // 2.a) if inside a table we want to access
                                	//std::cout << "[INFO] Inside a table we want to access" << std::endl;
                                        //s << line << std::endl;
                                        //std::cout << "DEBUG: " << line << std::endl;
					//std::cout << " [INFO] lineNumber: " << lineNumber << std::endl;
                        		std::vector<std::string> row_in_OPtable;
					std::string word = "";
					while(iss_op_line >> word){ //iss_op_line: way to access the line number as a string
						row_in_OPtable.push_back(word);
					}
					//std::cout << " [DEBUG 11]" << std::endl;
					//iss_op_line >> row_in_OPtable[0] >> row_in_OPtable[1] >> row_in_OPtable[2] >> row_in_OPtable[3] >> row_in_OPtable[4] >> row_in_OPtable[5] >> row_in_OPtable[6] >> row_in_OPtable[7] >> row_in_OPtable[8] >> row_in_OPtable[9] >> row_in_OPtable[10] >> row_in_OPtable[11] >> row_in_OPtable[12] >> row_in_OPtable[13] >> row_in_OPtable[14] >> row_in_OPtable[15] >> row_in_OPtable[16] >> row_in_OPtable[17] >> row_in_OPtable[18] >> row_in_OPtable[19];   
					//std::cout << " [DEBUG 12]" << std::endl;
					std::string OPtable_temp = row_in_OPtable[0];
					//std::cout << " [DEBUG 13]" << std::endl;
					//if(row_in_OPtable[0] == logT[std::stoi(T)])
					//std::cout << "[DEBUG 14] OPtable_temp: " << OPtable_temp << std::endl;
					if(std::stod(OPtable_temp) == logT[T]) {
						log_rad_opa = row_in_OPtable[R+1];
						std::cout << "[INFO] Found the best log_rad_opa for this row: " << log_rad_opa << std::endl;
						return std::stod(log_rad_opa);
					}
				}
                                else { continue;// 2.b) if inside a table we don't want to access
                                       //std::cout << "DEBUG: inside a table we don't want to access" << std::endl;
                                }
			}			
		}
	}
	return 0;
}
/*
	------  main comes here ------

	int main(){
	    // create a MyClass object
	    MyClass obj;
	    
	    for (auto el : obj.myNumbers){
	        std::cout << "Element: " << el << std::endl;
	    }
	    // destructor will be called automatically, since it's a normal object and not
	    // a pointer created by `new` or something similar
		    return 0;
	}

*/
