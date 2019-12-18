//Code to read the opacity files and extract the required tables from them to be accessed later
//Author: Arshia Ruina

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include "ReadOpacityFile.h"

int main(){
	ReadOpacityFile f;
	f.ReadFileName();
	f.GetTableIndex();
	f.ReadAndStoreTable();
	return 0;
}

//int readSolarModel::{
//}

void ReadOpacityFile::ReadFileName() {
	
}

 
int ReadOpacityFile::GetTableIndex() {

	int lineNumber = 0;
	int posX = 0;
	int posY = 0;

	std::cout << "size of list: " << input << std::endl;

	for(int i = 0; i < input; i++) {
		select.push_back(SELECT());
		select.back().H_massFrac = H_massFrac_list[i];
		select.back().He_massFrac = He_massFrac_list[i];
	}
	
	in_file.open(in_filename.c_str());

	if(!in_file.good()){
        	std::cout << "Problem with file!" << std::endl;
		return 1;
	}

	while(!in_file.eof()){
		std::getline(in_file, line);
                ++lineNumber;

		if(lineNumber > 63 && lineNumber < 189){		
			for(int i = 0; i < input; i++){
				//std::istringstream iss_line(line);
				std::string tableIndex;
				posX = line.find("X");
				posY = line.find("Y");
				std::string H_massFrac = line.substr(posX+2,6);
				std::string He_massFrac = line.substr(posY+2,6);
				if(H_massFrac == select[i].H_massFrac && He_massFrac == select[i].He_massFrac){ 
					tableIndex = line.substr(8,3);
                                        select[i].tableIndex = tableIndex;
					std::cout << "Selected mass fraction of H: " << select[i].H_massFrac << std::endl;
					std::cout << "Selected mass fraction of He: " << select[i].He_massFrac << std::endl;
					std::cout << "Selected table: " << select[i].tableIndex << std::endl;
				}	
			}
		}
	}
        in_file.close();
}		

int ReadOpacityFile::ReadAndStoreTable() {

        for(int i = 0; i < input; i++){

        	int lineNumber = 0;
        	int posTableIndex = 0;

        	in_file.open(in_filename.c_str());
	
	        if(!in_file.good()){
	                std::cout << "Problem with file!" << std::endl;
	                return 1;
	        }

                std::stringstream s;
                std::ofstream out_file;
                std::string out_filename = "../resources/extracted_tables/" + select[i].H_massFrac + "-" + select[i].He_massFrac + "-" + in_filename.substr(32,60) + ".stored";
                out_file.open(out_filename);

                while(!in_file.eof()){
                        std::getline(in_file, line);
                        ++lineNumber;
                        if(lineNumber > 240 ){

                                std::string tableIndex;
                                std::size_t posTableIndex = line.find("TABLE");

                                if (posTableIndex!=std::string::npos){ // 1. if at the first line of a table
					//std::cout << "DEBUG: at first line of a table " << line << std::endl;
				        tableIndex = line.substr(posTableIndex+7,3);
                                        if(tableIndex == select[i].tableIndex){	// 1.a) if at the start of a table we want to access
                                                //std::cout << "DEBUG: at first line of a table we want to access" << std::endl;
						//std::cout << "DEBUG: " << line << std::endl;
						s << line << std::endl;
                                                select[i].lineNumber = lineNumber; // storing the location of start of that table
                                                std::cout << "Table number " << tableIndex << " starts on line " << lineNumber << std::endl;
                                                //foundTable = true;
                                        }
					else { // 1.b) if at the start of a table we don't want to access
						//std::cout << "DEBUG: at the start of a table we don't want to access" << std::endl;
						continue;
                                	}
				}
                                else { // 2. if inside a table
					//std::cout << "DEBUG: inside a table" << std::endl;
					if(lineNumber < select[i].lineNumber+76) { // 2.a) if inside a table we want to access
                                		//std::cout << "DEBUG: inside a table we want to access" << std::endl;
						s << line << std::endl;
						//std::cout << "DEBUG: " << line << std::endl;
                                	}
					else { // 2.b) if inside a table we don't want to access
						//std::cout << "DEBUG: inside a table we don't want to access" << std::endl;
					}
				}
                        }
                }
                std::cout << "Finished parsing table" << std::endl;
                std::cout << "Storing the parsed table " << select[i].tableIndex << " in ouput file " << out_filename << std::endl;
                out_file << s.str() << std::endl;
                out_file.close();
                std::cout << "Output file closed" << std::endl;
                in_file.close();
                std::cout << "Input file closed" << std::endl;
        }

        std::cout << "Done!" << std::endl;
}
