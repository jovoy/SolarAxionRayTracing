#ifndef READFILE_H
#define READFILE_H

struct SELECT {
	std::string H_massFrac = "";
	std::string He_massFrac = "";
	std::string tableIndex = "";
	int lineNumber = 0;
};

struct FILENAME {
	double  
}

class ReadOpacityFile{
private:
    	std::vector<SELECT> select;
    	std::vector<std::string> H_massFrac_list = {"0.9980","0.9999"};
    	std::vector<std::string> He_massFrac_list = {"0.0000","0.0000"};
	int input = H_massFrac_list.size();//or size of He_massFrac_list
    	std::string in_filename = "../resources/opacity_tables/OP17.1.0-0.0-0.0-0.0-0.0-0.0-0.0-0.0-0.0-0.0-0.0-0.0-0.0-0.0-0.0.stored";
    	std::string line;
    	std::ifstream in_file;
public:
    	int GetTableIndex();
    	int ReadAndStoreTable();
};

#endif // READFILE_H
