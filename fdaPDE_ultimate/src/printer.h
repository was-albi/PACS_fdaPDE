#ifndef _PRINTER_HPP_
#define _PRINTER_HPP_

#include<fstream>
#include<string>
#include"fdaPDE.h"
#include <iomanip>

class printer{


// std::string directory = "/Users/giuliopn/PACSworkspace3/PACSworkspace/GAM_tests/debugging_output/CPP/";

// std::string directory = "/home/alb/Scrivania/PACS/Git_Folder/debugging_output/CPP/";

public:

static void SaveMatrixXr(std::string& name_txt, const MatrixXr & mat){
		const UInt numRows = mat.rows();
		const UInt numCols = mat.cols();
		constexpr UInt CONSTPREC = 16;


		std::string directory = "/home/alb/Scrivania/PACS/Git_Folder/debugging_output/CPP/";
		name_txt = directory + name_txt;

		std::ofstream FileDataMatrix(name_txt);
		if(FileDataMatrix.is_open()){
			for(int i=0; i<numRows; i++){
				for(int j=0; j<numCols; j++){
					FileDataMatrix<< std::scientific << mat(i,j);
					if( j< numCols -1 )FileDataMatrix<< ',';
				}
			FileDataMatrix<<'\n';
			}
		}
		FileDataMatrix.close();
}

template<typename T>
static void SaveDimension(std::string& name_txt, const std::vector<T> & vec){
		const UInt numelem = vec.size();
		constexpr UInt CONSTPREC = 16;


		std::string directory = "/home/alb/Scrivania/PACS/Git_Folder/debugging_output/CPP/";
		name_txt = directory + name_txt;

		std::ofstream FileDataMatrix(name_txt);
		if(FileDataMatrix.is_open()){
				FileDataMatrix<< vec.size();
		}
		FileDataMatrix.close();
}


static void saveVectorXr(std::string& name_txt, const VectorXr & vect){
   	const UInt size = vect.size();
		constexpr UInt CONSTPREC = 16;

	    std::string directory = "/home/alb/Scrivania/PACS/Git_Folder/debugging_output/CPP/";
		name_txt = directory + name_txt;

	std::ofstream File(name_txt);
	if(File.is_open()){
		for(int k=0; k<size; k++){
				File<< std::scientific <<vect[k]<<'\n';
		}
	}
	File.close();
}

};

#endif
