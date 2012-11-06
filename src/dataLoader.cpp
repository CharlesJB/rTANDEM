/**
 * dataLoader.cpp
 *
 * Author: Charles Joly Beauparlant
 * Creation Date: 18 June 2012
 *
 */

#include "dataLoader.h"

void dataLoader::convertSEXPToMap(SEXP RData, std::map<std::string, std::string>* ptrMap) {
	Rcpp::CharacterVector data(RData);

	for (int i = 0; i < data.size(); i = i + 2) {
		std::string key(data[i]);
		std::string value(data[i+1]);
		(*ptrMap)[key] = value;
	}
}

void dataLoader::convertSEXPToVector(SEXP RData, std::vector<std::string>* ptrVector) {
	Rcpp::CharacterVector data(RData);

	for (int i = 0; i < data.size(); i++) {
		std::string value(data[i]);
		ptrVector->push_back(value);
	}
}

void dataLoader::convertSEXPToDeque(SEXP RData, std::deque<std::string>* ptrDeque) {
	Rcpp::CharacterVector data(RData);

	for (int i = 0; i < data.size(); i = i + 2) {
		std::string value(data[i]);
		ptrDeque->push_back(value);
	}
}
