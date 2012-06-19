/**
 * dataLoader.cpp
 *
 * Author: Charles Joly Beauparlant
 * Creation Date: 18 June 2012
 *
 */

#include "dataLoader.h"

void dataLoader::convertSEXP(SEXP RData, std::map<std::string, std::string>* ptrMap) {
	Rcpp::CharacterVector data(RData);

	for (size_t i = 0; i < data.size(); i = i + 2) {
		std::string key(data[i]);
		std::string value(data[i+1]);
		(*ptrMap)[key] = value;
	}
}
