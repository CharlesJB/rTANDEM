/**
 * dataLoader.h
 *
 * Author: Charles Joly Beauparlant
 * Creation Date: 18 June 2012
 *
 */

#ifndef DATALOADER_H
#define DATALOADER_H

#include <Rcpp.h>
#include <map>
#include <string>

class dataLoader {
public:
	static void convertSEXPToMap(SEXP RData, std::map<std::string, std::string>* ptrMap);
	static void convertSEXPToVector(SEXP RData, std::vector<std::string>* ptrVector);
	static void convertSEXPToDeque(SEXP RData, std::deque<std::string>* ptrDeque);
	
private:
	dataLoader() { }
};

#endif /* DATALOADER_H */
