#include "Rcpp.h"
#include <vector>
#include <iostream>

using namespace Rcpp;

// [[Rcpp::export]]

NumericVector fgcv(const NumericMatrix & s) {
    std::vector<unsigned>   informa(0);
    // vector of columns that have more than one 0 and more than one 1.
    // std::cout << "hi\n" ;
    for(unsigned i=0; i<(unsigned) s.ncol(); ++i) {
	std::vector<unsigned> c={0, 0};
	for(unsigned j=0; j<(unsigned) s.nrow(); ++j) {
	    unsigned a=s(j, i);
	    if(!std::isnan(a)) {
		++(c[a>0]);
	    }
	}
	if(c[0] > 1 && c[1] > 1) {
	    informa.push_back(i);
	}
    }
    //   std::cout << "bla\n" ;
    unsigned viola=0; // no. of violations of the 4-gamete condition (4gc)
    if(informa.size()>1) {
	for(unsigned i=0; i<informa.size()-1; ++i) {
	    for(unsigned j=i+1; j<informa.size(); ++j) {
		// check whether colums informa[i] and informa[j] violate 4gc
		unsigned vi = informa[i], vj = informa[j];
		std::vector<std::vector<unsigned>> v={{0,0}, {0,0}};
		for(unsigned k=0; k<(unsigned) s.nrow(); k+=2) {
		    unsigned ki=(s(k, vi)>0),
			kj=(s(k, vj)>0),
			kip=(s(k+1, vi)>0),
			kjp=(s(k+1, vj)>0);
		    // std::cout << k << " " << s.nrow() << " " << vi << " " << vj << " " << s.ncol() << " "  << ki << " " << kj << " "  << kip << " " << kjp << "\n";
		    if(!std::isnan(ki) && !std::isnan(kip) && ki==kip) {
			if(!std::isnan(kj)) {
			    v[ki][kj] = 1;
			}
			if(!std::isnan(kjp)) {
			     v[kip][kjp] = 1;
			}
		    }
		    else if(!std::isnan(kj) && !std::isnan(kjp) && kj==kjp) {
			if(!std::isnan(ki)) {
			    v[ki][kj] = 1;
			}
			if(!std::isnan(kip)) {
			     v[kip][kjp] = 1;
			}
		    }
		}
		// std::cout << vi << '\t' << vj << '\t' <<
		//     v[0][0] << '\t' <<  v[1][0] << '\t' << 
		//     v[0][1] << '\t' <<  v[1][1] << '\n'; 
		if(v[0][0] == 1 && v[0][1] == 1 &&
		   v[1][0] == 1 && v[1][1] == 1) {
		    ++viola;
		}
	    }
	}
    }
    NumericVector res(2);
    res[0] = viola;
    res[1] = informa.size();
    // std::cout << "fin\n";
    return res;
}

// [[Rcpp::export]]

NumericVector jfgcv(const NumericMatrix & s, int samp_1_size) {
    // returns the number of pairs of sites that violate the four gamete condition only if both
    // populations are taken together. (The first 2*samp_1_size lines of the matrix are from population 1, the others from population 2)
    // Also returns the number of informative sites, that is, sites that have at least two 0s and at least two 1s
    std::vector<unsigned>   informa(0);
    // vector of columns that have more than one 0 and more than one 1.
    for(unsigned i=0; i<(unsigned) s.ncol(); ++i) {
	std::vector<unsigned> c={0, 0};
	for(unsigned j=0; j<(unsigned) s.nrow(); ++j) {
	    unsigned a=s(j, i);
	    if(!std::isnan(a)) {
		++(c[a>0]);
	    }
	}
	if(c[0] > 1 && c[1] > 1) {
	    informa.push_back(i);
	}
    }
    unsigned viola=0; // no. of violations of the 4-gamete condition (4gc)
    if(informa.size()>1) {
	for(unsigned i=0; i<informa.size()-1; ++i) {
	    for(unsigned j=i+1; j<informa.size(); ++j) {
		// check whether colums informa[i] and informa[j] informate 4gc
		unsigned vi = informa[i], vj = informa[j];
		std::vector<std::vector<unsigned>> v1={{0,0}, {0,0}},
		    v2={{0,0}, {0,0}}, vb={{0,0}, {0,0}};
		for(unsigned k=0; k<2*(unsigned)samp_1_size; k+=2) {
		    unsigned ki=(s(k, vi)>0),
			kj=(s(k, vj)>0),
			kip=(s(k+1, vi)>0),
			kjp=(s(k+1, vj)>0);
		    if(!std::isnan(ki) && !std::isnan(kip) && ki==kip) {
			if(!std::isnan(kj)) {
			    v1[ki][kj] = 1;
			    vb[ki][kj] = 1;
			}
			if(!std::isnan(kjp)) {
			     v1[kip][kjp] = 1;
			     vb[kip][kjp] = 1;
			}
		    }
		    else if(!std::isnan(kj) && !std::isnan(kjp) && kj==kjp) {
			if(!std::isnan(ki)) {
			    v1[ki][kj] = 1;
			    vb[ki][kj] = 1;
			}
			if(!std::isnan(kip)) {
			     v1[kip][kjp] = 1;
			     vb[kip][kjp] = 1;
			}
		    }
		}
		for(unsigned k=2*samp_1_size; k<(unsigned) s.nrow()-1; k+=2) {
		    unsigned ki=(s(k, vi)>0),
			kj=(s(k, vj)>0),
			kip=(s(k+1, vi)>0),
			kjp=(s(k+1, vj)>0);
		    if(!std::isnan(ki) && !std::isnan(kip) && ki==kip) {
			if(!std::isnan(kj)) {
			    v2[ki][kj] = 1;
			    vb[ki][kj] = 1;
			}
			if(!std::isnan(kjp)) {
			     v2[kip][kjp] = 1;
			     vb[kip][kjp] = 1;
			}
		    }
		    else if(!std::isnan(kj) && !std::isnan(kjp) && kj==kjp) {
			if(!std::isnan(ki)) {
			    v2[ki][kj] = 1;
			    vb[ki][kj] = 1;
			}
			if(!std::isnan(kip)) {
			     v2[kip][kjp] = 1;
			     vb[kip][kjp] = 1;
			}
		    }
		}
		// std::cout << vi << '\t' << vj << '\t' <<
		//     v[0][0] << '\t' <<  v[1][0] << '\t' << 
		//     v[0][1] << '\t' <<  v[1][1] << '\n'; 
		if(vb[0][0] == 1 && vb[0][1] == 1 &&
		   vb[1][0] == 1 && vb[1][1] == 1 &&
		   !(v1[0][0] == 1 && v1[0][1] == 1 &&
		     v1[1][0] == 1 && v1[1][1] == 1) &&
		   !(v2[0][0] == 1 && v2[0][1] == 1 &&
		     v2[1][0] == 1 && v2[1][1] == 1)
		   ) {
		    ++viola;
		}
	    }
	}
    }
    NumericVector res(2);
    res[0] = viola;
    res[1] = informa.size();
    return res;
}
