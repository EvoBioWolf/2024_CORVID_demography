
#include <random>
#include <vector>
#include <iostream>

class NumericMatrix {
public:
    std::vector<std::vector<int>> v;
    NumericMatrix(unsigned i, unsigned j) : v(i, std::vector<int>(j,0)) {};
    unsigned nrow() {return v.size();}
    unsigned ncol() {return v[0].size();}
    int get(unsigned i, unsigned j) {return v[i][j];}
    void set(unsigned i, unsigned j, int x) {v[i][j]=x;}
};

std::vector<double> jfgcv(NumericMatrix & s, int samp_1_size) {
    // returns the number of pairs of sites that violate the four gamete condition only if both
    // populations are taken together. (The first 2*samp_1_size lines of the matrix are from population 1, the others from population 2)
    // Also returns the number of informative sites, that is, sites that have at least two 0s and at least two 1s
    std::vector<unsigned>   informa(0);
    // vector of columns that have more than one 0 and more than one 1.
    for(unsigned i=0; i<(unsigned) s.ncol(); ++i) {
	std::vector<unsigned> c={0, 0};
	for(unsigned j=0; j<(unsigned) s.nrow(); ++j) {
	    unsigned a=s.get(j, i);
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
		    unsigned ki=s.get(k, vi),
			kj=s.get(k, vj),
			kip=s.get(k+1, vi),
			kjp=s.get(k+1, vj);
		    if( ki==kip) {
			    v1[ki][kj] = 1;
			    vb[ki][kj] = 1;
			     v1[kip][kjp] = 1;
			     vb[kip][kjp] = 1;
		    }
		    else if( kj==kjp) {
			    v1[ki][kj] = 1;
			    vb[ki][kj] = 1;
			     v1[kip][kjp] = 1;
			     vb[kip][kjp] = 1;
		    }
		}
		for(unsigned k=2*samp_1_size; k<(unsigned) s.nrow(); k+=2) {
		    unsigned ki=s.get(k, vi),
			kj=s.get(k, vj),
			kip=s.get(k+1, vi),
			kjp=s.get(k+1, vj);
		    if( ki==kip) {
			    v2[ki][kj] = 1;
			    vb[ki][kj] = 1;
			     v2[kip][kjp] = 1;
			     vb[kip][kjp] = 1;
		    }
		    else if( kj==kjp) {
			    v2[ki][kj] = 1;
			    vb[ki][kj] = 1;
			     v2[kip][kjp] = 1;
			     vb[kip][kjp] = 1;
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
    std::vector<double> res(2);
    res[0] = viola;
    res[1] = informa.size();
    return res;
}

int main() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::bernoulli_distribution d(0.25);
    
    NumericMatrix s(200,200);
    for(unsigned i=0; i<200; ++i)
	for(unsigned j=0; j<200; ++j) {
	    s.set(i,j,d(gen));
	}
    std::vector<double> v = jfgcv(s,15);
    for(unsigned k=0; k<v.size(); ++k) std::cout << v[k] << ' ';
    std::cout << '\n';
}
