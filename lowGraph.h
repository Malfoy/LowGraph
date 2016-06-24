#ifndef lowgraph
#define lowgraph

#include <stdio.h>
#include <fstream>
#include <vector>
#include <atomic>
#include <mutex>
#include <unordered_map>
#include <algorithm>
#include "BooPHF.h"


using namespace std;


#define uNumber int32_t
#define kmer uint64_t


typedef boomphf::SingleHashFunctor<uint64_t>  hasher;
typedef boomphf::mphf<uint64_t, hasher  > MPHF;


uint64_t transform_to_size_t(__uint128_t& n){
	return (uint64_t)n;
}


namespace std { template <> struct hash<__uint128_t> {
	typedef __uint128_t argument_type;
	typedef uint64_t result_type; uint64_t operator()(__uint128_t key) const { return transform_to_size_t(key); } };
}


struct unitigIndices{
	kmer overlap;
	uint32_t indice1;
	uint32_t indice2;
	uint32_t indice3;
	uint32_t indice4;
};


kmer rcb(kmer min,uint n){
	kmer res(0);
	kmer offset(1);
	offset<<=(2*n-2);
	for(uint i(0); i<n;++i){
		res+=(3-(min%4))*offset;
		min>>=2;
		offset>>=2;
	}
	return res;
}


kmer str2num(const string& str){
	kmer res(0);
	for(uint i(0);i<str.size();i++){
		res<<=2;
		switch (str[i]){
			case 'A':res+=0;break;
			case 'C':res+=1;break;
			case 'G':res+=2;break;
			default:res+=3;break;
		}
	}
	return res;
}


char revCompChar(char c) {
	switch (c) {
		case 'A': return 'T';
		case 'C': return 'G';
		case 'G': return 'C';
	}
	return 'A';
}


string reverseComplements(const string& s){
	string rc(s.size(),0);
	for (int i((int)s.length() - 1); i >= 0; i--){
		rc[s.size()-1-i]= revCompChar(s[i]);
		// rc[s.size()-1-i]=char2int[(uint)s[i]];
	}
	return rc;
}


class lowGraph{
public:	
	uint coreNumber;
	uint gammaFactor;
	uint k;
	MPHF leftMPHF,rightMPHF;
	vector<unitigIndices> leftIndices,rightIndices;
	vector<string> unitigs;
	
	lowGraph(const string& fileUnitigName, uint kinput,uint gammaFactorinput=1,uint coreNumberinput=1){
		k=kinput;
		coreNumber=coreNumberinput;
		gammaFactor=gammaFactorinput;
		indexUnitigs(fileUnitigName);
	}

vector<pair<string,uNumber>> getEnd(kmer bin){
	vector<pair<string,uNumber>> result;
	kmer rc(rcb(bin,k-1));
	string unitig;
	unitigIndices indices;
	uNumber num;
	bool go(false);
	if(bin<=rc){
		uint64_t hash=rightMPHF.lookup(bin);
		if(hash!=ULLONG_MAX){
			indices=rightIndices[hash];
			if(indices.overlap==bin){
				go=true;
			}
		}
	}else{
		uint64_t hash=leftMPHF.lookup(rc);
		if(hash!=ULLONG_MAX){
			indices=leftIndices[hash];
			if(indices.overlap==rc){
				go=true;
			}
		}
	}
	if(go){
		if(indices.indice1!=0){
			unitig=unitigs[indices.indice1];
			if(str2num(unitig.substr(unitig.size()-k+1,k-1))==bin){
				result.push_back({unitig,indices.indice1});
			}else{
				result.push_back({reverseComplements(unitig),-indices.indice1});
			}
			if(indices.indice2!=0){
				unitig=unitigs[indices.indice2];
				if(str2num(unitig.substr(unitig.size()-k+1,k-1))==bin){
					result.push_back({unitig,indices.indice2});
				}else{
					result.push_back({reverseComplements(unitig),-indices.indice2});
				}
				if(indices.indice3!=0){
					unitig=unitigs[indices.indice3];
					if(str2num(unitig.substr(unitig.size()-k+1,k-1))==bin){
						result.push_back({unitig,indices.indice3});
					}else{
						result.push_back({reverseComplements(unitig),-indices.indice3});
					}
					if(indices.indice4!=0){
						unitig=unitigs[indices.indice4];
						if(str2num(unitig.substr(unitig.size()-k+1,k-1))==bin){
							result.push_back({unitig,indices.indice4});
						}else{
							result.push_back({reverseComplements(unitig),-indices.indice4});
						}
					}
				}
			}
		}
	}
	return result;
}


vector<pair<string,uNumber>> getBegin(kmer bin){
	vector<pair<string,uNumber>> result;
	kmer rc(rcb(bin,k-1));
	string unitig;
	unitigIndices indices;
	bool go(false);
	if(bin<=rc){
		uint64_t hash=leftMPHF.lookup(bin);
		if(hash!=ULLONG_MAX){
			indices=leftIndices[hash];
			if(indices.overlap==bin){
				go=true;
			}
		}
	}else{
		uint64_t hash=rightMPHF.lookup(rc);
		if(hash!=ULLONG_MAX){
			indices=rightIndices[hash];
			if(indices.overlap==rc){
				go=true;
			}
		}
	}
	if(go){
		if(indices.indice1!=0){
			unitig=unitigs[indices.indice1];
			if(str2num(unitig.substr(0,k-1))==bin){
				result.push_back({unitig,indices.indice1});
			}else{
				result.push_back({reverseComplements(unitig),-indices.indice1});
			}
			if(indices.indice2!=0){
				unitig=unitigs[indices.indice2];
				if(str2num(unitig.substr(0,k-1))==bin){
					result.push_back({unitig,indices.indice2});
				}else{
					result.push_back({reverseComplements(unitig),-indices.indice2});
				}
				if(indices.indice3!=0){
					unitig=unitigs[indices.indice3];
					if(str2num(unitig.substr(0,k-1))==bin){
						result.push_back({unitig,indices.indice3});
					}else{
						result.push_back({reverseComplements(unitig),-indices.indice3});
					}
					if(indices.indice4!=0){
						unitig=unitigs[indices.indice4];
						if(str2num(unitig.substr(0,k-1))==bin){
							result.push_back({unitig,indices.indice4});
						}else{
							result.push_back({reverseComplements(unitig),-indices.indice4});
						}
					}
				}
			}
		}
	}
	return result;
}



void indexUnitigs(const string& unitigFileName){
	ifstream unitigFile(unitigFileName);
	unitigs.push_back("");
	string line;
	unitigIndices indices;
	uint leftsize,rightsize,anchorSize;
	vector<kmer>* leftOver=new vector<kmer>;
	vector<kmer>* rightOver=new vector<kmer>;
	vector<kmer>* anchors=new vector<kmer>;
	while(!unitigFile.eof()){
		getline(unitigFile,line);
		getline(unitigFile,line);
		if(line.size()<k){
			break;
		}else{
			unitigs.push_back(line);
			kmer beg(str2num(line.substr(0,k-1))),rcBeg(rcb(beg,k-1));
			if(beg<=rcBeg){
				leftOver->push_back(beg);
			}else{
				rightOver->push_back(rcBeg);
			}
			kmer end(str2num(line.substr(line.size()-k+1,k-1))),rcEnd(rcb(end,k-1));
			if(end<=rcEnd){
				rightOver->push_back(end);
			}else{
				leftOver->push_back(rcEnd);
			}
		}
	}
	sort( leftOver->begin(), leftOver->end() );
	leftOver->erase( unique( leftOver->begin(), leftOver->end() ), leftOver->end() );
	sort( rightOver->begin(), rightOver->end() );
	rightOver->erase( unique( rightOver->begin(), rightOver->end() ), rightOver->end() );
	auto data_iterator = boomphf::range(static_cast<const kmer*>(&((*leftOver)[0])), static_cast<const kmer*>((&(*leftOver)[0])+leftOver->size()));
	leftMPHF= boomphf::mphf<kmer,hasher>(leftOver->size(),data_iterator,coreNumber,gammaFactor,false);
	leftsize=leftOver->size();
	delete leftOver;
	auto data_iterator2 = boomphf::range(static_cast<const kmer*>(&(*rightOver)[0]), static_cast<const kmer*>((&(*rightOver)[0])+rightOver->size()));
	rightMPHF= boomphf::mphf<kmer,hasher>(rightOver->size(),data_iterator2,coreNumber,gammaFactor,false);
	rightsize=rightOver->size();
	delete rightOver;
	anchorSize=anchors->size();
	delete anchors;
	leftIndices.resize(leftsize,{0,0,0,0,0});
	rightIndices.resize(rightsize,{0,0,0,0,0});
	for(uint i(1);i<unitigs.size();++i){
		line=unitigs[i];
		kmer beg(str2num(line.substr(0,k-1))),rcBeg(rcb(beg,k-1));
		if(beg<=rcBeg){
			indices=leftIndices[leftMPHF.lookup(beg)];
			indices.overlap=beg;
			if(indices.indice1==0){
				indices.indice1=i;
			}else if(indices.indice2==0){
				indices.indice2=i;
			}else if(indices.indice3==0){
				indices.indice3=i;
			}else{
				indices.indice4=i;
			}
			leftIndices[leftMPHF.lookup(beg)]=indices;
		}else{
			indices=rightIndices[rightMPHF.lookup(rcBeg)];
			indices.overlap=rcBeg;
			if(indices.indice1==0){
				indices.indice1=i;
			}else if(indices.indice2==0){
				indices.indice2=i;
			}else if(indices.indice3==0){
				indices.indice3=i;
			}else{
				indices.indice4=i;
			}
			rightIndices[rightMPHF.lookup(rcBeg)]=indices;
		}
		kmer end(str2num(line.substr(line.size()-k+1,k-1))),rcEnd(rcb(end,k-1));
		if(end<=rcEnd){
			indices=rightIndices[rightMPHF.lookup(end)];
			indices.overlap=end;
			if(indices.indice1==0){
				indices.indice1=i;
			}else if(indices.indice2==0){
				indices.indice2=i;
			}else if(indices.indice3==0){
				indices.indice3=i;
			}else{
				indices.indice4=i;
			}
			rightIndices[rightMPHF.lookup(end)]=indices;
		}else{
			indices=leftIndices[leftMPHF.lookup(rcEnd)];
			indices.overlap=rcEnd;
			if(indices.indice1==0){
				indices.indice1=i;
			}else if(indices.indice2==0){
				indices.indice2=i;
			}else if(indices.indice3==0){
				indices.indice3=i;
			}else{
				indices.indice4=i;
			}
			leftIndices[leftMPHF.lookup(rcEnd)]=indices;
		}
	}
}


};



#endif

