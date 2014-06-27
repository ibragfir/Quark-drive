#ifndef __PROGTEST__
#include "common.h"
using namespace std;
#endif /* __PROGTEST__ */


pthread_mutex_t mutex;
pthread_mutex_t mutexReact;

sem_t disSem;
sem_t thrSem;

struct Quark{
	unsigned int energy[6]{0};
	CReactNode* root[6]{NULL};//tree that shows how quark was created

	void AddNode(uint8_t formedQuark, unsigned int energy, CReactNode* l, CReactNode* r){
		if(root[formedQuark]){
			root[formedQuark]->m_L = root[formedQuark]->m_R = NULL;
			delete root[formedQuark];
		}
		root[formedQuark] = new CReactNode(formedQuark, energy, l, r);
	}
	void DeleteQuarks(){
		for(int i = 0; i < 6; ++i){
			if(root[i])
				root[i]->m_L = root[i]->m_R = NULL;
			delete root[i];
		}
	}
};
struct Job{
	const TGenerator* gens;
	int gNr;
	int genToSolve;
	void (*engines)(const TRequest*, TSetup *);
	deque<const TRequest*> req_deque;
	bool isFinished;
	TSetup* setup;
	bool *finishedThreads;

	Job(const TGenerator* gens, int gNr, void (*en)(const TRequest*, TSetup *)){
		this->gens = gens;
		this->gNr = gNr;
		engines = en;
		isFinished = false;
		genToSolve = 0;
		setup = new TSetup();
		setup->m_Root = NULL;
		finishedThreads = new bool[gNr];
		for(int i = 0; i < gNr; ++i)
			finishedThreads[i] = false;
	}

	bool AllThreadsFinished(){
		for( int i = 0; i < gNr; ++i)
			if(finishedThreads[i] == false)
				return false;
		return true;
	}
	void SetAllThreadsStatus(){
		for(int i = 0; i < gNr; ++i)
			finishedThreads[i] = false;
	}
};



Quark** CreateTable(const TRequest *request){
	Quark** table;

	table = new Quark*[request->m_FuelNr];
	for(int i = 0; i < request->m_FuelNr; i++){
		table[i] = new Quark[request->m_FuelNr];
		table[0][i].root[request->m_Fuel[i]] = new CReactNode(request->m_Fuel[i], 0);
	}

	return table;
}

void DeleteTable(Quark** table, const TRequest *request){

	for(int i = request->m_FuelNr -1; i >=0; --i){
		for(int j = request->m_FuelNr -1; j >= 0; --j)
			table[i][j].DeleteQuarks();
		delete[] table[i];
	}
   delete[] table;
}

void Reaction(int gen, const TGenerator *gens, const TRequest *req, int tableRow, int tableCol, Quark** table){
	unsigned int prevEnergy;
	unsigned int energy;
	uint8_t ind1;
	uint8_t ind2;
	int ctr = tableRow;
	

	if(tableRow == 1){
		ind1 = req->m_Fuel[tableCol];
		ind2 = req->m_Fuel[ (tableCol+1) % (req->m_FuelNr) ];
		for(uint8_t formedQuark = 0; formedQuark < 6; ++formedQuark){
			energy = gens[gen].m_Energy[ind1][ind2][formedQuark];
			if(energy   &&   energy > table[tableRow][tableCol].energy[formedQuark]){
				table[tableRow][tableCol].energy[formedQuark] = energy;
				CReactNode *left = table[0][tableCol].root[ind1];
				CReactNode *right = table[0][(tableCol+1) % (req->m_FuelNr)] .root[ind2];
				table[tableRow][tableCol].AddNode(formedQuark, energy, left, right);
			}
		}
		return;
	}
	for(int i = 1; i <= tableRow; ++i){
		for(uint8_t j = 0; j < 6; ++j){

			if(i == 1){
				prevEnergy = table[tableRow-i][tableCol].energy[j];
				if(!prevEnergy)
					continue;
				ind1 = j;
				ind2 = req->m_Fuel[ (tableCol+tableRow) % (req->m_FuelNr) ];
				for(uint8_t formedQuark = 0; formedQuark < 6; formedQuark++){
					energy = gens[gen].m_Energy[ind1][ind2][formedQuark];
					if(energy   &&   energy+prevEnergy > table[tableRow][tableCol].energy[formedQuark]){
						table[tableRow][tableCol].energy[formedQuark] = energy+prevEnergy;
						CReactNode *left = table[tableRow-i][tableCol].root[ind1];
						CReactNode *right = table[0][(tableCol+tableRow) % (req->m_FuelNr)].root[ind2];
						table[tableRow][tableCol].AddNode(formedQuark, energy+prevEnergy, left, right);
					}
				}
			}
			else if(i == tableRow){
				prevEnergy = table[tableRow-ctr][(tableCol+ctr) % (req->m_FuelNr)].energy[j];
				if(!prevEnergy)
					continue;
				ind2 = j;
				ind1 = req->m_Fuel[tableCol];
				for(uint8_t formedQuark = 0; formedQuark < 6; formedQuark++){
					energy = gens[gen].m_Energy[ind1][ind2][formedQuark];
					if(energy   &&   energy+prevEnergy > table[tableRow][tableCol].energy[formedQuark]){
						table[tableRow][tableCol].energy[formedQuark] = energy+prevEnergy;
						CReactNode *left = table[0][tableCol].root[ind1];
						CReactNode *right = table[tableRow-ctr][(tableCol+ctr) % (req->m_FuelNr)].root[ind2];
						table[tableRow][tableCol].AddNode(formedQuark, energy+prevEnergy, left, right);
					}
				}
			}
			else{
				prevEnergy = table[tableRow-i][tableCol].energy[j];
				if(!prevEnergy)
					continue;
				ind1 = j;
				for(uint8_t k = 0; k < 6; ++k){
					energy = table[tableRow-ctr][(tableCol+ctr) % (req->m_FuelNr)].energy[k];
					if(!energy)
						continue;
					ind2 = k;
					for(uint8_t formedQuark = 0; formedQuark < 6; formedQuark++){
						int tmp = gens[gen].m_Energy[ind1][ind2][formedQuark];
						if(tmp   &&   tmp+prevEnergy+energy > table[tableRow][tableCol].energy[formedQuark]){
							table[tableRow][tableCol].energy[formedQuark] = tmp+prevEnergy+energy;
							CReactNode *left = table[tableRow-i][tableCol].root[ind1];
							CReactNode *right = table[tableRow-ctr][(tableCol+ctr) % (req->m_FuelNr)].root[ind2];
							table[tableRow][tableCol].AddNode(formedQuark, tmp+prevEnergy+energy, left, right);
						}
					}
				}
			}
		}
		ctr--;
	}
}

CReactNode* GetMaxEnergyRoot(int& startPos, const TRequest *req, Quark** table){
	unsigned int max = 0;
	CReactNode* maxRoot = NULL;
	int fuelNr = req->m_FuelNr;
	uint8_t finProd = req->m_FinalProduct;
	
	for(int i = 0; i < fuelNr; i++){
		if(table[fuelNr-1][i].energy[finProd] > max){
			max = table[fuelNr-1][i].energy[finProd];
			maxRoot = table[fuelNr-1][i].root[finProd];
			startPos = i;
		}
	}
	return maxRoot;
}

CReactNode* CopyTree(CReactNode* root){
	CReactNode* copyNode = NULL;
	
	if(root){
		copyNode = new CReactNode(root->m_Product, root->m_Energy);
		copyNode->m_L = CopyTree(root->m_L);
		copyNode->m_R = CopyTree(root->m_R);
	}
	return copyNode;
}

void CopyMaxRootToSetupRoot(CReactNode* maxRoot, TSetup *setup){
	delete setup->m_Root;
	setup->m_Root = CopyTree(maxRoot);
}

void SolveOneGeneratorOneFuel(int gen, const TGenerator *gens, const TRequest *req, TSetup *setup){
	Quark** table;
	int startPos;

	table = CreateTable(req);
	for(int tableRow = 1; tableRow < req->m_FuelNr; ++tableRow)
		for(int tableCol = 0; tableCol < req->m_FuelNr; ++tableCol)
			Reaction(gen, gens, req, tableRow, tableCol, table);
	CReactNode* root = GetMaxEnergyRoot(startPos, req, table);
	pthread_mutex_lock(&mutexReact);
	if(root  &&  (!setup->m_Root || root->m_Energy > setup->m_Root->m_Energy)){
		CopyMaxRootToSetupRoot(root, setup);
		setup->m_Generator = gen;
		setup->m_Energy = root->m_Energy;
		setup->m_StartPos = startPos;
	}
	pthread_mutex_unlock(&mutexReact);
	DeleteTable(table, req);
}

void optimizeEnergySeq (const TGenerator *generators, int generatorsNr, const TRequest *request, TSetup *setup){
	setup->m_Root = NULL;
	if(!request->m_Fuel || request->m_FuelNr <= 1){
		return;
	}
	for(int g = 0; g < generatorsNr; ++g)
		SolveOneGeneratorOneFuel(g, generators, request, setup);
}                                                              

void *ThrJob1(void* arg){
	Job* job = (Job*) arg;
	const TRequest* req;
	sem_post(&disSem);
	sem_wait(&thrSem);
	TSetup* setup = new TSetup();
	setup->m_Root = NULL;

	while(job->isFinished == false){
		pthread_mutex_lock(&mutex);
			if(job->req_deque.empty()){
				pthread_mutex_unlock(&mutex);
				break;
			}
			req = job->req_deque.back();
			job->req_deque.pop_back();
		pthread_mutex_unlock(&mutex);
		SolveOneGeneratorOneFuel(0, job->gens, req, setup);
		job->engines(req,setup);
		sem_post(&disSem);
		sem_wait(&thrSem);
	}
	delete setup;
	return NULL;
}
void *ThrJob2(void* arg){
	int gen;
	Job* job = (Job*) arg;
	const TRequest* req;
	const TRequest* req2;
	
	while(1){
		sem_wait(&thrSem);	
		pthread_mutex_lock(&mutex);
			if(job->req_deque.empty()  &&  job->isFinished){
				pthread_mutex_unlock(&mutex);
				return NULL;
			}
			req2 = req = job->req_deque.back();
			gen = job->genToSolve;
			job->genToSolve++;
		pthread_mutex_unlock(&mutex);
		if(gen < job->gNr){
			SolveOneGeneratorOneFuel(gen, job->gens, req, job->setup);
			job->finishedThreads[gen] = true;
			if(job->AllThreadsFinished()){
				job->SetAllThreadsStatus();
				job->req_deque.clear();
				
				TSetup* setup2;
				setup2 = job->setup;
				job->setup = new TSetup();
				job->setup->m_Root = NULL;
				sem_post(&disSem);
				job->engines(req2, setup2);
			}
		}
	}
	return NULL;
}

void optimizeEnergy(int threads, const TGenerator* generators, int generatorsNr, 
						const TRequest *(* dispatcher)( void ), void (*engines)(const TRequest*, TSetup *)){
	pthread_t* thrIDs = new pthread_t[threads];
	pthread_attr_t attr;
	Job* job = new Job(generators, generatorsNr, engines);
	const TRequest* req;
	pthread_mutex_init ( &mutex, NULL);
	pthread_mutex_init ( &mutexReact, NULL);
	pthread_attr_init ( &attr );
  	pthread_attr_setdetachstate ( &attr, PTHREAD_CREATE_JOINABLE );
  	sem_init(&disSem, 0, 0);
  	sem_init(&thrSem, 0, 0);

  	if(generatorsNr == 1){
  		for(int i = 0; i < threads; ++i)
  			pthread_create(&thrIDs[i], &attr, ThrJob1, (void*)job);
  		while(1){
  			sem_wait(&disSem);
  			req = dispatcher();
  			if(req == NULL){
  				job->isFinished = true;
  				for(int i = 0; i < threads; ++i)
  					sem_post(&thrSem);
  				break;
  			}
  			pthread_mutex_lock(&mutex);
  			job->req_deque.push_back(req);
  			pthread_mutex_unlock(&mutex);
			sem_post(&thrSem);
  		}
  	}
  	else{
  		for(int i = 0; i < threads; ++i)
  			pthread_create(&thrIDs[i], &attr, ThrJob2, (void*)job);
  		while(job->isFinished == false){
  			req = dispatcher();
  			if(!req){
  				job->isFinished = true;
  				for(int i = 0; i < threads; ++i)
  					sem_post(&thrSem);
  				break;
  			}
  			else{
  				pthread_mutex_lock(&mutex);
  					job->req_deque.push_back(req);
  					job->genToSolve = 0;
  				pthread_mutex_unlock(&mutex);
  			}
  			for(int i = 0; i < generatorsNr; ++i)
  				sem_post(&thrSem);
  			sem_wait(&disSem);
  		}
  	}
  	for(int i = 0; i < threads; ++i)
  		pthread_join(thrIDs[i], NULL);

  	pthread_attr_destroy(&attr);
  	sem_destroy(&disSem);
  	sem_destroy(&thrSem);
  	pthread_mutex_destroy(&mutex);
  	pthread_mutex_destroy(&mutexReact);
  	delete[] job->finishedThreads;
  	delete job->setup;
  	delete job;
  	delete[] thrIDs;
}
