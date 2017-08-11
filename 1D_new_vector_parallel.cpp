#include <iostream>
#include <cstdlib>
#include <vector>
#include <ctime>
#include <array>
#include <cmath>
#include <trng/yarn2.hpp>
#include <trng/uniform01_dist.hpp>
#include <omp.h>
#include <sys/time.h>
#include <string>
using namespace std;

double TL = 1.0;
double TR = 2.0;

enum Rate_Func {
	rf1,
	rf2
};

struct thread_info {
	int current_rank = -1;
	double running_time = 0.0;
	int count = 0;
};

Rate_Func rate_func = rf1;

class interaction {
	public: 
		double time;
		int location;
		
		bool operator< (const interaction &right) {
			return (this -> time < right.time);
		}
		
		bool operator<= (const interaction &right) {
			return (this -> time <= right.time);
		}
		
		bool operator> (const interaction &right) {
			return (this -> time > right.time);
		}
		
		bool operator>= (const interaction &right) {
			return (this -> time >= right.time);
		}
		
		bool operator== (const interaction &right) {
			return (this -> time == right.time && this -> location == right.location);
		}
		
		bool operator< (const double &right) {
			return (this -> time < right);
		}
		
		bool operator<= (const double &right) {
			return (this -> time <= right);
		}
		
		bool operator> (const double &right) {
			return (this -> time > right);
		}
		
		bool operator>= (const double &right) {
			return (this -> time >= right);
		}
		
		void peek_info() {
			//cout<<"Time: "<<this -> time << " Location: "<< this -> location<<endl;
			printf("Time: %.5lf, Location: %d\n", this->time, this->location);
		}
};

template <class T>
class Vector : public vector<T> {
	
	private:
		bool remove_elt(unsigned index) {
			if(this -> empty() || index < 0 || index >= this -> size()) {
				cout<<"Empty"<<endl;
				return false;
			}
			
			if(index == this->size() - 1) {
				this->pop_back();
			}
			else {
				T elt = (*this)[(*this).size() - 1];
				(*this)[index] = elt;
				(*this).pop_back();
			}
			return true;
		}
	
	public:
		
		bool remove(const T& element) {
			unsigned index = 0;
			for(auto&& item : (*this)) {
				if(item == element) {
					return remove_elt(index);
				}
				index++;
			}
			return false;
		}
		
		bool find(const T &value) {
			for(auto it = this -> begin(); it != this -> end() ; it ++) {
				if((*it) == value) {
					return true;
				}
			}
			return false;
		}
		
		int find_min() {
			if(!this -> empty()) {
				auto min = (*this)[0];
				auto index = 0;
				auto iter = 0;
				for(auto&& item : (*this)) {
					if(item < min) {
						min = item;
						index = iter;
					}
					iter++;
				}
				return index;
			}
			return -1;
		}
		
		void print_vector() {
			for(auto it = (*this).begin() ; it != (*this).end() ; it++) {
				(*it).peek_info();
			}
			cout<<endl;
		}
		
		bool check_value(double value) {
			bool ret = true;
			for(auto&& item : (*this)) {
				if(item <= value) {
					ret = false;
					break;
				}
			}
			
			return ret;
		}
}; 

inline double rate_function(double x, double y) {
	
	return (rate_func == rf1) ? sqrt(x*y/(x+y)) : sqrt(x+y);
	
}

//inline int find_level(double curr_time, int ratio, double small_tau, int Step) {
//	return int((curr_time - ratio*small_tau*Step)/small_tau);
//}

template <class T, class P>
void big_step_distribute(T &clock_time_in_step, P &time_array, const int N, const double small_tau, const int ratio, const int Step) {
	int tmp;
	for(auto&& item : time_array) {
		//Here is problem
		if(item.time >= (Step + 1)*ratio*small_tau)
		{
			tmp = ratio;
		}
		else
		{
			tmp = int((item.time - ratio*small_tau*Step)/small_tau);
		}
		
		clock_time_in_step[tmp].push_back(item);
	}
}

template <class T>
void move_interaction(T &clock_time_in_step, interaction *pt, const double small_tau, const int ratio, const int Step, const double new_time, const int level, const string invoke)
//move the interaction pointed by *pt from old bucket to new bucket
//n_move1: move without relinking pointers
//n_move2: move with relinking pointers
{
	//printf("(Move_Interaction) current Step: %d\n", Step);
//	cout<<"(Move_interaction) Time: "<<pt.time<<" Location: "<<pt.location<<endl;
	double old_time = pt->time;
	int old_level, new_level;
	
	if (old_time >= (Step + 1)*ratio*small_tau )
	{
		old_level = ratio;
	}
	else
	{
		old_level = int( (old_time - Step*ratio*small_tau)/small_tau );
	}
	//adjcent clocks might not be in the level where the min_loc clock is.
	if ( new_time >= (Step + 1)*ratio*small_tau )
	{
		new_level = ratio;
	}
	else
	{
		new_level = int( (new_time - Step*ratio*small_tau)/small_tau );
		if( new_level < 0 || new_level >= ratio)
		{
			printf("Error! Moving %d from %d to %d ; new time is %.5lf ; current step is %d*********\n", pt->location, old_level, new_level, new_time, Step);
			if(old_level >= 0 && old_level <= ratio) {
				//Check if pt->loction is in old_level
				cout<<"Bad!!!!!!\n";
//				for(interaction* tmp_pt = clock_time_in_step[old_level] ; tmp_pt != NULL ; tmp_pt = tmp_pt -> right) {
//					if(tmp_pt -> location == pt -> location) {
//						cout<<"Found!********** \n";
//						if(old_level < level) {
//							cout<<"Uh Oh! Not good*********\n";
//						}
//						break;
//					}
//				}
			}
			//Bug1: Step is negative number;
		}
	}
	
//	cout<<"(Move_interaction) Old_Level: "<<old_level<<" New_Level: "<<new_level<<endl;
	bool found = clock_time_in_step[old_level].remove((*pt));
	if(!found) {
		cout<<invoke<<endl;
		printf("Check Result: %.5lf, Location: %d\n", pt->time, pt->location);
		printf("Current Step: %d, Boundary: %.5lf\n", Step, Step * small_tau * ratio);
		printf("Value before taking flooring function: %.5lf\n", (pt->time - Step * small_tau * ratio)/small_tau);
		//cout<<<<pt->time<<", "<<pt->location<<endl;
		printf("Old Time: %.5lf, New Time: %.5lf, pt->time: %.5lf\n", old_time, new_time, pt->time);
		cout<<"Level: "<<old_level<<" Size: "<<clock_time_in_step[old_level].size()<<"************************\n";
		//clock_time_in_step[old_level].print_vector();
		
		cout<<"All buckets Info\n";
		for(int i = 0 ; i < ratio + 1 ; i++) {
			//cout<<"bucket: "<<i<<" ("<<(i + 1) * Step * small_tau * ratio<<")"<<endl;
			printf("bucket: %d (%.5lf)", i , Step * small_tau * ratio + (i * small_tau));
			cout<<"------------"<<endl;
			clock_time_in_step[i].print_vector();
			cout<<endl;
		}
		
		cout<<"*************************";
		exit(0);
	}
	pt->time = new_time;
	clock_time_in_step[new_level].push_back((*pt));
}

template <class T, class P, class Q>
void update(T &clock_time_in_step, P &time_array, Q &energy_array, const int N, const double small_tau, const int ratio, const int Step, trng::uniform01_dist<> &u, trng::yarn2 &r, int &count, const int level, double &energy_integration) {
//update clock_time_in_step[level]
	//cout<<"****************Level: "<<level<<"****************"<<endl;
	int min_index_in_current_bucket = clock_time_in_step[level].find_min();
	int min_loc = clock_time_in_step[level].at(min_index_in_current_bucket).location;
	interaction *pointer = &time_array[min_loc];
	double current_time = pointer->time;
	double previous_energy = 0.0;
	double current_energy = 0.0;
	//int a = 0;
	while(true) {//(current_time < next_time)
		
		count++;
		
		//Step 1: update min interaction and energy
		double total_energy = energy_array[min_loc] + energy_array[min_loc + 1];
		double tmp_double;
		double tmp_rate;
		double rnd;
		do{
			rnd = u(r);
		}while(rnd < 1e-16 || rnd > 1 -  1e-16);
		tmp_rate = rate_function(energy_array[min_loc], energy_array[min_loc + 1]);
		if(tmp_rate < 1e-16|| std::isnan(tmp_rate))
	  	{
	    	tmp_rate = 1e-16;
	    	cout<<"energy ="<< energy_array[min_loc]<<" "<<energy_array[min_loc + 1]<<endl;
		exit(0);
	  	}
		tmp_double = - log(rnd)/tmp_rate;

		//        cout<<"added time = " << tmp_double << endl;
		double old_e_left = energy_array[min_loc];
		double old_e_right = energy_array[min_loc + 1];
		double tmp_rnd_uni;
		do{
		  tmp_rnd_uni = u(r);
		}while(tmp_rnd_uni < 1e-16 || tmp_rnd_uni > 1 - 1e-16);
		
		previous_energy = energy_array[min_loc];
		if(min_loc == 0)
		{
			previous_energy = energy_array[min_loc + 1];
			total_energy = old_e_right  -log(1 - u(r))*TL;
		}
		if(min_loc == N)
		{
			total_energy = old_e_left  -log(1 - u(r))*TR;
		}
		if(min_loc != 0)
		{
			energy_array[min_loc] = tmp_rnd_uni*total_energy;
			//            E_avg[min_loc - 1] += (current_time - last_update[min_loc - 1])*old_e_left;
			//            last_update[min_loc - 1] = current_time;
		}
		if(min_loc != N)
		{
			energy_array[min_loc + 1] = (1 - tmp_rnd_uni)*total_energy;//update energy
			//            E_avg[min_loc] += (current_time - last_update[min_loc])*old_e_right;
			//            last_update[min_loc] = current_time;
		}
		move_interaction(clock_time_in_step, pointer, small_tau,ratio, Step,  current_time + tmp_double, level, "Min_loc Call");
		//clock_time_in_step[level].print_queue();
		if(time_array[min_loc].time != current_time + tmp_double) {
			cout<<"Not Updated Properly!"<<endl;
			cout<<"==================\n";
			exit(0);
		}
		current_energy = energy_array[min_loc];
		if(min_loc == 0) {
			current_energy = energy_array[min_loc + 1];
		}
		
		energy_integration += (1 - 2*int(min_loc == 0))*(current_energy - previous_energy);
		//cout<<energy_integration<<endl;
		//Step 2: update other interactions
		if(min_loc != 0)
		{
			pointer = &time_array[min_loc - 1];
			tmp_double = (pointer->time - current_time)*rate_function(energy_array[min_loc - 1], old_e_left)/rate_function(energy_array[min_loc - 1], energy_array[min_loc]) + current_time;
	    		if(tmp_double > 1e16 || tmp_double < current_time || std::isnan(tmp_double))
	     	{
		    //cout<<"tmp_double = "<<tmp_double<<endl;		
	          //cout<<"energy ="<< energy_array[min_loc - 1]<<" "<<energy_array[min_loc]<<endl;
				if(tmp_double < current_time) {
					//printf("Smaller Tmp_double");
					if(pointer->time - current_time < 0) {
						printf("Min_loc != 0; Error! pt->time: %.5lf, current_time: %.5lf*********\n", pointer->time, current_time);
					}
				}
				printf("tmp_double = %.5lf, energy_array[min_loc - 1] = %.5lf, energy_array[min_loc] = %.5lf*********\n",tmp_double, energy_array[min_loc - 1], energy_array[min_loc]);
				exit(0);
	      	}
			if(pointer->time < current_time)
			{
				//cout<<"Error with time"<<pt->time<<"  "<<current_time<<endl;
				printf("Error! pt->time: %.5lf, current_time: %.5lf*********\n", pointer->time, current_time);
				exit(0);
			}
			//tmp_double = (pt->time - current_time)*sqrt(energy_array[min_loc - 1] + old_e_left)/sqrt(energy_array[min_loc - 1] + energy_array[min_loc]) + current_time;
			move_interaction(clock_time_in_step, pointer ,small_tau,ratio, Step, tmp_double, level, "Left Adjcent");
			if(time_array[min_loc - 1].time != tmp_double) {
				cout<<"Not Updated Properly!"<<endl;
				cout<<"==================\n";
				exit(0);
			}
		}
		if(min_loc != N)
		{
			pointer = &time_array[min_loc + 1];
			tmp_double = (pointer->time - current_time)*rate_function(energy_array[min_loc + 2], old_e_right)/rate_function(energy_array[min_loc + 2], energy_array[min_loc + 1]) + current_time;
	   		if(tmp_double > 1e16 || tmp_double < current_time || std::isnan(tmp_double)) {
		    		
				if(tmp_double < current_time) {
					//printf("Smaller Tmp_double");
					if(pointer->time - current_time < 0) {
						printf("Min_loc != N; Error! pt->time: %.5lf, current_time: %.5lf*********\n", pointer->time, current_time);
					}
				}
				
				printf("tmp_double = %.5lf, energy_array[min_loc + 1] = %.5lf, energy_array[min_loc + 2] = %.5lf*********\n", tmp_double, energy_array[min_loc + 1], energy_array[min_loc + 2]);
				exit(0);
			//cout<<"tmp_double = "<<tmp_double<<endl;
		    //cout<<"energy ="<< energy_array[min_loc + 1]<<" "<<energy_array[min_loc + 2]<<endl;
	      	}
			if(pointer->time < current_time)
			{
				printf("Error! pt->time: %.5lf, current_time: %.5lf*********\n", pointer->time, current_time);
				exit(0);
				//cout<<"Error with time"<<pt->time<<"  "<<current_time<<endl;
			}
			//tmp_double = (pt->time - current_time)*sqrt(energy_array[min_loc + 2] + old_e_right)/sqrt(energy_array[min_loc + 2] + energy_array[min_loc + 1]) + current_time;
			move_interaction(clock_time_in_step, pointer,small_tau,ratio, Step, tmp_double, level, "Right Adjcent");
			if(time_array[min_loc + 1].time != tmp_double) {
				cout<<"Not Updated Properly!"<<endl;
				cout<<"==================\n";
				exit(0);
			}
			
		}
		
		//Step 3: update current time
		if(!clock_time_in_step[level].empty())
		{
			min_index_in_current_bucket = clock_time_in_step[level].find_min();
			min_loc = clock_time_in_step[level].at(min_index_in_current_bucket).location;
			//pt = time_array[min_loc];
			pointer = &time_array[min_loc];
			current_time = pointer->time;
			
			for(int i = 0 ; i < level ; i++) {
				if(!clock_time_in_step[i].empty()) {
					cout<<"(Update Func) Oops! Not Correct!\n";
					exit(0);
				}
			}
		}
		else
		{
			break;
		}
		
		
	}
	
	//cout<<"****************End Level: "<<level<<"****************"<<endl<<endl;
	
}

int main(int argc, char *argv[]) {
	
	timeval t1, t2;
	const int N = 60;
	const int Step = 1000000;
	const double big_tau = 0.2;
	const int ratio = int(N/10);
	double small_tau = big_tau/double(ratio);
	const int size_of_time_array = N + 1;
	const int N_thread = 8;
	array<double, N_thread> energy_integration;
	
	if(argc > 1) {
		int func_type = stoi(argv[argc - 1]);
		if(func_type == 1) {
			rate_func = rf1;
			cout<<"Current Rate Function: sqrt(xy/(x + y))"<<endl;
		}
		else {
			cout<<"Current Rate Function: sqrt(x + y)"<<endl;
			rate_func = rf2;
		}
	}
	
	energy_integration.fill(0.0);
	
	array<thread_info, N_thread> threads;
	
	#pragma omp parallel num_threads(N_thread)
	{
		const int rank = omp_get_thread_num();
		
		trng::yarn2 r;
		r.seed(time(NULL));
		trng::uniform01_dist<> u;
		r.split(N_thread, rank);
		
		int count = 0;
		array<double, N + 2> energy_array;
		energy_array[0] = TL;
		energy_array[N+1] = TR;
		for(int n = 1; n < N+1; n++)
			energy_array[n] = 1;
		int n = 0;
		
		array<interaction, size_of_time_array> time_array;

		for(auto&& item : time_array) {
			item.time = -log(1 - u(r))/rate_function(energy_array[n], energy_array[n+1]);
			item.location = n;
			n++;
		}
		
		array<Vector<interaction>, ratio + 1> clock_time_in_step;
		
		gettimeofday(&t1,NULL);
		
		for(int i = 0 ; i < Step ; i++) {
					
			big_step_distribute(clock_time_in_step, time_array, N, small_tau, ratio, i);
			
			bool found = false;
			
			for(int j = 0 ; j < ratio; j ++) {
				for(auto&& item : clock_time_in_step[j]) {
					if(item.time < (i * small_tau * ratio + j * small_tau) || item.time >= (i * small_tau * ratio + (j + 1) * small_tau)) {
						cout<<"Check Big_Step_Distribute()\n";
						printf("Time: %.5lf, Level: %d\n", item.time, j);
						//printf("Check Result: %.5lf, Location: %d\n", pt->time, pt->location);
						printf("Current Step: %d, Boundary: %.5lf\n", i, i * small_tau * ratio);
						//printf("Value before taking flooring function: %.5lf\n", (pt->time - Step * small_tau * ratio)/small_tau);
						//cout<<<<pt->time<<", "<<pt->location<<endl;
						//printf("Old Time: %.5lf, New Time: %.5lf, pt->time: %.5lf\n", old_time, new_time, pt->time);
						//cout<<"Level: "<<old_level<<" Size: "<<clock_time_in_step[old_level].size()<<"************************\n";
						//clock_time_in_step[old_level].print_vector();
						
						cout<<"*******************************\n";
						cout<<"All buckets Info\n";
						for(int k = 0 ; k < ratio + 1 ; k++) {
							//cout<<"bucket: "<<i<<" ("<<(i + 1) * Step * small_tau * ratio<<")"<<endl;
							printf("bucket: %d (%.5lf)", k , i * small_tau * ratio + (k * small_tau));
							cout<<"------------"<<endl;
							clock_time_in_step[k].print_vector();
							cout<<endl;
						}
						
						cout<<"*************************\n\n";
						
						found = true;
						break;
						//exit(0);
					}
				}
				if(found) {
					break;
				}
			}
			
			for(int j = 0 ; j < ratio; j ++) {
				if(!clock_time_in_step[j].empty()) {
					update(clock_time_in_step, time_array, energy_array, N, small_tau, ratio, i, u, r, count, j, energy_integration[rank]);
				}
			}
			
			if(clock_time_in_step[ratio].size() != N + 1) {
				cout<<"No, Check Algorithm\n";
				exit(0);
			}
			bool res = clock_time_in_step[ratio].check_value((i + 1) * small_tau * ratio);

			if(!res) {
				cout<<"Uh Oh, Check Algorithm"<<endl;
				exit(0);
			}
			
			clock_time_in_step[ratio].clear();
			
		}
		
		gettimeofday(&t2, NULL);
		
		double delta = ((t2.tv_sec  - t1.tv_sec) * 1000000u + t2.tv_usec - t1.tv_usec) / 1.e6;
		
		threads[rank].current_rank = rank;
		threads[rank].running_time = 1000000*delta/double(count);
		threads[rank].count = count;
		
	}
	
	
	
	auto rank = 0;
	auto sum = 0.0;
	
	for(auto&& item : energy_integration) {
		cout<<"Rank: "<<rank<<" Flux: "<<item<<endl;
		sum += item;
		rank++;
	}
	cout<<endl;
	cout<<"Mean: "<<sum/(N_thread*Step*big_tau)<<endl<<endl;
	
	
	for(auto&& item : threads) {
		cout<<"Thread: "<<item.current_rank<<" Running Time: "<<item.running_time<<" Count: "<<item.count<<endl;
	}
	
	cout<<endl;
	cout<<"================================================"<<endl;
	cout<<endl;
	
	return 0;
}