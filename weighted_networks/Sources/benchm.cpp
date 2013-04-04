



/*
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                               *
 *	This program is free software; you can redistribute it and/or modify         *
 *  it under the terms of the GNU General Public License as published by         *
 *  the Free Software Foundation; either version 2 of the License, or            *
 *  (at your option) any later version.                                          *
 *                                                                               *
 *  This program is distributed in the hope that it will be useful,              *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of               *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                *
 *  GNU General Public License for more details.                                 *
 *                                                                               *
 *  You should have received a copy of the GNU General Public License            *
 *  along with this program; if not, write to the Free Software                  *
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA    *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                               *
 *  Created by Andrea Lancichinetti on 7/01/09 (email: arg.lanci@gmail.com)      *
 *	Modified on 28/05/09                                                         *
 *	Collaborators: Santo Fortunato												 *
 *  Location: ISI foundation, Turin, Italy                                       *
 *	Project: Benchmarking community detection programs                           *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 */








#include "./standard_include.cpp"
#define unlikely -214741

#include "set_parameters.cpp"




bool they_are_mate(int a, int b, const deque<deque<int> > & member_list) {


	for(int i=0; i<member_list[a].size(); i++) {
		
		if(binary_search(member_list[b].begin(), member_list[b].end(), member_list[a][i]))
			return true;
	
	}

	return false;

}


#include "cc.cpp"






// it computes the sum of a deque<int>

int deque_int_sum(const deque<int> & a) {
	
	int s=0;
	for(int i=0; i<a.size(); i++)
		s+=a[i];

	return s;
}

// it computes the integral of a power law
double integral (double a, double b) {

	
	if (fabs(a+1.)>1e-10)
		return (1./(a+1.)*pow(b, a+1.));
	
	
	else
		return (log(b));

}


// it returns the average degree of a power law
double average_degree(const double &dmax, const double &dmin, const double &gamma) {

	return (1./(integral(gamma, dmax)-integral(gamma, dmin)))*(integral(gamma+1, dmax)-integral(gamma+1, dmin));

}


//bisection method to find the inferior limit, in order to have the expected average degree
double solve_dmin(const double& dmax, const double &dmed, const double &gamma) {
	
	double dmin_l=1;
	double dmin_r=dmax;
	double average_k1=average_degree(dmin_r, dmin_l, gamma);
	double average_k2=dmin_r;
	
	
	if ((average_k1-dmed>0) || (average_k2-dmed<0)) {
		
		cerr<<"\n***********************\nERROR: the average degree is out of range:";
		
		if (average_k1-dmed>0) {
			cerr<<"\nyou should increase the average degree (bigger than "<<average_k1<<")"<<endl; 
			cerr<<"(or decrease the maximum degree...)"<<endl;
		}
		
		if (average_k2-dmed<0) {
			cerr<<"\nyou should decrease the average degree (smaller than "<<average_k2<<")"<<endl; 
			cerr<<"(or increase the maximum degree...)"<<endl;
		}
		
		return -1;	
	}
	
		
	while (fabs(average_k1-dmed)>1e-7) {
		
		double temp=average_degree(dmax, ((dmin_r+dmin_l)/2.), gamma);
		if ((temp-dmed)*(average_k2-dmed)>0) {
			
			average_k2=temp;
			dmin_r=((dmin_r+dmin_l)/2.);
		
		}
		else {
			
			average_k1=temp;
			dmin_l=((dmin_r+dmin_l)/2.);
			
		
		}
			
		

	
	}
	
	return dmin_l;
}



// it computes the correct (i.e. discrete) average of a power law
double integer_average (int n, int min, double tau) {
	
	double a=0;

	for (double h=min; h<n+1; h++)
		a+= pow((1./h),tau);
	
	
	double pf=0;
	for(double i=min; i<n+1; i++)
		pf+=1/a*pow((1./(i)),tau)*i;
	
	return pf;

}



// this function changes the community sizes merging the smallest communities
int change_community_size(deque<int> &seq) {

	
			
	if (seq.size()<=2)
		return -1;
	
	int min1=0;
	int min2=0;
	
	for (int i=0; i<seq.size(); i++)		
		if (seq[i]<=seq[min1])
			min1=i;
	
	if (min1==0)
		min2=1;
	
	for (int i=0; i<seq.size(); i++)		
		if (seq[i]<=seq[min2] && seq[i]>seq[min1])
			min2=i;
	

	
	seq[min1]+=seq[min2];
	
	int c=seq[0];
	seq[0]=seq[min2];
	seq[min2]=c;
	seq.pop_front();
	
	
	return 0;
}








int build_bipartite_network(deque<deque<int> >  & member_matrix, const deque<int> & member_numbers, const deque<int> &num_seq) {

	
	
	// this function builds a bipartite network with num_seq and member_numbers which are the degree sequences. in member matrix links of the communities are stored
	// this means member_matrix has num_seq.size() rows and each row has num_seq[i] elements
	
	
	
	deque<set<int> > en_in;			// this is the Ein of the subgraph
	deque<set<int> > en_out;		// this is the Eout of the subgraph
	
	
	{
		set<int> first;
		for(int i=0; i<member_numbers.size(); i++) {
			en_in.push_back(first);
		}
	}
	
	{
		set<int> first;
		for(int i=0; i<num_seq.size(); i++) {
			en_out.push_back(first);
		}
	}

	
	
	multimap <int, int> degree_node_out;
	deque<pair<int, int> > degree_node_in;
	
	for(int i=0; i<num_seq.size(); i++)
		degree_node_out.insert(make_pair(num_seq[i], i));
	
	for(int i=0; i<member_numbers.size(); i++)
		degree_node_in.push_back(make_pair(member_numbers[i], i));
	
	
	sort(degree_node_in.begin(), degree_node_in.end());
	
	
	
	deque<pair<int, int> >::iterator itlast = degree_node_in.end();
	
	/*
	for (int i=0; i<degree_node_in.size(); i++)
		cout<<degree_node_in[i].first<<" "<<degree_node_in[i].second<<endl;
	*/
	
	
	
	while (itlast != degree_node_in.begin()) {
		
		itlast--;
		
		
		multimap <int, int>::iterator itit= degree_node_out.end();
		deque <multimap<int, int>::iterator> erasenda;
		
		for (int i=0; i<itlast->first; i++) {
			
			if(itit!=degree_node_out.begin()) {
				
				itit--;
				
					
				en_in[itlast->second].insert(itit->second);
				en_out[itit->second].insert(itlast->second);
					
				erasenda.push_back(itit);
				
			}
			
			else
				return -1;
		
		}
		
		
		//cout<<"degree node out before"<<endl;
		//prints(degree_node_out);
		
		for (int i=0; i<erasenda.size(); i++) {
			
			
			if(erasenda[i]->first>1)
				degree_node_out.insert(make_pair(erasenda[i]->first - 1, erasenda[i]->second));
	
			
			degree_node_out.erase(erasenda[i]);

			
		
		}
		
		//cout<<"degree node out after"<<endl;
		//prints(degree_node_out);
		
	}
	
	
	// this is to randomize the subgraph -------------------------------------------------------------------

	
	for(int node_a=0; node_a<num_seq.size(); node_a++) for(int krm=0; krm<en_out[node_a].size(); krm++) {
				
				
		int random_mate=irand(member_numbers.size()-1);
		
				
		if (en_out[node_a].find(random_mate)==en_out[node_a].end()) {
			
			deque <int> external_nodes;
			for (set<int>::iterator it_est=en_out[node_a].begin(); it_est!=en_out[node_a].end(); it_est++)
				external_nodes.push_back(*it_est);
						
										
			int	old_node=external_nodes[irand(external_nodes.size()-1)];
					
			
			deque <int> not_common;
			for (set<int>::iterator it_est=en_in[random_mate].begin(); it_est!=en_in[random_mate].end(); it_est++)
				if (en_in[old_node].find(*it_est)==en_in[old_node].end())
					not_common.push_back(*it_est);
					
			
			if (not_common.empty())
				break;
				
			int node_h=not_common[irand(not_common.size()-1)];
			
			
			en_out[node_a].insert(random_mate);
			en_out[node_a].erase(old_node);
			
			en_in[old_node].insert(node_h);
			en_in[old_node].erase(node_a);
			
			en_in[random_mate].insert(node_a);
			en_in[random_mate].erase(node_h);
			
			en_out[node_h].erase(random_mate);
			en_out[node_h].insert(old_node);
			

		}
	
	
	}

	
	
	member_matrix.clear();
	deque <int> first;
	
	for (int i=0; i<en_out.size(); i++) { 
		
		member_matrix.push_back(first);
		for (set<int>::iterator its=en_out[i].begin(); its!=en_out[i].end(); its++)
			member_matrix[i].push_back(*its);
	
	}
	
	
	return 0;


}



int internal_degree_and_membership (double mixing_parameter, int overlapping_nodes, int max_mem_num, int num_nodes, deque<deque<int> >  & member_matrix, 
bool excess, bool defect,  deque<int> & degree_seq, deque<int> &num_seq, deque<int> &internal_degree_seq, bool fixed_range, int nmin, int nmax, double tau2) {
	
	
	
	
	
	if(num_nodes< overlapping_nodes) {
		
		cerr<<"\n***********************\nERROR: there are more overlapping nodes than nodes in the whole network! Please, decrease the former ones or increase the latter ones"<<endl;
		return -1;
	}
	
	
	// 
	member_matrix.clear();
	internal_degree_seq.clear();
	
	deque<double> cumulative;
	
	// it assigns the internal degree to each node -------------------------------------------------------------------------
	int max_degree_actual=0;		// maximum internal degree

	for (int i=0; i<degree_seq.size(); i++) {
		
		double interno=(1-mixing_parameter)*degree_seq[i];
		int int_interno=int(interno);
		
		
		if (ran4()<(interno-int_interno))
			int_interno++;
		
		if (excess) {
			
			while (   (  double(int_interno)/degree_seq[i] < (1-mixing_parameter) )  &&   (int_interno<degree_seq[i])   )
				int_interno++;
				
		
		}
		
		
		if (defect) {
			
			while (   (  double(int_interno)/degree_seq[i] > (1-mixing_parameter) )  &&   (int_interno>0)   )
				int_interno--;
				
		
		}

		
		
		
		internal_degree_seq.push_back(int_interno);
		
		
		if (int_interno>max_degree_actual)
			max_degree_actual=int_interno;
		
			
	}
	
	
	// it assigns the community size sequence -----------------------------------------------------------------------------
	
	powerlaw(nmax, nmin, tau2, cumulative);
	
	
	if (num_seq.empty()) {
		
		int _num_=0;
		if (!fixed_range && (max_degree_actual+1)>nmin) {
		
			_num_=max_degree_actual+1;			// this helps the assignment of the memberships (it assures that at least one module is big enough to host each node)
			num_seq.push_back(max_degree_actual+1);
		
		}
		
		
		while (true) {
			
			
			int nn=lower_bound(cumulative.begin(), cumulative.end(), ran4())-cumulative.begin()+nmin;
			
			if (nn+_num_<=num_nodes + overlapping_nodes * (max_mem_num-1) ) {
				
				num_seq.push_back(nn);				
				_num_+=nn;
			
			}
			else
				break;
			
			
		}
		
		num_seq[min_element(num_seq.begin(), num_seq.end()) - num_seq.begin()]+=num_nodes + overlapping_nodes * (max_mem_num-1) - _num_;
		
	}
	
	
	//cout<<"num_seq"<<endl;
	//prints(num_seq);
	
	int ncom=num_seq.size();
	
	//cout<<"\n----------------------------------------------------------"<<endl;

	/*
	cout<<"community sizes"<<endl;
	for (int i=0; i<num_seq.size(); i++)
		cout<<num_seq[i]<<" ";
	cout<<endl<<endl;
	//*/
	

	/*
	deque <int> first;
	for (int i=0; i<ncom; i++)
		member_matrix.push_back(first);
	
	
	
	// it puts the overlapping_nodes inside
	cout<<ncom<<endl;
	for (int i=degree_seq.size() - overlapping_nodes; i<degree_seq.size(); i++) {
		
		cout<<i<<endl;
		set<int> members;
		int hh=0;
			
		while(members.size()<max_mem_num) {
				
			int random_module=irand(ncom-1);
				
			if(member_matrix[random_module].size()!=num_seq[random_module])
				members.insert(random_module);
				
			hh++;
				
			if(hh>3*num_nodes) {
				cerr<<"it seems that the overlapping nodes need more communities that those I provided. Please increase the number of communities or decrease the number of overlapping nodes"<<endl;
				return -1;				
			}
	
		}
			
				
		for (set<int>::iterator its=members.begin(); its!=members.end(); its++)
			member_matrix[*its].push_back(i);
				
	}
	
	
	
	// it decides the memberships for the not overlapping nodes		
	
	int moment_module=0;
	for (int i=0; i<num_nodes - overlapping_nodes; i++) {
	
		while(member_matrix[moment_module].size()==num_seq[moment_module])
			 moment_module++;

		member_matrix[moment_module].push_back(i);
		
	}
		
		
	
	*/
	
	// I have to assign the degree to the nodes
	
	
	deque<int> member_numbers;
	for(int i=0; i<overlapping_nodes; i++)
		member_numbers.push_back(max_mem_num);
	for(int i=overlapping_nodes; i<degree_seq.size(); i++)
		member_numbers.push_back(1);
	
	//prints(member_numbers);
	//prints(num_seq);
	
	if(build_bipartite_network(member_matrix, member_numbers, num_seq)==-1) {
		
		cerr<<"it seems that the overlapping nodes need more communities that those I provided. Please increase the number of communities or decrease the number of overlapping nodes"<<endl;
		return -1;			
	
	}
		
	//printm(member_matrix);
	
	//cout<<"degree_seq"<<endl;
	//prints(degree_seq);
	
	//cout<<"internal_degree_seq"<<endl;
	//prints(internal_degree_seq);

	deque<int> available;
	for (int i=0; i<num_nodes; i++)
		available.push_back(0);
	
	for (int i=0; i<member_matrix.size(); i++) {
		for (int j=0; j<member_matrix[i].size(); j++)
			available[member_matrix[i][j]]+=member_matrix[i].size()-1;
	}
	
	//cout<<"available"<<endl;
	//prints(available);
	
	
	deque<int> available_nodes;
	for (int i=0; i<num_nodes; i++)
		available_nodes.push_back(i);
	
	
	deque<int> map_nodes;				// in the position i there is the new name of the node i
	for (int i=0; i<num_nodes; i++)
		map_nodes.push_back(0);

	
	for (int i=degree_seq.size()-1; i>=0; i--) {
		
		int & degree_here=internal_degree_seq[i];
		int try_this = irand(available_nodes.size()-1);
		
		int kr=0;
		while (internal_degree_seq[i] > available[available_nodes[try_this]]) {
		
			kr++;
			try_this = irand(available_nodes.size()-1);
			if(kr==3*num_nodes) {
			
				if(change_community_size(num_seq)==-1) {
					
					cerr<<"\n***********************\nERROR: this program needs more than one community to work fine"<<endl;
					return -1;
				
				}
				
				cout<<"it took too long to decide the memberships; I will try to change the community sizes"<<endl;

				cout<<"new community sizes"<<endl;
				for (int i=0; i<num_seq.size(); i++)
					cout<<num_seq[i]<<" ";
				cout<<endl<<endl;
				
				return (internal_degree_and_membership(mixing_parameter, overlapping_nodes, max_mem_num, num_nodes, member_matrix, excess, defect, degree_seq, num_seq, internal_degree_seq, fixed_range, nmin, nmax, tau2));
			
			
			}
			
			
		}
		
		
		
		map_nodes[available_nodes[try_this]]=i;
		
		available_nodes[try_this]=available_nodes[available_nodes.size()-1];
		available_nodes.pop_back();
		
	
	
	}
	
	
	for (int i=0; i<member_matrix.size(); i++) {
		for (int j=0; j<member_matrix[i].size(); j++)
			member_matrix[i][j]=map_nodes[member_matrix[i][j]];	
	}
	
	
	
	for (int i=0; i<member_matrix.size(); i++)
		sort(member_matrix[i].begin(), member_matrix[i].end());

		
	return 0;

}



int compute_internal_degree_per_node(int d, int m, deque<int> & a) {
	
	
	// d is the internal degree
	// m is the number of memebership 
	
	a.clear();
	int d_i= d/m;
	for (int i=0; i<m; i++)
		a.push_back(d_i);
		
	for(int i=0; i<d%m; i++)
		a[i]++;
	
	
	

	return 0;

}



/*
int check_link_list(const deque<deque<int> > & link_list, const deque<int> & degree_seq) {

	
	for (int i=0; i<link_list.size(); i++) {
	
		int s=0;
		for (int j=0; j<link_list[i].size(); j++)
			s+=link_list[i][j];
		
		if(s!=degree_seq[i]) {
			
			int ok;
			cerr<<"wrong link list"<<endl;
			cin>>ok;
		
		}
		
		
	
	
	}




}

*/


int build_subgraph(deque<set<int> > & E, const deque<int> & nodes, const deque<int> & degrees) {
	
	
	/*
	cout<<"nodes"<<endl;
	prints(nodes);
	
	cout<<"degrees"<<endl;
	prints(degrees);
	*/
	

	
	if(degrees.size()<3) {
		
		cerr<<"it seems that some communities should have only 2 nodes! This does not make much sense (in my opinion) Please change some parameters!"<<endl;
		return -1;
	
	}
	
	// this function is to build a network with the labels stored in nodes and the degree seq in degrees (correspondence is based on the vectorial index)
	// the only complication is that you don't want the npdes to have neighbors they already have
	
	
	
	// labels will be placed in the end
	deque<set<int> > en; // this is the E of the subgraph

	{
		set<int> first;
		for(int i=0; i<nodes.size(); i++) 
			en.push_back(first);
	}
	
	
	
	multimap <int, int> degree_node;
	
	for(int i=0; i<degrees.size(); i++)
		degree_node.insert(degree_node.end(), make_pair(degrees[i], i));
	
	int var=0;

	while (degree_node.size() > 0) {
		
		multimap<int, int>::iterator itlast= degree_node.end();
		itlast--;
		
		multimap <int, int>::iterator itit= itlast;
		deque <multimap<int, int>::iterator> erasenda;
		
		int inserted=0;
		
		for (int i=0; i<itlast->first; i++) {
			
			if(itit!=degree_node.begin()) {
			
				itit--;
				
				
				en[itlast->second].insert(itit->second);
				en[itit->second].insert(itlast->second);
				inserted++;
				
				erasenda.push_back(itit);				
				
			}
			
			else
				break;
		
		}
		
		
		for (int i=0; i<erasenda.size(); i++) {
			
			
			if(erasenda[i]->first>1)
				degree_node.insert(make_pair(erasenda[i]->first - 1, erasenda[i]->second));
	
			degree_node.erase(erasenda[i]);
		
		}

		
		var+= itlast->first - inserted;
		degree_node.erase(itlast);
		
	}

	
	
	// this is to randomize the subgraph -------------------------------------------------------------------
	
	for(int node_a=0; node_a<degrees.size(); node_a++) for(int krm=0; krm<en[node_a].size(); krm++) {
	
					
				
		int random_mate=irand(degrees.size()-1);
		while (random_mate==node_a)
			random_mate=irand(degrees.size()-1);
				
		
		if (en[node_a].insert(random_mate).second) {
			
			deque <int> out_nodes;
			for (set<int>::iterator it_est=en[node_a].begin(); it_est!=en[node_a].end(); it_est++) if ((*it_est)!=random_mate)
				out_nodes.push_back(*it_est);
						
										
					
			int old_node=out_nodes[irand(out_nodes.size()-1)];
					
			en[node_a].erase(old_node);
			en[random_mate].insert(node_a);
			en[old_node].erase(node_a);

										
			deque <int> not_common;
			for (set<int>::iterator it_est=en[random_mate].begin(); it_est!=en[random_mate].end(); it_est++)
				if ((old_node!=(*it_est)) && (en[old_node].find(*it_est)==en[old_node].end()))
					not_common.push_back(*it_est);
					
						
			int node_h=not_common[irand(not_common.size()-1)];
			
			en[random_mate].erase(node_h);
			en[node_h].erase(random_mate);
			en[node_h].insert(old_node);
			en[old_node].insert(node_h);
			
			
		}
		
		
	}

	
	
	// now I try to insert the new links into the already done network. If some multiple links come out, I try to rewire them
	
	deque < pair<int, int> > multiple_edge;
	for (int i=0; i<en.size(); i++) {
		
		for(set<int>::iterator its=en[i].begin(); its!=en[i].end(); its++) if(i<*its ) {
		
			bool already = !(E[nodes[i]].insert(nodes[*its]).second);		// true is the insertion didn't take place
			if (already)
				multiple_edge.push_back(make_pair(nodes[i], nodes[*its]));			
			else
				E[nodes[*its]].insert(nodes[i]);
			
		
		}
	
	
	}
	
	
	for (int i=0; i<multiple_edge.size(); i++) {
		
		
		int &a = multiple_edge[i].first;
		int &b = multiple_edge[i].second;
		
	
		// now, I'll try to rewire this multiple link among the nodes stored in nodes.
		int stopper_ml=0;
		
		while (true) {
					
			stopper_ml++;
			
			int random_mate=nodes[irand(degrees.size()-1)];
			while (random_mate==a || random_mate==b)
				random_mate=nodes[irand(degrees.size()-1)];
			
			if(E[a].find(random_mate)==E[a].end()) {
				
				deque <int> not_common;
				for (set<int>::iterator it_est=E[random_mate].begin(); it_est!=E[random_mate].end(); it_est++)
					if ((b!=(*it_est)) && (E[b].find(*it_est)==E[b].end()) && (binary_search(nodes.begin(), nodes.end(), *it_est)))
						not_common.push_back(*it_est);
				
				if(not_common.size()>0) {
				
					int node_h=not_common[irand(not_common.size()-1)];
					
					
					
					E[random_mate].insert(a);
					E[random_mate].erase(node_h);
					
					E[node_h].erase(random_mate);
					E[node_h].insert(b);
					
					E[b].insert(node_h);
					E[a].insert(random_mate);
					
					break;

				
			
				}
			
			}
			
			if(stopper_ml==2*E.size()) {
	
				cout<<"sorry, I need to change the degree distribution a little bit (one less link)"<<endl;
				break;
	
			}
			
			
			
		}
	
	
	}
	
	


	
	
	return 0;

}





int build_subgraphs(deque<set<int> > & E, const deque<deque<int> > & member_matrix, deque<deque<int> > & member_list, 
	deque<deque<int> > & link_list, const deque<int> & internal_degree_seq, const deque<int> & degree_seq, const bool excess, const bool defect) {
	
	
	
	E.clear();
	member_list.clear();
	link_list.clear();
	
	int num_nodes=degree_seq.size();
	
	//printm(member_matrix);

	
	
	
	{
		
		deque<int> first;
		for (int i=0; i<num_nodes; i++)
			member_list.push_back(first);
	
	}
	
	
	
	for (int i=0; i<member_matrix.size(); i++)
		for (int j=0; j<member_matrix[i].size(); j++)
			member_list[member_matrix[i][j]].push_back(i);
	
		
	for (int i=0; i<member_list.size(); i++) {
		
		deque<int> liin;

		
		for (int j=0; j<member_list[i].size(); j++) {
			
			compute_internal_degree_per_node(internal_degree_seq[i], member_list[i].size(), liin);
			liin.push_back(degree_seq[i] - internal_degree_seq[i]);
		
		}
		
		link_list.push_back(liin);
		
	}
			
	
	

	// now there is the check for the even node (it means that the internal degree of each group has to be even and we want to assure that, otherwise the degree_seq has to change) ----------------------------

	
	
	
	// ------------------------ this is done to check if the sum of the internal degree is an even number. if not, the program will change it in such a way to assure that. 
	
	
			
	for (int i=0; i<member_matrix.size(); i++) {
	
		
		int internal_cluster=0;
		for (int j=0; j<member_matrix[i].size(); j++) {
			
			int right_index= lower_bound(member_list[member_matrix[i][j]].begin(), member_list[member_matrix[i][j]].end(), i) - member_list[member_matrix[i][j]].begin();
			
				
			
			internal_cluster+=link_list[member_matrix[i][j]][right_index];
		}
		
		
		if(internal_cluster % 2 != 0) {
			
			
			
			bool default_flag=false;
			
			if(excess)
				default_flag=true;
			else if(defect)
				default_flag=false;
			else if (ran4()>0.5)
				default_flag=true;
			
			if(default_flag) {
				
				// if this does not work in a reasonable time the degree sequence will be changed
				
				for (int j=0; j<member_matrix[i].size(); j++) {		
					
					
					int random_mate=member_matrix[i][irand(member_matrix[i].size()-1)];
					
					int right_index= lower_bound(member_list[random_mate].begin(), member_list[random_mate].end(), i) - member_list[random_mate].begin();
					
					
					if ((link_list[random_mate][right_index]<member_matrix[i].size()-1) && (link_list[random_mate][link_list[random_mate].size()-1] > 0 )) {
						
						link_list[random_mate][right_index]++;
						link_list[random_mate][link_list[random_mate].size()-1]--;

						break;
						
					}
				
				}
			
			
			}
			
			
			else {
				
				for (int j=0; j<member_matrix[i].size(); j++) {

					int random_mate=member_matrix[i][irand(member_matrix[i].size()-1)];

					int right_index= lower_bound(member_list[random_mate].begin(), member_list[random_mate].end(), i) - member_list[random_mate].begin();
					
					if (link_list[random_mate][right_index] > 0 ) {
						
						link_list[random_mate][right_index]--;
						link_list[random_mate][link_list[random_mate].size()-1]++;

						break;
						
					}
				
				}
			
			
			}
			
			
			
					
		}
			
	
	}
	
	
	// ------------------------ this is done to check if the sum of the internal degree is an even number. if not, the program will change it in such a way to assure that. 
	
	
	
	
	{
	
		set<int> first;
		for(int i=0; i<num_nodes; i++)
			E.push_back(first);
	
	}
	
	for (int i=0; i<member_matrix.size(); i++) {
		
		
		deque<int> internal_degree_i;
		for (int j=0; j<member_matrix[i].size(); j++) {
		
		
			int right_index= lower_bound(member_list[member_matrix[i][j]].begin(), member_list[member_matrix[i][j]].end(), i) - member_list[member_matrix[i][j]].begin();
		
			internal_degree_i.push_back(link_list[member_matrix[i][j]][right_index]);
		
		}
		
		
		
		
		if(build_subgraph(E, member_matrix[i], internal_degree_i)==-1)
			return -1;
	
	
	}




	return 0;
	
}




int connect_all_the_parts(deque<set<int> > & E, const deque<deque<int> > & member_list, const deque<deque<int> > & link_list) {

	
	deque<int> degrees;
	for(int i=0; i<link_list.size(); i++)
		degrees.push_back(link_list[i][link_list[i].size()-1]);
	
	
		
	
	deque<set<int> > en; // this is the en of the subgraph

	{
		set<int> first;
		for(int i=0; i<member_list.size(); i++) 
			en.push_back(first);
	}
	
	
	
	multimap <int, int> degree_node;
	
	for(int i=0; i<degrees.size(); i++)
		degree_node.insert(degree_node.end(), make_pair(degrees[i], i));
	
	int var=0;

	while (degree_node.size() > 0) {
		
		multimap<int, int>::iterator itlast= degree_node.end();
		itlast--;
		
		multimap <int, int>::iterator itit= itlast;
		deque <multimap<int, int>::iterator> erasenda;
		
		int inserted=0;
		
		for (int i=0; i<itlast->first; i++) {
			
			if(itit!=degree_node.begin()) {
			
				itit--;
				
				
				en[itlast->second].insert(itit->second);
				en[itit->second].insert(itlast->second);
				inserted++;
				
				erasenda.push_back(itit);				
				
			}
			
			else
				break;
		
		}
		
		
		for (int i=0; i<erasenda.size(); i++) {
			
			
			if(erasenda[i]->first>1)
				degree_node.insert(make_pair(erasenda[i]->first - 1, erasenda[i]->second));
	
			degree_node.erase(erasenda[i]);
		
		}

		
		var+= itlast->first - inserted;
		degree_node.erase(itlast);
		
	}

		
	// this is to randomize the subgraph -------------------------------------------------------------------
	
	for(int node_a=0; node_a<degrees.size(); node_a++) for(int krm=0; krm<en[node_a].size(); krm++) {
				
		
				
		int random_mate=irand(degrees.size()-1);
		while (random_mate==node_a)
			random_mate=irand(degrees.size()-1);
		
		
		
		if (en[node_a].insert(random_mate).second) {
			
			deque <int> out_nodes;
			for (set<int>::iterator it_est=en[node_a].begin(); it_est!=en[node_a].end(); it_est++) if ((*it_est)!=random_mate)
				out_nodes.push_back(*it_est);
						
										
					
			int old_node=out_nodes[irand(out_nodes.size()-1)];
					
			en[node_a].erase(old_node);
			en[random_mate].insert(node_a);
			en[old_node].erase(node_a);

										
			deque <int> not_common;
			for (set<int>::iterator it_est=en[random_mate].begin(); it_est!=en[random_mate].end(); it_est++)
				if ((old_node!=(*it_est)) && (en[old_node].find(*it_est)==en[old_node].end()))
					not_common.push_back(*it_est);
					
						
			int node_h=not_common[irand(not_common.size()-1)];
			
			en[random_mate].erase(node_h);
			en[node_h].erase(random_mate);
			en[node_h].insert(old_node);
			en[old_node].insert(node_h);
			
			
		}
	
	
	

	}

	
	
	
	// now there is a rewiring process to avoid "mate nodes" (nodes with al least one membership in common) to link each other
	
	int var_mate=0;
	for(int i=0; i<degrees.size(); i++) for(set<int>::iterator itss= en[i].begin(); itss!=en[i].end(); itss++) if(they_are_mate(i, *itss, member_list)) {
		var_mate++;
	}
	
	//cout<<"var mate = "<<var_mate<<endl;
	
	int stopper_mate=0;
	int mate_trooper=10;
	
	while(var_mate>0) {
	
		
		//cout<<"var mate = "<<var_mate<<endl;

		
		int best_var_mate=var_mate;
	
		// ************************************************  rewiring
		
		
		for(int a=0; a<degrees.size(); a++) for(set<int>::iterator its= en[a].begin(); its!=en[a].end(); its++) if(they_are_mate(a, *its, member_list)) {
				
			
			
			int b=*its;
			int stopper_m=0;
			
			while (true) {
						
				stopper_m++;
				
				int random_mate=irand(degrees.size()-1);
				while (random_mate==a || random_mate==b)
					random_mate=irand(degrees.size()-1);
				
				
				if(!(they_are_mate(a, random_mate, member_list)) && (en[a].find(random_mate)==en[a].end())) {
					
					deque <int> not_common;
					for (set<int>::iterator it_est=en[random_mate].begin(); it_est!=en[random_mate].end(); it_est++)
						if ((b!=(*it_est)) && (en[b].find(*it_est)==en[b].end()))
							not_common.push_back(*it_est);
					
					if(not_common.size()>0) {
					
						int node_h=not_common[irand(not_common.size()-1)];
						
						
						en[random_mate].erase(node_h);
						en[random_mate].insert(a);
						
						en[node_h].erase(random_mate);
						en[node_h].insert(b);
						
						en[b].erase(a);
						en[b].insert(node_h);
						
						en[a].insert(random_mate);
						en[a].erase(b);
						
											
						
						if(!they_are_mate(b, node_h, member_list))
							var_mate-=2;
						
						
						if(they_are_mate(random_mate, node_h, member_list))
							var_mate-=2;
						
						break;

					
				
					}
				
				}
				
				if(stopper_m==en[a].size())
					break;
				
				
				
			}
				
			
			break;		// this break is done because if you erased some link you have to stop this loop (en[i] changed)
	
	
		}

		// ************************************************  rewiring
		
		
				

		if(var_mate==best_var_mate) {
			
			stopper_mate++;
			
			if(stopper_mate==mate_trooper)
				break;

		}
		else
			stopper_mate=0;
		
		
		
		//cout<<"var mate = "<<var_mate<<endl;

	
	}
	
	
	
	//cout<<"var mate = "<<var_mate<<endl;

	for (int i=0; i<en.size(); i++) {
		
		for(set<int>::iterator its=en[i].begin(); its!=en[i].end(); its++) if(i<*its) {
		
			E[i].insert(*its);
			E[*its].insert(i);
			
		
		}
	
	
	}
	
	
	
	return 0;

}

int internal_kin(deque<set<int> > & E, const deque<deque<int> > & member_list, int i) {
	
	int var_mate2=0;
	for(set<int>::iterator itss= E[i].begin(); itss!=E[i].end(); itss++) if(they_are_mate(i, *itss, member_list)) 
		var_mate2++;	

	return var_mate2;
	
}




int internal_kin_only_one(set<int> & E, const deque<int> & member_matrix_j) {		// return the overlap between E and member_matrix_j
	
	int var_mate2=0;
	
	for(set<int>::iterator itss= E.begin(); itss!=E.end(); itss++) {
	
		if(binary_search(member_matrix_j.begin(), member_matrix_j.end(), *itss))
			var_mate2++;
	
	}
	
	return var_mate2;
	
}

int erase_links(deque<set<int> > & E, const deque<deque<int> > & member_list, const bool excess, const bool defect, const double mixing_parameter) {

	
	int num_nodes= member_list.size();
	
	int eras_add_times=0;
	
	if (excess) {
		
		for (int i=0; i<num_nodes; i++) {
			
			
			while ( (E[i].size()>1) &&  double(internal_kin(E, member_list, i))/E[i].size() < 1 - mixing_parameter) {
			
			//---------------------------------------------------------------------------------
				
				
				cout<<"degree sequence changed to respect the option -sup ... "<<++eras_add_times<<endl;
				
				deque<int> deqar;
				for (set<int>::iterator it_est=E[i].begin(); it_est!=E[i].end(); it_est++)
					if (!they_are_mate(i, *it_est, member_list))
						deqar.push_back(*it_est);
				
				
				if(deqar.size()==E[i].size()) {	// this shouldn't happen...
				
					cerr<<"sorry, something went wrong: there is a node which does not respect the constraints. (option -sup)"<<endl;
					return -1;
				
				}
				
				int random_mate=deqar[irand(deqar.size()-1)];
				
				E[i].erase(random_mate);
				E[random_mate].erase(i);
				
		
			}
		}
	
	}
	
	
	
	if (defect) {
			
		for (int i=0; i<num_nodes; i++)
			while ( (E[i].size()<E.size()) &&  double(internal_kin(E, member_list, i))/E[i].size() > 1 - mixing_parameter) {
				
				//---------------------------------------------------------------------------------
					
				
				cout<<"degree sequence changed to respect the option -inf ... "<<++eras_add_times<<endl;


				int stopper_here=num_nodes;
				int stopper_=0;
				
				int random_mate=irand(num_nodes-1);
				while ( (    (they_are_mate(i, random_mate, member_list)) || E[i].find(random_mate)!=E[i].end())      &&      (stopper_<stopper_here) ) {
					
					random_mate=irand(num_nodes-1);
					stopper_++;
				
				
				}
				
				if(stopper_==stopper_here) {	// this shouldn't happen...
				
					cerr<<"sorry, something went wrong: there is a node which does not respect the constraints. (option -inf)"<<endl;
					return -1;
				
				}
				
				
				
				E[i].insert(random_mate);
				E[random_mate].insert(i);
				
								
		
			}
			
		
	}

	//------------------------------------ Erasing links   ------------------------------------------------------

	


	return 0;
	
}


int print_network(deque<set<int> > & E, const deque<deque<int> > & member_list, const deque<deque<int> > & member_matrix, 
	deque<int> & num_seq, deque<map <int, double > > & neigh_weigh, double beta, double mu, double mu0) {

	
	int edges=0;

	
	int num_nodes=member_list.size();
	
	deque<double> double_mixing;
	for (int i=0; i<E.size(); i++) {
		
		double one_minus_mu = double(internal_kin(E, member_list, i))/E[i].size();
		
		double_mixing.push_back(1.- one_minus_mu);
				
		edges+=E[i].size();
		
	}
	
	
	//cout<<"\n----------------------------------------------------------"<<endl;
	//cout<<endl;
	
		
	double density=0; 
	double sparsity=0;
	
	for (int i=0; i<member_matrix.size(); i++) {

		double media_int=0;
		double media_est=0;
		
		for (int j=0; j<member_matrix[i].size(); j++) {
			
			
			double kinj = double(internal_kin_only_one(E[member_matrix[i][j]], member_matrix[i]));
			media_int+= kinj;
			media_est+=E[member_matrix[i][j]].size() - double(internal_kin_only_one(E[member_matrix[i][j]], member_matrix[i]));
					
		}
		
		double pair_num=(member_matrix[i].size()*(member_matrix[i].size()-1));
		double pair_num_e=((num_nodes-member_matrix[i].size())*(member_matrix[i].size()));
		
		if(pair_num!=0)
			density+=media_int/pair_num;
		if(pair_num_e!=0)
			sparsity+=media_est/pair_num_e;
		
		
	
	}
	
	density=density/member_matrix.size();
	sparsity=sparsity/member_matrix.size();
	
	
	


	ofstream out1("network.dat");
	for (int u=0; u<E.size(); u++) {

		set<int>::iterator itb=E[u].begin();
	
		while (itb!=E[u].end())
			out1<<u+1<<"\t"<<*(itb++)+1<<"\t"<<neigh_weigh[u][*(itb)]<<endl;
		
		

	}
		

	
	ofstream out2("community.dat");

	for (int i=0; i<member_list.size(); i++) {
		
		out2<<i+1<<"\t";
		for (int j=0; j<member_list[i].size(); j++)
			out2<<member_list[i][j]+1<<" ";
		out2<<endl;
	
	}

	cout<<"\n\n---------------------------------------------------------------------------"<<endl;
	
	
	cout<<"network of "<<num_nodes<<" vertices and "<<edges/2<<" edges"<<";\t average degree = "<<double(edges)/num_nodes<<endl;
	cout<<"\naverage mixing parameter (topology): "<< average_func(double_mixing)<<" +/- "<<sqrt(variance_func(double_mixing))<<endl;
	cout<<"p_in: "<<density<<"\tp_out: "<<sparsity<<endl;

	
	
	ofstream statout("statistics.dat");
	
	deque<int> degree_seq;
	for (int i=0; i<E.size(); i++)
		degree_seq.push_back(E[i].size());
	
	statout<<"degree distribution (probability density function of the degree in logarithmic bins) "<<endl;
	log_histogram(degree_seq, statout, 10);
	statout<<"\ndegree distribution (degree-occurrences) "<<endl;
	int_histogram(degree_seq, statout);
	statout<<endl<<"--------------------------------------"<<endl;

		
	statout<<"community distribution (size-occurrences)"<<endl;
	int_histogram(num_seq, statout);
	statout<<endl<<"--------------------------------------"<<endl;

	statout<<"mixing parameter (topology)"<<endl;
	not_norm_histogram(double_mixing, statout, 20, 0, 0);
	statout<<endl<<"--------------------------------------"<<endl;
	
	
	//*
	
	deque<double> inwij;
	deque<double> outwij;
	//deque<double> inkij;
	//deque<double> outkij;
	
	double csi=(1. - mu) / (1. - mu0);
	double csi2=mu /mu0;
	
	
	double tstrength=0;
	deque<double> one_minus_mu2;
	
	for(int i=0; i<neigh_weigh.size(); i++) {
		
		
		double internal_strength_i=0;
		double strength_i=0;
		
		for(map<int, double>::iterator itm = neigh_weigh[i].begin(); itm!=neigh_weigh[i].end(); itm++) {
			
			
			if(they_are_mate(i, itm->first, member_list)) {
				
				inwij.push_back(itm->second);
				//inkij.push_back(csi * pow(E[i].size(), beta-1));
				internal_strength_i+=itm->second;

				
			}
			else {
				
				outwij.push_back(itm->second);
				//outkij.push_back(csi2 * pow(E[i].size(), beta-1));
			
			
			}
			
			tstrength+=itm->second;
			strength_i+=itm->second;
		
		
		}
		
		one_minus_mu2.push_back(1 - internal_strength_i/strength_i);
		
	}
	
	
	//cout<<"average strength "<<tstrength / E.size()<<"\taverage internal strenght: "<<average_internal_strenght<<endl;
	cout<<"\naverage mixing parameter (weights): "<<average_func(one_minus_mu2)<<" +/- "<<sqrt(variance_func(one_minus_mu2))<<endl;	
	statout<<"mixing parameter (weights)"<<endl;
	not_norm_histogram(one_minus_mu2, statout, 20, 0, 0);
	statout<<endl<<"--------------------------------------"<<endl;

	
	//cout<<" expected internal "<<tstrength * (1 - mu) / E.size()<<endl;
	//cout<<"internal links: "<<inwij.size()<<" external: "<<outwij.size()<<endl;
	
	/*
	ofstream hout1("inwij.dat");
	not_norm_histogram(inwij, hout1, 20, 0, 0);
	ofstream hout2("outwij.dat");
	not_norm_histogram(outwij, hout2, 20, 0, 0);
	ofstream hout3("corrin.dat");
	not_norm_histogram_correlated(inkij, inwij, hout3, 20, 0, 0);
	ofstream hout4("corrout.dat");
	not_norm_histogram_correlated(outkij, outwij, hout4, 20, 0, 0);
	
	//*/
	
	//*/
	
	cout<<"average weight of an internal link "<<average_func(inwij)<<" +/- "<<sqrt(variance_func(inwij))<<endl;
	cout<<"average weight of an external link "<<average_func(outwij)<<" +/- "<<sqrt(variance_func(outwij))<<endl;


	//cout<<"average weight of an internal link expected "<<tstrength / edges * (1. - mu) / (1. - mu0)<<endl;
	//cout<<"average weight of an external link expected "<<tstrength / edges * (mu) / (mu0)<<endl;
	
	
	statout<<"internal weights (weight-occurrences)"<<endl;
	not_norm_histogram(inwij, statout, 20, 0, 0);
	statout<<endl<<"--------------------------------------"<<endl;
	
	
	statout<<"external weights (weight-occurrences)"<<endl;
	not_norm_histogram(outwij, statout, 20, 0, 0);

	


	cout<<endl<<endl;

	return 0;

}



//*
int check_weights(deque<map <int, double > > & neigh_weigh, const deque<deque<int> > & member_list, 
	deque<deque<double> > & wished, deque<deque<double> > & factual, const double tot_var, double *strs) {

	
	
	
	double d1t=0;
	double d2t=0;
	double d3t=0;
	
	double var_check=0;
	
	for (int k=0; k<member_list.size(); k++) {
	
		
		double in_s=0;
		double out_s=0;
		
		
		
			
		for(map<int, double>::iterator itm=neigh_weigh[k].begin(); itm!=neigh_weigh[k].end(); itm++) {
			
			if(itm->second<0)
				cherr(333);
			
			if(fabs(itm->second - neigh_weigh[itm->first][k]) > 1e-7)
				cherr(111);
			
			if (they_are_mate(k, itm->first, member_list))
				in_s+=itm->second;
			else
				out_s+=itm->second;
		
		}
		
		
					
		if (fabs(in_s - factual[k][0]) > 1e-7)
			cherr(in_s - factual[k][0]);
		//cout<<"ok1"<<endl;
		
		if (fabs(out_s - factual[k][1]) > 1e-7)
			cherr(out_s - factual[k][1]);		
		//cout<<"ok2"<<endl;


		if (fabs(in_s + out_s + factual[k][2] - strs[k]) > 1e-7)
			cherr(in_s + out_s + factual[k][2] - strs[k]);
		//cout<<"ok3"<<endl;
		
		double d1=(in_s - wished[k][0]);
		double d2=(out_s - wished[k][1]);
		double d3=(strs[k] - in_s - out_s);
		var_check+= d1*d1 + d2*d2 + d3*d3;
		
		d1t+=d1*d1;
		d2t+=d2*d2;
		d3t+=d3*d3;
		
		
	
	}
	
		
	cout<<"tot_var "<<tot_var<<"\td1t "<<d1t<<"\td2t "<<d2t<<"\td3t "<<d3t<<endl;
	if (fabs(var_check - tot_var) > 1e-5)
		cerr<<"found this difference in check "<<fabs(var_check - tot_var)<<endl;
	else
		cout<<"ok: check passed"<<endl;

	// ----------------------------------------------------------------------- CHECK -----------------------------------------------
	// ----------------------------------------------------------------------- CHECK -----------------------------------------------
	
	


	return 0;

}


//*/



int propagate(deque<map <int, double > > & neigh_weigh, const deque<deque<int> > & member_list, 
	deque<deque<double> > & wished, deque<deque<double> > & factual, int i, double & tot_var, double *strs, const deque<int> & internal_kin_top) {
	
	
	
	
		
	
	
	{		// in this case I rewire the idle strength

		
		double change=factual[i][2]/neigh_weigh[i].size();
		
		
		
		
		
		double oldpartvar=0;
		for(map<int, double>::iterator itm=neigh_weigh[i].begin(); itm!=neigh_weigh[i].end(); itm++) if(itm->second + change > 0)
			for (int bw=0; bw<3; bw++) 
				oldpartvar+= (factual[itm->first][bw] - wished[itm->first][bw]) * (factual[itm->first][bw] - wished[itm->first][bw]);
	
	
		for (int bw=0; bw<3; bw++)
			oldpartvar+= (factual[i][bw] - wished[i][bw]) * (factual[i][bw] - wished[i][bw]);

		
		double newpartvar=0;

		for(map<int, double>::iterator itm=neigh_weigh[i].begin(); itm!=neigh_weigh[i].end(); itm++) if(itm->second + change > 0) {
		
			
			if (they_are_mate(i, itm->first, member_list)) {
					
				factual[itm->first][0]+=change;
				factual[itm->first][2]-=change;
					
				factual[i][0]+=change;
				factual[i][2]-=change;
					
			}
				
			else {
					
				factual[itm->first][1]+=change;
				factual[itm->first][2]-=change;
					
				factual[i][1]+=change;
				factual[i][2]-=change;

					
			}
			
			for (int bw=0; bw<3; bw++)
				newpartvar+= (factual[itm->first][bw] - wished[itm->first][bw]) * (factual[itm->first][bw] - wished[itm->first][bw]);
			
			
			
			itm->second+= change;
			neigh_weigh[itm->first][i]+=change;
			
	
		}
		
	
	
		for (int bw=0; bw<3; bw++)
			newpartvar+= (factual[i][bw] - wished[i][bw]) * (factual[i][bw] - wished[i][bw]);
		
		
		
		tot_var+= newpartvar - oldpartvar;
			
	
		
		
			
	
	}
	
	
	int internal_neigh=internal_kin_top[i];


	if(internal_neigh!=0) {		// in this case I rewire the difference strength

		
		
		double changenn=(factual[i][0] - wished[i][0]);
				
		
		double oldpartvar=0;
		for(map<int, double>::iterator itm=neigh_weigh[i].begin(); itm!=neigh_weigh[i].end(); itm++) {
				
			
			if(they_are_mate(i, itm->first, member_list)) {
				
				double change = changenn/internal_neigh;
				
				if(itm->second - change > 0)
					for (int bw=0; bw<3; bw++) 
						oldpartvar+= (factual[itm->first][bw] - wished[itm->first][bw]) * (factual[itm->first][bw] - wished[itm->first][bw]);				
			
			} 
			
			else {
				
				double change = changenn/(neigh_weigh[i].size() - internal_neigh);

				
				if(itm->second + change > 0)
					for (int bw=0; bw<3; bw++) 
						oldpartvar+= (factual[itm->first][bw] - wished[itm->first][bw]) * (factual[itm->first][bw] - wished[itm->first][bw]);

			}
		
		}
		
	
		for (int bw=0; bw<3; bw++)
			oldpartvar+= (factual[i][bw] - wished[i][bw]) * (factual[i][bw] - wished[i][bw]);

		
		double newpartvar=0;

		for(map<int, double>::iterator itm=neigh_weigh[i].begin(); itm!=neigh_weigh[i].end(); itm++) {
		
			
			if (they_are_mate(i, itm->first, member_list)) {
				
				double change = changenn/internal_neigh;

				
				if(itm->second - change > 0) {
					
					factual[itm->first][0]-=change;
					factual[itm->first][2]+=change;
					
					factual[i][0]-=change;
					factual[i][2]+=change;
					
					for (int bw=0; bw<3; bw++)
						newpartvar+= (factual[itm->first][bw] - wished[itm->first][bw]) * (factual[itm->first][bw] - wished[itm->first][bw]);
					
					
					itm->second-= change;
					neigh_weigh[itm->first][i]-=change;


				
				}
					
			}
				
			else {
				
				
				double change = changenn/(neigh_weigh[i].size() - internal_neigh);
				
				if(itm->second + change > 0) {
					
					factual[itm->first][1]+=change;
					factual[itm->first][2]-=change;
					
					factual[i][1]+=change;
					factual[i][2]-=change;
					
					for (int bw=0; bw<3; bw++)
						newpartvar+= (factual[itm->first][bw] - wished[itm->first][bw]) * (factual[itm->first][bw] - wished[itm->first][bw]);
					
					
					itm->second+= change;
					neigh_weigh[itm->first][i]+=change;


					
				}
					
			}
			
							
	
		}
		
	
	
		for (int bw=0; bw<3; bw++)
			newpartvar+= (factual[i][bw] - wished[i][bw]) * (factual[i][bw] - wished[i][bw]);
		
		
		tot_var+=newpartvar - oldpartvar;
			
		
				

	
	
	
	}
	
	
	//check_weights(neigh_weigh, member_list, wished, factual, tot_var, strs);
	
	
	return 0;
	
}
	
	

int weights(deque<set<int> > & en, const deque<deque<int> > & member_list, const double beta, const double mu, deque<map <int, double > > & neigh_weigh) {

	
	
	double tstrength=0;
	
	
	for(int i=0; i<en.size(); i++)
		tstrength+=pow(en[i].size(), beta);
		
	
	
	double strs[en.size()]; // strength of the nodes
	// build a matrix like this: deque < map <int, double > > each row corresponds to link - weights 

	
	
	for(int i=0; i<en.size(); i++) {
		
		map<int, double> new_map;
		neigh_weigh.push_back(new_map);
		
		for (set<int>::iterator its=en[i].begin(); its!=en[i].end(); its++)
			neigh_weigh[i].insert(make_pair(*its, 0.));
		
		strs[i]=pow(double(en[i].size()), beta);
		
	}
	
	
	
	
	deque<double> s_in_out_id_row(3);
	s_in_out_id_row[0]=0;
	s_in_out_id_row[1]=0;
	s_in_out_id_row[2]=0;
	
	
	deque<deque<double> > wished;	// 3 numbers for each node: internal, idle and extra strength. the sum of the three is strs[i]. wished is the theoretical, factual the factual one.
	deque<deque<double> > factual;
	
	
	
	for (int i=0; i<en.size(); i++) {
		
		wished.push_back(s_in_out_id_row);
		factual.push_back(s_in_out_id_row);
		
	}
	
	
	
	double tot_var=0;
	

	for (int i=0; i<en.size(); i++) {
		
		wished[i][0]=(1. -mu)*strs[i];
		wished[i][1]=mu*strs[i];
		
		factual[i][2]=strs[i];
		
		tot_var+= wished[i][0] * wished[i][0] + wished[i][1] * wished[i][1] + strs[i] * strs[i];
	
	}
	
	
	deque<int> internal_kin_top;
	for(int i=0; i<en.size(); i++)
		internal_kin_top.push_back(internal_kin(en, member_list, i));

	
	
	double precision = 1e-9;
	double precision2 = 1e-2;
	double not_better_than = pow(tstrength, 2) * precision;
	//cout<<"tot_var "<<tot_var<<";\tnot better "<<not_better_than<<endl;
	
	
	
	int step=0;
	
	while (true) {
	
		
		

		time_t t0=time(NULL);
		
		double pre_var=tot_var;
		
		
		for (int i=0; i<en.size(); i++) 
			propagate(neigh_weigh, member_list, wished, factual, i, tot_var, strs, internal_kin_top);
		
		
		//check_weights(neigh_weigh, member_list, wished, factual, tot_var, strs);
		
		double relative_improvement=double(pre_var - tot_var)/pre_var;		
		//cout<<"tot_var "<<tot_var<<"\trelative improvement: "<<relative_improvement<<endl;

		if (tot_var<not_better_than)
				break;
		
		if (relative_improvement < precision2)
			break;
		
		time_t t1= time(NULL);
		int deltat= t1 - t0;
		
		/*
		if(step%2==0 && deltat !=0)		
			cout<<"About "<<cast_int((log(not_better_than) -log(tot_var)) / log(1. - relative_improvement)) * deltat<<" secs..."<<endl;
		*/
		step++;
		
	}
	
	
	//check_weights(neigh_weigh, member_list, wished, factual, tot_var, strs);

	

	
	

	
	
	return 0;

}


int benchmark(bool excess, bool defect, int num_nodes, double  average_k, int  max_degree, double  tau, double  tau2, 
	double  mixing_parameter, double  mixing_parameter2, double  beta, int  overlapping_nodes, int  overlap_membership, int  nmin, int  nmax, bool  fixed_range, double ca) {	

	
	
	
	

	double dmin=solve_dmin(max_degree, average_k, -tau);
	if (dmin==-1)
		return -1;
	
	int min_degree=int(dmin);
	
	
	double media1=integer_average(max_degree, min_degree, tau);
	double media2=integer_average(max_degree, min_degree+1, tau);
	
	if (fabs(media1-average_k)>fabs(media2-average_k))
		min_degree++;
	
	
	
	
	
	
		
	// range for the community sizes
	if (!fixed_range) {
	
		nmax=max_degree;
		nmin=max(int(min_degree), 3);
		cout<<"-----------------------------------------------------------"<<endl;
		cout<<"community size range automatically set equal to ["<<nmin<<" , "<<nmax<<"]"<<endl;
		

	}
	
	
	//----------------------------------------------------------------------------------------------------
	
	
	deque <int> degree_seq ;		//  degree sequence of the nodes
	deque <double> cumulative;
	powerlaw(max_degree, min_degree, tau, cumulative);
	
	for (int i=0; i<num_nodes; i++) {
		
		int nn=lower_bound(cumulative.begin(), cumulative.end(), ran4())-cumulative.begin()+min_degree;
		degree_seq.push_back(nn);
	
	}
	
	
	
	sort(degree_seq.begin(), degree_seq.end());
		
	if(deque_int_sum(degree_seq) % 2!=0)
		degree_seq[max_element(degree_seq.begin(), degree_seq.end()) - degree_seq.begin()]--;
	
	
	deque<deque<int> >  member_matrix;
	deque<int> num_seq;
	deque<int> internal_degree_seq;
	
	// ********************************			internal_degree and membership			***************************************************

	
	

	if(internal_degree_and_membership(mixing_parameter, overlapping_nodes, overlap_membership, num_nodes, member_matrix, excess, defect, degree_seq, num_seq, internal_degree_seq, fixed_range, nmin, nmax, tau2)==-1)
		return -1;
	
	
	
		
	
	
	
	deque<set<int> > E;					// E is the adjacency matrix written in form of list of edges
	deque<deque<int> > member_list;		// row i cointains the memberships of node i
	deque<deque<int> > link_list;		// row i cointains degree of the node i respect to member_list[i][j]; there is one more number that is the external degree

	
	
	cout<<"building communities... "<<endl;
	if(build_subgraphs(E, member_matrix, member_list, link_list, internal_degree_seq, degree_seq, excess, defect)==-1)
		return -1;	
	




	cout<<"connecting communities... "<<endl;
	connect_all_the_parts(E, member_list, link_list);
	



	if(erase_links(E, member_list, excess, defect, mixing_parameter)==-1)
		return -1;
	

	if(ca!=unlikely) {
		cout<<"trying to approach an average clustering coefficient ... "<<ca<<endl;
		cclu(E, member_list, member_matrix, ca);
	}

	
	deque<map <int, double > > neigh_weigh;

	cout<<"inserting weights..."<<endl;
	weights(E, member_list, beta, mixing_parameter2, neigh_weigh);
	
	
	

	cout<<"recording network..."<<endl;	
	print_network(E, member_list, member_matrix, num_seq, neigh_weigh, beta, mixing_parameter2, mixing_parameter);

	
	
	
		
	return 0;
	
}





void erase_file_if_exists(string s) {

	char b[100];
	cast_string_to_char(s, b);
	
	
	ifstream in1(b);
	
	if(in1.is_open()) {
		
		char rmb[120];
		sprintf(rmb, "rm %s", b);

		int erase= system(rmb);
	}


}



int main(int argc, char * argv[]) {
	
		
	srand_file();
	Parameters p;
	if(set_parameters(argc, argv, p)==false) {
		
		if (argc>1)
			cerr<<"Please, look at ReadMe.txt..."<<endl;
		
		return -1;
	
	}
	
	
	erase_file_if_exists("network.dat");
	erase_file_if_exists("community.dat");
	erase_file_if_exists("statistics.dat");

	
	benchmark(p.excess, p.defect, p.num_nodes, p.average_k, p.max_degree, p.tau, p.tau2, p.mixing_parameter,  p.mixing_parameter2,  p.beta, p.overlapping_nodes, p.overlap_membership, p.nmin, p.nmax, p.fixed_range, p.clustering_coeff);	
		
	return 0;
	
}


