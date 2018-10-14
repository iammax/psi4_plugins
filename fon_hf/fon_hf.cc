///SDF - FDS commutator

#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>
#include <typeinfo>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <cmath>
INIT_PLUGIN

using namespace boost;

namespace psi{ namespace fon_hf {

double largest_element(SharedMatrix mat, int nso) //Finds the largest single element in a matrix
{
	double max = mat->get(0,0);
	for (int i = 0; i < nso; i++)
	{
		for (int j = 0; j < nso; j++)
		{
			if (mat->get(i,j) > max)
			{
				max = mat->get(i, j);
			}
		}
	}
	return max;
}


//for ethylene test case: 7 orbitals, want bottom 4 closed, 5 and 6 open, 7 closed. at ground state 1-5 are filled	
		
double populate_orbitals(SharedVector ni, SharedVector ei, int nso, float T, int num_electrons, double ef, int active, int closed)
{
	double K = 1; //boltsmann = 1 in AU?
	double kt = K*T;
	double denom;
	double this_energy;
	double this_pop;
	double total_pop = 0;
	
	for (int i = 0; i < closed; i++)
	{
//		fprintf(outfile, "closing %d\n", i);
		ni->set(i, 2);
		total_pop += 2;
	}
	for (int i = closed; i < closed + active; i++)
	{
//		fprintf(outfile, "active %d\n", i);
		this_energy = ei->get(i);
		denom = 1 + exp((this_energy - ef)/kt);
		this_pop = 2/denom;
		ni->set(i, this_pop);
		total_pop += this_pop;
	}
	for (int i = closed + active; i < nso; i++)
	{
//		fprintf(outfile, "top %d\n", i);
		ni->set(i, 0);	
	}

	return total_pop;
}


double fermi_energy(SharedVector ni, SharedVector ei, int nso, float T, int num_electrons, float tolerance, int active, int closed)
{
	int target = num_electrons;
	bool correct = 0;
	int tries = 0;
	double diff;
	double sum;
	//start by finding sum at ef = 0
	double prev = populate_orbitals(ni, ei, nso, T, num_electrons, 0, active, closed);
//	fprintf(outfile, "Try # %d. ef = %f sum = %f diff = %f\n", tries, 0.0, prev, fabs(target-sum));
	double bigstep = 125;
	tries += 1;
//	fprintf(outfile, "Starting fermi energy finder.\n");
	double ef = 250; //next, find at 250 and after that, it'll decide whether to bisect or expand
	double prev_en = 0;
	double displacement;
	while (correct == 0)
	{
		sum = populate_orbitals(ni, ei, nso, T, num_electrons, ef, active, closed);
		diff = (sum - target);
//		fprintf(outfile, "Try # %d. ef = %f prev_ef = %f sum = %f diff = %f\n", tries, ef, prev_en, sum, diff);
		if (sum > target)
		{
			if (prev > target) 
			{
				prev_en = ef;
				ef -= bigstep;
				bigstep *= .9;
			}
			else if (prev < target) 
			{
				displacement = fabs(ef - prev_en)/2;
				prev_en = ef;
				ef -= displacement;

			}
		}
		else if (sum < target)
		{
			if (prev < target) 
			{
				prev_en = ef;	
				ef += bigstep;
				bigstep *= .9;
			}
			else if (prev > target) 
			{
				displacement = fabs(ef - prev_en)/2;
				prev_en = ef;		
				ef += displacement;
			}
		}
		if (fabs(diff) < tolerance)
		{
//			fprintf(outfile, "Found close enough fermi energy.\n");
			correct = 1;
		}
		tries += 1;
		prev = sum;
		if (tries > 5000)
		{
			fprintf(outfile, "Warning! Could not properly find fermi energy to populate orbitals!\n");
			exit (EXIT_FAILURE);
		}
	}
	return ef;

}

void matrix_multiply_by_factor(SharedMatrix mat, int dim, double factor)
{
	for (int i = 0; i<dim; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			double mij = mat->get(i, j);
			mat->set(i, j, mij*factor);
		}
	}
}

double frobenius(SharedMatrix A, SharedMatrix B, int nso)
{
	double total = 0;
	for (int i = 0; i < nso; i++)
	{
		for (int j = 0; j < nso; j++)
		{
			total += (A->get(i, j)*B->get(i, j));
		}
	}
	return total;		
}

void make_density_matrix(SharedMatrix Dij, SharedMatrix Cij, int num_electrons, int num_orbitals, SharedVector Ni)
{
	double nm;
	for (int i = 0; i < num_orbitals; i++)
	{
		for (int j = 0; j < num_orbitals; j++)
		{
		
			double total = 0;
			for (int m = 0; m < num_orbitals; m++)
			{
				nm = Ni->get(m)/2;
//				nm = 1;
				double first = Cij->get(i, m);
				double second = Cij->get(j, m);
				total += (nm*first*second);
			}
			Dij->set(i, j, total);
		}
	}
}
	
double electronic_energy(SharedMatrix Dij, SharedMatrix Fij, SharedMatrix Hij, int num_orbitals)
{
	double energy = 0;
	for (int i = 0; i < num_orbitals; i++)
	{
		for (int j = 0; j < num_orbitals; j++)
		{
			double densitypart = Dij->get(i, j);	
			double hpart = Hij->get(i, j);
			double fpart = Fij->get(i, j);
			double term = densitypart*(hpart + fpart);
//			fprintf(outfile, "i, j, pij, hij, fij, term: %d %d %f %f %f %f \n", i, j, densitypart, hpart, fpart, term);

			energy += term;
		}
	}
	return energy;
}

double compute_energy(SharedMatrix F_fock,SharedMatrix H_core,SharedMatrix D_density,int nso){
	SharedMatrix HplusF(new Matrix("HplusF",nso,nso));
	HplusF->copy(H_core);
	HplusF->add(F_fock);
	return HplusF->vector_dot(D_density);
}



extern "C" 
int read_options(std::string name, Options& options)
{
    if (name == "FON_HF"|| options.read_globals()) {
        /*- The amount of information printed to the output file -*/
        options.add_int("PRINT", 1);
    }

    return true;
}

extern "C" 
PsiReturnType fon_hf(Options& options)
{
    int print = options.get_int("PRINT");

    /* Your code goes here */
fprintf(outfile, "\n\n--------------Beginning of customized output------------\n\n");
fprintf(outfile, "Pre-iteration info:\n");

//Basis set and molecular info
shared_ptr<Molecule> mol = Process::environment.molecule();
shared_ptr<Wavefunction> wfn = Process::environment.reference_wavefunction(); 
int nso = wfn->nso(); 
fprintf(outfile, "Number of orbital: %d\n", nso);
int nmo = wfn->nmo();
int numatoms =mol->natom();
fprintf(outfile, "Number of atoms: %i\n", numatoms);
double eSCF = wfn->reference_energy();
SharedVector _Evals = wfn->epsilon_a();
double *Evals = _Evals->pointer();
shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
shared_ptr<BasisSet> aoBasis = BasisSet::construct(parser, mol, "BASIS");
//aoBasis->print_detail();
//fprintf(outfile, "==>> START BASISSET <<==\n" );
//for (unsigned int n=0;n<aoBasis->nshell();n++)
//{
//fprintf(outfile, "SHELL %d is: \n",n );
//aoBasis->shell(n).print(outfile);
//}  
//fprintf(outfile, "==>> FINISH BASISSET <<==\n\n" );
fprintf(outfile, "Reference energy: %f\n",eSCF);
int num_electrons = 0; 
for (int i = 0; i < numatoms; i++) { //calculates total number of electrons
 num_electrons += mol->Z(i);
}
fprintf (outfile, "Number of electrons: %d\n", num_electrons);
double enuc = mol->nuclear_repulsion_energy();
fprintf (outfile, "Nuclear repulsion energy: %f\n", enuc);

//configurable stuff is here

int active = 7; //Number of active orbitals in active space
//int closed = (num_electrons)/2-2; //Number of closed orbitals in active space
int closed = 0;
fprintf(outfile, "Automatically setting 4/3 active space\n");
float T = 0.35; //Temperature to use for FON
bool print1e = 0; //to print one electron matrices pre-iterating to outfile
bool printeri = 0; //1 to put ERIs in outfile, 0 to not do so
bool print1e_in_iter = 0; //to print matrices in each iteration
float tolerance = 1e-7; //tolerance for fermi energy to force sum of electron occupations to = total num
int max_iterations = 50000;
double convthre = 1e-8;
double commutator_convthre = 1e-8;


//construct + print one-electron matrices
shared_ptr<IntegralFactory> integral = shared_ptr<IntegralFactory>(new IntegralFactory(aoBasis, aoBasis, aoBasis, aoBasis));
shared_ptr<TwoBodyAOInt> eri = shared_ptr<TwoBodyAOInt>(integral->eri(0));
shared_ptr<MatrixFactory> factory = shared_ptr<MatrixFactory>(new MatrixFactory);
factory->init_with(1, &nso, &nso);
shared_ptr<OneBodyAOInt> Soverlap(integral->ao_overlap());
const double *buffer;
shared_ptr<Matrix> Sij(factory->create_matrix("Overlap Sij"));
Soverlap->compute(Sij); //this line seems to kill it?
shared_ptr<OneBodyAOInt> tOBI(integral->ao_kinetic());
shared_ptr<OneBodyAOInt> vOBI(integral->ao_potential());
SharedMatrix Tij(factory->create_matrix("Kinetic Tij"));
SharedMatrix Vij(factory->create_matrix("Potential Vij"));
SharedMatrix Hij(factory->create_matrix("Hamiltonian Hij = Tij + Vij"));
tOBI->compute(Tij);
vOBI->compute(Vij);
Hij->copy(Tij);
Hij->add(Vij);
if (print1e){
	fprintf(outfile, "Now printing matrices...\n\n");
	Sij->print();
	Tij->print();
	Vij->print();
	Hij->print();
}


//if 1, outputs Sjj, Vij, Tij matrices to binary files which are readable by numpy with numpy.fromfile
if (0){
FILE*  sout;
sout = fopen("Sij.dat", "w");
fwrite(Sij->pointer()[0], sizeof(double), nso * nso, sout) ;
fclose(sout);

FILE*  vout;
vout = fopen("Vij.dat", "w");
fwrite(Vij->pointer()[0], sizeof(double), nso * nso, vout) ;
fclose(vout);

FILE*  tout;
tout = fopen("Tij.dat", "w");
fwrite(Tij->pointer()[0], sizeof(double), nso * nso, tout) ;
fclose(tout);
}



//construct + print two-electron tensor

if (printeri){ fprintf(outfile,"Two-Electron Integrals: \n"); }



buffer = eri->buffer();
double eri_tensor [nso][nso][nso][nso];


for(int MU=0; MU<aoBasis->nshell(); MU++) 
{
 int nummu = aoBasis->shell(MU).nfunction();
 for(int NU=0; NU<aoBasis->nshell(); NU++) 
 {
  int numnu = aoBasis->shell(NU).nfunction();
  for(int RHO=0; RHO<aoBasis->nshell(); RHO++) 
  {
   int numrho = aoBasis->shell(RHO).nfunction();
   for(int SIG=0; SIG<aoBasis->nshell(); SIG++) 
   {  
    int numsig = aoBasis->shell(SIG).nfunction();
    eri->compute_shell(MU,NU,RHO,SIG);
    for(int mu=0, munu=0, index=0; mu<nummu; mu++) 
    {
     int omu = aoBasis->shell(MU).function_index() + mu;
     for(int nu=0; nu<numnu; nu++, munu++) 
     {
      int onu = aoBasis->shell(NU).function_index() + nu;
      for(int rho=0; rho<numrho; rho++) 
      {
       int orho = aoBasis->shell(RHO).function_index() + rho;
       for(int sig=0; sig<numsig; sig++, index++) 
       {
        int osig = aoBasis->shell(SIG).function_index() + sig; 
         if (printeri) {fprintf(outfile,"( %2d %2d | %2d %2d ) = %20.14lf\n", omu,onu,orho,osig,buffer[index]);}
      	 eri_tensor[omu][onu][orho][osig] = buffer[index];
       }
      }
     }
    }
   } 
  } 
 }
}


//if 1, print out the ERI tensor, as a test
if (0) {
	for (int i=0; i < nso; i++) {
		for (int j = 0; j < nso; j++) {
			for (int k = 0; k < nso; k++) {
				for (int l = 0; l < nso; l++) {
					std::cout << "i, j, k, l, ERI: " << i << " " << j << " " << k << " " << l << " " << eri_tensor[i][j][k][l] << std::endl;
				}
			}
		}
	}
}


//Prepration of 1 and 2 electron integrals is finished; now, begin preparing things for the SCF procedure

//Symmetric orthoganalization matrix
SharedMatrix Xij(new Matrix("Transformation matrix Xij", nso, nso));
Xij->copy(Sij);
Xij->power(-0.5); //Xij = S^-1/2

//Fock matrix
SharedMatrix F0ij(new Matrix("F0ij used to get C", nso, nso));
SharedMatrix HX(new Matrix("H times X", nso, nso));
HX->gemm(false, false, 1.0, Hij, Xij, 0.0); 
F0ij->gemm(true, false, 1.0, Xij, HX, 0.0);

//coeff matrix and initial energies
SharedMatrix Cij_Transformed(new Matrix("Temporary orthonormalized coefficient matrix", nso, nso));
SharedVector Ei(new Vector("Orbital energies", nso));
SharedVector Ni(new Vector("Orbital populations", nso));
SharedMatrix Cij(new Matrix("Coefficient matrix Cij, in original basis", nso, nso));
F0ij->diagonalize(Cij_Transformed, Ei);
Cij->gemm(false, false, 1.0, Xij, Cij_Transformed, 0);


//initial guess for density matrix + calculate initial energy. Initial fock is just the hamiltonian, becuase Gij = 0
SharedMatrix Gij(new Matrix("Gij from ERIs", nso, nso));
SharedMatrix Dij(new Matrix("Density matrix", nso, nso));
SharedMatrix Fij(new Matrix("Fock matrix", nso, nso));
Fij->copy(Hij);
Fij->add(Gij);
if (print1e){
Xij->print();
Fij->print();
F0ij->print();
Cij_Transformed->print();
Cij->print();
Dij->print();
}


//FON stuff

//fprintf(outfile, "Initial orbital energies and populations: \n");
double ef = fermi_energy(Ni, Ei, nso, T, num_electrons, tolerance, active, closed);
double total_pop = populate_orbitals(Ni,Ei, nso, T, num_electrons, ef, active, closed);
Ei->print();
Ni->print();
fprintf(outfile, "Total population: %f\n", total_pop);
//fprintf(outfile, "Fermi energy which gives the above total pop: %f\n", ef);


//end FON for pre-iteration

make_density_matrix(Dij, Cij, num_electrons, nso, Ni);
double eelec = electronic_energy(Dij, Fij, Hij, nso);
double prev_en = eelec;
bool converged = 0;
int iter_count = 0;
double etotal = enuc + eelec;
double delta = etotal;


//iterate
fprintf(outfile, "Beginning iterating. Convergence thresholds: Delta: %e S-D-F Commutator: %e \n", convthre, commutator_convthre);
fprintf(outfile, "Using DIIS for convergence acceleration\n");
fprintf(outfile, "Temperature used for FON: %f atomic units\n", T);
fprintf(outfile, "Tolerance for forcing the sum of fractional occupations to sum to the correct electron #: %e\n\n", tolerance);
fprintf(outfile, "Iter eelec                etotal              delta                  max_commutator\n");

//list for DIIS error matrices

int dn = nso; //dn = diis num = number of matrices to keep for DIIS
SharedMatrix error_mat_list [dn]; //error matrix list for DIIS
SharedMatrix fock_mat_list [dn]; //fock matrix list for DIIS
double largest_in_error;

while (not converged && iter_count < max_iterations)
{	
	//
//	fprintf(outfile, "gij for iteration %d\n", iter_count);
//	Gij->print();
	//	

	
	prev_en = etotal;
	SharedMatrix Fij(new Matrix("Fock Matrix", nso, nso));
	for (int mu = 0; mu < nso; mu++)
	{
		for (int nu = 0; nu < nso; nu++)
		{
		double total = 0;
			for (int lam = 0; lam < nso; lam++)
			{		
				for (int sig = 0; sig < nso; sig++)
				{
					double densitypart = Dij->get(lam, sig);	
					double firsteri = eri_tensor[mu][nu][lam][sig];
					double seconderi = eri_tensor[mu][lam][nu][sig];
					total += densitypart*(firsteri + firsteri - seconderi);
			
				}
			}
		Gij->set(mu, nu, total);
		}
	}

	
	Fij->copy(Gij);
	Fij->add(Hij);	
//	fprintf(outfile, "Normal fock matrix, before DIIS part, for iter %d\n:", iter_count);
//	Fij->print();
//	fprintf(outfile, "Density matrix used to make the above fock: \n");
//	Dij->print();
	etotal = eelec + enuc;




	


//DIIS attempts
	SharedMatrix F_D(new Matrix("Fock times Density", nso, nso));
	SharedMatrix error_left(new Matrix("Error matrix left side", nso, nso));
	SharedMatrix S_D(new Matrix("Overlap times Density", nso, nso));
	SharedMatrix error_right(new Matrix("Error matrix right side", nso, nso));
	SharedMatrix error_matrix(new Matrix("Error matrix", nso, nso));
	F_D->gemm(false, false, 1.0, Fij, Dij, 0);
	error_left->gemm(false, false, 1.0, F_D, Sij, 0);
	S_D->gemm(false, false, 1.0, Sij, Dij, 0);
	error_right->gemm(false, false, 1.0, S_D, Fij, 0);
	error_matrix->copy(error_left);
	error_matrix->subtract(error_right);
	SharedMatrix tempfij(new Matrix("temp Fock", nso, nso));
	SharedMatrix tempfij2(new Matrix("temp Fock 2", nso, nso));
	largest_in_error = largest_element(error_matrix, nso);
//	fprintf(outfile, "For iteration %d, the largest value in the error: %e\n", iter_count, largest_in_error);


//Print results of the previous iteration

	if(iter_count == 0)
	{
		delta = etotal;
	}
	else
	{
		delta = etotal - prev_en;
	}
	fprintf(outfile, "%d %20.14f %20.14f %20.14f %20.14f\n", iter_count, eelec, etotal, delta, largest_in_error);

	
	if (fabs(largest_in_error) < commutator_convthre)
	{
		fprintf(outfile, "Converged via commutator\n");
		converged = 1;
	}


//DIIS: Prepare error matrix list and fock m	atrix list
//	fprintf(outfile, "Error matrix for iteration %d\n:", iter_count);
//	error_matrix->print();
	if (iter_count < dn)
	{
		error_mat_list[iter_count] = error_matrix->clone();
		fock_mat_list[iter_count] = Fij->clone();
	}
	else
	{
		for (int i = 0; i < dn-1; i++)
		{
//			fprintf(outfile, "%d\n", i);
			error_mat_list[i] = error_mat_list[i+1]->clone();
			fock_mat_list[i] = fock_mat_list[i+1]->clone();
	
	
		}
		error_mat_list[dn-1] = error_matrix->clone();
		fock_mat_list[dn-1] = Fij->clone();
	}		
	
	


	int bmatsize;
	if (iter_count < dn)
	{
		bmatsize = iter_count + 1;
	}
	else
	{
		bmatsize = dn;
	}	
	double target [bmatsize+1];

//	for (int q = 0; q < bmatsize; q++)	
//	{
//		fprintf(outfile, "fock matrix %d for iteration %d: \n", q, iter_count);
//		fock_mat_list[q]->print();
//	}
//	for (int q = 0; q < bmatsize; q++)	
//	{
//		fprintf(outfile, "error matrix %d for iteration %d: \n", q, iter_count);
//		error_mat_list[q]->print();
//	}


//	fprintf(outfile, "check 1\n");
	if (bmatsize>0 && iter_count > 0)
	{

		SharedMatrix Bmat(new Matrix("DIIS B Matrix", bmatsize+1, bmatsize+1));
		double bijresult;
		for (int i = 0; i <= bmatsize; i++)
		{
			target[i] = 0;
			if ( i==bmatsize)
			{
				target[i] = -1;
			}
			for (int j = i; j <= bmatsize; j++)
			{
//				fprintf(outfile, "check for i, j = %d %d\n", i, j);
				if (i==bmatsize ^ j == bmatsize)
				{
//					fprintf(outfile, "xor\n");
					Bmat->set(i, j, -1);	
					Bmat->set(j, i, -1);
					bijresult = -1;
				}
				else if (i==bmatsize & j == bmatsize)
				{
//					fprintf(outfile, "both\n");
					Bmat->set(i, j, 0);
					bijresult = 1;
				}
				else
				{
//					fprintf(outfile, "neither\n");
					double Bij = frobenius(error_mat_list[i], error_mat_list[j], nso);
					Bmat->set(i, j, Bij);
					Bmat->set(j, i, Bij);
					bijresult = Bij;
				}
//				fprintf(outfile, "i, j, Bij: %d %d %f\n", i, j, bijresult);
			}
		}
	
		int output_array [bmatsize];
//		Bmat->print();
//		fprintf(outfile, "Target vector: \n");
//		for (int i = 0; i < bmatsize+1; i++)
//		{
//			fprintf(outfile, "%f\n", target[i]);
//		}


		int a = C_DGESV(bmatsize+1, 1, Bmat->pointer()[0], bmatsize+1, output_array, target, bmatsize+1);
	//	after the above line, "target" should be a list of the C's with the last one being the lagrange multiplier
		if (a != 0)
		{
			fprintf(outfile, "for iteration %d, system of linear equations for DIIS seems to have failed", iter_count);	
		}
	
		double c_total = 0;
	
		for (int i = 0; i < bmatsize; i++)
		{
//			fprintf(outfile, "c sub %d = %f\n", i, target[i]);
			c_total += target[i];
		}
//		fprintf(outfile, "sum of C's = %f\n", c_total);
//		fprintf(outfile, "lagrange multiplier (end of C list): %f\n", target[bmatsize]);
//		fprintf(outfile, "for iteration %d, system of linear equations for DIIS return value: %d\n", iter_count, a);
//		
//		
//		fprintf(outfile, "Solution vector: \n");
//		for (int i = 0; i < bmatsize +1 ; i++)
//		{
//			fprintf(outfile, "%f\n", target[i]);
//		}
		
	}
	
//////
//	fprintf(outfile, "now printing fock matrix list...\n");	
//	for (int i = 0; i < bmatsize; i++)
//	{
//		fprintf(outfile, "%d\n", i);
//		fock_mat_list[i]->print();
//	}

	SharedMatrix fock_i(new Matrix("temporary fock matrix for DIIS", nso, nso));
	double c_i;
	SharedMatrix DIISFij(new Matrix("DIIS Fock", nso, nso));
	if (iter_count > 0)
	{
//		fprintf(outfile, "DIIS fock. Should be blank here");	
//		DIISFij->print();
		if (iter_count < dn)
		{
			for (int i = 0; i < iter_count+1; i++)
			{
//				fprintf(outfile, "DIIS fock so far: \n");
//				DIISFij->print();
				c_i = target[i];
				fock_i = fock_mat_list[i]->clone();
//				fprintf(outfile, "fock matrix # %d for iteration # %d\n", i, iter_count);
//				fock_i->print();
				matrix_multiply_by_factor(fock_i, nso, c_i);
//				fprintf(outfile, "Factor to multiply by: %f\n", c_i);
				DIISFij->add(fock_i);
//				fprintf(outfile, "fock matrix after multiplying:\n");
//				fock_i->print();
//				fprintf(outfile, "DIIS fock after adding: \n");
//				DIISFij->print();
			}
		}
		else
		{
			for (int i = 0; i < dn; i++)
			{
//				fprintf(outfile, "DIIS fock so far: \n");
//				DIISFij->print();
				c_i = target[i];
				fock_i = fock_mat_list[i] ->clone();
//				fprintf(outfile, "fock matrix # %d for iteration # %d\n", i, iter_count);
//				fock_i->print();
				matrix_multiply_by_factor(fock_i, nso, c_i);
//				fprintf(outfile, "Factor to multiply by: %f\n", c_i);
				DIISFij->add(fock_i);
//				fprintf(outfile, "fock matrix after multiplying:\n");
//				fock_i->print();
//				fprintf(outfile, "DIIS fock after adding: \n");
//				DIISFij->print();
			}
			error_mat_list[dn-1] = error_matrix;
			fock_mat_list[dn-1] = Fij;
		}	
//	fprintf(outfile, "Regular fock: \n");
//	Fij->print();	
//	fprintf(outfile, "DIIS Fock: \n");
//	DIISFij->print();
	Fij = (DIISFij)->clone();
//	fprintf(outfile, "DIIS-adjusted fock for iteration %d\n:", iter_count);
//	DIISFij->print();
	}
	
///////	


	
//end DIIS

	SharedMatrix F_times_X(new Matrix("F times X", nso, nso));
	F_times_X->gemm(false, false, 1.0, Fij, Xij, 0.0);
	F0ij->gemm(true, false, 1.0, Xij, F_times_X, 0.0);
	F0ij->diagonalize(Cij_Transformed, Ei);
	Cij->gemm(false, false, 1.0, Xij, Cij_Transformed, 0);
//FON stuff for iterations
//	fprintf(outfile, "Orbital energies and populations for iteration %d: \n", iter_count);
	double ef = fermi_energy(Ni, Ei, nso, T, num_electrons, tolerance, active, closed);
	double total_pop = populate_orbitals(Ni,Ei, nso, T, num_electrons, ef, active, closed);
	Ei->print();
	Ni->print();
//	fprintf(outfile, "Total population: %f\n", total_pop);
//	fprintf(outfile, "Total population 1: (1e-7)%f\n", total_pop);
//	fprintf(outfile, "Fermi energy which gives the above total pop: %f\n", ef);

//end FON stuff
	make_density_matrix(Dij, Cij, num_electrons, nso, Ni);
	eelec = electronic_energy(Dij, Fij, Hij, nso);

	if (print1e_in_iter)
	{
		fprintf(outfile, "Iter num: %d\n", iter_count);
		Gij->print();
		Fij->print();
		F0ij->print();
		Cij->print();
		Dij->print();
	}
	if (iter_count>1) //regular iterative convergence
	{
		if (fabs(delta) < convthre)
		{
			converged = 1;
			fprintf(outfile, "Iterations converged delta_energy in iterations\n");

		}
	}	
	iter_count++;
	if (iter_count == max_iterations)
	{
		fprintf(outfile, "Max iterations reached\n");
	}
}

fprintf(outfile, "Final orbital energies and populations: \n");
Ei->print();
Ni->print();

//fprintf(outfile, "Populations (1 = doubly occupied, etc): \n");
//for (int i = 0; i < nso; i++)
//{
//	fprintf(outfile, "      %d: %f\n", i+1, Ni->get(i)/2);
//}

fprintf(outfile, "Total population: %f\n", total_pop);


fprintf(outfile, "Reference energy from built in SCF: %f\n",eSCF);
if (converged){
double diff = etotal - eSCF;
fprintf(outfile, "Difference between your SCF and reference SCF: %f\n", diff);
fprintf(outfile, "Finished with run!\n");
}


    return Success;
}

}} // End namespaces

