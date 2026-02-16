/*
 * GENERAL REMARKS                                                                          
 *                                                                                          
 *  This code is freely available under the following conditions:                           
 *                                                                                          
 *  1) The code is to be used only for non-commercial purposes.                             
 *  2) No changes and modifications to the code without prior permission of the developer.  
 *  3) No forwarding the code to a third party without prior permission of the developer.   
 *                                                                                          
 *                       
 *  This file contains implementation of the COCR solver with folded preconditioner and minimal residual smoothing               
 *  The matix stored in the sparse format  
 *                                                                                      
 *  Earlier uploaded: https://github.com/kiselev2013/MTCalc_with_DFP_COCR/blob/master/Code/COCR_FP/pcocr.cpp                                                       
 *  Novosibirsk State Technical University,                                                 
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                       
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                                                                       
 *                                                                
*/                                                                                          

#include "stdafx.h"
#include "pcocr.h"
#include "pcocr_rci.h"
#include "block_2x2_solver.h"
#include "MRS.h"
#include "FoldedPreconditioner.h"
extern ofstream logfile;
#include <chrono>

//------------------------------------------------------------------------
PCOCR::PCOCR()
{
}
//------------------------------------------------------------------------
PCOCR::~PCOCR()
{
}
//------------------------------------------------------------------------
int PCOCR::PCOCR_2x2_Folded(int n, int *ig, int *jg, int *idi, int *ijg, double *di, double *gg, double *pr,
					 double *x, double eps, int maxiter, double *y_omp,
					 int kpar, int n_edges_c, int n_nodes_c,
					 double (*xyz)[3], int *nvkat, int (*nver)[14], double *sigma3d, int (*edges)[2],
					 int *ig_t, int *jg_t, double *gg_t, int *is_node_bound)
{
	int nb = n/2;
	int req;
	double *a=NULL, *b=NULL;
	Block_2x2_solver s;
	FoldedPreconditioner p;

	cout << "PCOCR_2x2_Folded...\n";
	logfile << "PCOCR_2x2_Folded...\n";

	__time64_t time_total, time_beg, time_end; 
	__time64_t timestruct;

	_time64(&timestruct);
	time_beg = timestruct;

	auto tStart = std::chrono::high_resolution_clock::now();


#if defined(__PROGRESS_INDICATOR_INCLUDED)
	CString sbuf;
	GUI::ProgressIndicator * pin=NULL;
	pin=new GUI::ProgressIndicator("Solving SLAE using PCOCR_Fold_MRS...", int(-log10(eps)*1000));
#endif	

	PCOCR_RCI pcocr(n, maxiter, eps, x, pr, &a, &b);

	p.BuildGMatrix(kpar, n_edges_c, n_nodes_c, nvkat, nver, sigma3d, edges, ig_t, jg_t, gg_t, is_node_bound);
	p.BuildDiagComplex(ig, jg, idi, ijg, di, gg);

	MRS mrs(n);

	do 
	{
		req = pcocr.Run();
	
		switch(req)
		{
		case REQ_MULT_MV:
			s.Mult_MV_block_2x2((int)nb, (int*)ig, (int*)jg, (int*)idi, (int*)ijg, di, gg, a, b, y_omp);
			break;

		case REQ_PRECOND:
			p.ApplyPreconditionerComplex(a, b);
			break;	

		case REQ_STOP_TEST:
			mrs.Mrs(x, a);
			pcocr.DoStopTest(mrs.GetResidualNorm(), &req);
			if (pcocr.iter % 10 == 0) {
				pcocr.PrintIterResidual();
			}
#if defined(__PROGRESS_INDICATOR_INCLUDED)
			sbuf.Format("PCOCR_Fold_MRS: eps=%5.2e, eps current=%e, iter=%d", eps, pcocr.residualRel, pcocr.iter);
			pin->SetText(sbuf);
			pin->SetPos(int(-log10(pcocr.residualRel)*1000));
#endif 
			break;

		case REQ_X0_TEST:
			mrs.Mrs(x, a);
			pcocr.DoX0Test(a, &req);
			break;
		}
	}
	while(req > 0);

#if defined(__PROGRESS_INDICATOR_INCLUDED)
	delete pin;
#endif

	_time64(&timestruct);
	time_end = timestruct;
	time_total = time_end - time_beg;


	auto tEnd = std::chrono::high_resolution_clock::now();
	Base_solver bs;
	bs.WriteKitChrono("kit", pcocr.residualRel, eps, pcocr.iter, std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1>>>(tEnd - tStart).count());

	return req;
}
