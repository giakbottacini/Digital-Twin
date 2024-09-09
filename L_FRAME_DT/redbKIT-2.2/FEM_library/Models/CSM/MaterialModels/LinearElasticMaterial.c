/*   This file is part of redbKIT.
 *   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
 *   Author: Federico Negri <federico.negri@epfl.ch>
 */

#include "LinearElasticMaterial.h"

/*************************************************************************/
void LinearElasticMaterial_forces(mxArray* plhs[], const mxArray* prhs[])
{

    double* dim_ptr = mxGetPr(prhs[0]);
    int dim     = (int)(dim_ptr[0]);
    int noe     = mxGetN(prhs[4]);
    double* nln_ptr = mxGetPr(prhs[5]);
    int nln     = (int)(nln_ptr[0]);
    int numRowsElements  = mxGetM(prhs[4]);
    int nln2    = nln*nln;

    plhs[0] = mxCreateDoubleMatrix(nln*noe*dim,1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(nln*noe*dim,1, mxREAL);

    double* myRrows    = mxGetPr(plhs[0]);
    double* myRcoef    = mxGetPr(plhs[1]);

    int k,l;
    int q;
    int NumQuadPoints     = mxGetN(prhs[6]);
    int NumNodes          = (int)(mxGetM(prhs[3]) / dim);

	/*FILE * fp;
	FILE * fp2;*/

    double* U_h   = mxGetPr(prhs[3]);
    double* w   = mxGetPr(prhs[6]);
    double* invjac = mxGetPr(prhs[7]);
    double* detjac = mxGetPr(prhs[8]);
    double* phi = mxGetPr(prhs[9]);
    double* gradrefphi = mxGetPr(prhs[10]);

    double gradphi[dim][nln][NumQuadPoints];
    double* elements  = mxGetPr(prhs[4]);

    double GradV[dim][dim];
    double GradU[dim][dim];
    double GradUh[dim][dim][NumQuadPoints];

    double Id[dim][dim];
    int d1,d2;
    for (d1 = 0; d1 < dim; d1 = d1 + 1 )
    {
        for (d2 = 0; d2 < dim; d2 = d2 + 1 )
        {
            Id[d1][d2] = 0;
            if (d1==d2)
            {
                Id[d1][d2] = 1;
            }
        }
    }

    double F[dim][dim];
    double EPS[dim][dim];
    double dP[dim][dim];
    double P_Uh[dim][dim];

    double* material_param = mxGetPr(prhs[2]);
    double mu = material_param[0];
    double lambda = material_param[1];

    /* Assembly: loop over the elements */
    int ie;

	/*fp = fopen("C:/Users/Rosafalco/Documents/Polimi/3_Dottorato/Codici/Matlab/Quarteroni_Manzoni_Negri/pippoLE.out", "w");
	fp2 = fopen("C:/Users/Rosafalco/Documents/Polimi/3_Dottorato/Codici/Matlab/Quarteroni_Manzoni_Negri/comparisonLE.out", "w");*/

	#pragma omp parallel for shared(invjac,detjac,elements,myRrows,myRcoef,U_h) private(gradphi,F,EPS,dP,P_Uh,GradV,GradU,GradUh,ie,k,l,q,d1,d2) firstprivate(phi,gradrefphi,w,numRowsElements,nln2,nln,NumNodes,Id,mu,lambda)

    for (ie = 0; ie < noe; ie = ie + 1 )
    {
        for (k = 0; k < nln; k = k + 1 )
        {
            for (q = 0; q < NumQuadPoints; q = q + 1 )
            {
                for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                {
                    gradphi[d1][k][q] = 0;
                    for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                    {
                        gradphi[d1][k][q] = gradphi[d1][k][q] + INVJAC(ie,d1,d2)*GRADREFPHI(k,q,d2);
                    }
                }
            }
        }

        for (q = 0; q < NumQuadPoints; q = q + 1 )
        {
            for (d1 = 0; d1 < dim; d1 = d1 + 1 )
            {
                for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                {
                    GradUh[d1][d2][q] = 0;
                    for (k = 0; k < nln; k = k + 1 )
                    {
                        int e_k;
                        e_k = (int)(elements[ie*numRowsElements + k] + d1*NumNodes - 1);
                        GradUh[d1][d2][q] = GradUh[d1][d2][q] + U_h[e_k] * gradphi[d2][k][q];
                    }
                }
            }
        }

        int iii = 0;
        int ii = 0;
        int a, b, i_c, j_c;

        /* loop over test functions --> a */
        for (a = 0; a < nln; a = a + 1 )
        {
            /* loop over test components --> i_c */
            for (i_c = 0; i_c < dim; i_c = i_c + 1 )
            {
                /* set gradV to zero*/
                for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                {
                    for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                    {
                        GradV[d1][d2] = 0;
                    }
                }

                double rloc = 0;
                for (q = 0; q < NumQuadPoints; q = q + 1 )
                {

                    for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                    {
                        GradV[i_c][d2] = gradphi[d2][a][q];
                    }

                    for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                    {
                        for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                        {
                            F[d1][d2] = Id[d1][d2] + GradUh[d1][d2][q];
                        }
                    }

                    for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                    {
                        for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                        {
                            EPS[d1][d2] = 0.5 * ( F[d1][d2] + F[d2][d1] ) - Id[d1][d2];
                        }
                    }

                    double trace = Trace(dim, EPS);
                    for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                    {
                        for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                        {
                            P_Uh[d1][d2] = 2 * mu * EPS[d1][d2] + lambda * trace * Id[d1][d2];
							/*fprintf(fp, "indici ie=%u q=%u d1=%u d2=%u  \n", ie, q, d1, d2);
							fprintf(fp, "lambda = %f \n", lambda);
							fprintf(fp2, "P_Uh[d1][d2]=%f \n", P_Uh[d1][d2]);*/
                        }
                    }
                    rloc  = rloc + Mdot( dim, GradV, P_Uh) * w[q];
                }

                myRrows[ie*nln*dim+ii] = elements[a+ie*numRowsElements] + i_c * NumNodes;
                myRcoef[ie*nln*dim+ii] = rloc*detjac[ie];
                ii = ii + 1;
            }
        }
    }
}
/*fclose(fp);
fclose(fp2);*/
/*************************************************************************/

void LinearElasticMaterial_jacobian(mxArray* plhs[], const mxArray* prhs[])
{

    double* dim_ptr = mxGetPr(prhs[0]);
    int dim     = (int)(dim_ptr[0]);
    int noe     = mxGetN(prhs[4]);
    double* nln_ptr = mxGetPr(prhs[5]);
    int nln     = (int)(nln_ptr[0]);
    int numRowsElements  = mxGetM(prhs[4]);
    int nln2    = nln*nln;

    plhs[0] = mxCreateDoubleMatrix(nln2*noe*dim*dim,1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(nln2*noe*dim*dim,1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(nln2*noe*dim*dim,1, mxREAL);

    double* myArows    = mxGetPr(plhs[0]);
    double* myAcols    = mxGetPr(plhs[1]);
    double* myAcoef    = mxGetPr(plhs[2]);

    int k,l;
    int q;
    int NumQuadPoints     = mxGetN(prhs[6]);
    int NumNodes          = (int)(mxGetM(prhs[3]) / dim);

	/*FILE * fp3;
	FILE * fp4;*/

    double* U_h   = mxGetPr(prhs[3]);
    double* w   = mxGetPr(prhs[6]);
    double* invjac = mxGetPr(prhs[7]);
    double* detjac = mxGetPr(prhs[8]);
    double* phi = mxGetPr(prhs[9]);
    double* gradrefphi = mxGetPr(prhs[10]);

    double gradphi[dim][nln][NumQuadPoints];
    double* elements  = mxGetPr(prhs[4]);

    double GradV[dim][dim];
    double GradU[dim][dim];
    double GradUh[dim][dim][NumQuadPoints];

    double Id[dim][dim];
    int d1,d2;
    for (d1 = 0; d1 < dim; d1 = d1 + 1 )
    {
        for (d2 = 0; d2 < dim; d2 = d2 + 1 )
        {
            Id[d1][d2] = 0;
            if (d1==d2)
            {
                Id[d1][d2] = 1;
            }
        }
    }

    double F[dim][dim];
    double EPS[dim][dim];
    double dP[dim][dim];
    double P_Uh[dim][dim];

    double* material_param = mxGetPr(prhs[2]);
    double mu = material_param[0];
    double lambda = material_param[1];

    /* Assembly: loop over the elements */
    int ie;

	/*fp3 = fopen("C:/Users/Rosafalco/Documents/Polimi/3_Dottorato/Codici/Matlab/Quarteroni_Manzoni_Negri/pippo_grad.out", "w");
	fp4 = fopen("C:/Users/Rosafalco/Documents/Polimi/3_Dottorato/Codici/Matlab/Quarteroni_Manzoni_Negri/comparisonLEMS_grad.out", "w");*/

#pragma omp parallel for shared(invjac,detjac,elements,myAcols,myArows,myAcoef,U_h) private(gradphi,F,EPS,dP,P_Uh,GradV,GradU,GradUh,ie,k,l,q,d1,d2) firstprivate(phi,gradrefphi,w,numRowsElements,nln2,nln,NumNodes,Id,mu,lambda)

    for (ie = 0; ie < noe; ie = ie + 1 )
    {
        for (k = 0; k < nln; k = k + 1 )
        {
            for (q = 0; q < NumQuadPoints; q = q + 1 )
            {
                for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                {
                    gradphi[d1][k][q] = 0;
                    for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                    {
                        gradphi[d1][k][q] = gradphi[d1][k][q] + INVJAC(ie,d1,d2)*GRADREFPHI(k,q,d2);
                    }
                }
            }
        }

        for (q = 0; q < NumQuadPoints; q = q + 1 )
        {
            for (d1 = 0; d1 < dim; d1 = d1 + 1 )
            {
                for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                {
                    GradUh[d1][d2][q] = 0;
                    for (k = 0; k < nln; k = k + 1 )
                    {
                        int e_k;
                        e_k = (int)(elements[ie*numRowsElements + k] + d1*NumNodes - 1);
                        GradUh[d1][d2][q] = GradUh[d1][d2][q] + U_h[e_k] * gradphi[d2][k][q];
                    }
                }
            }
        }

        int iii = 0;
        int ii = 0;
        int a, b, i_c, j_c;

        /* loop over test functions --> a */
        for (a = 0; a < nln; a = a + 1 )
        {
            /* loop over test components --> i_c */
            for (i_c = 0; i_c < dim; i_c = i_c + 1 )
            {
                /* set gradV to zero*/
                for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                {
                    for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                    {
                        GradV[d1][d2] = 0;
                    }
                }

                /* loop over trial functions --> b */
                for (b = 0; b < nln; b = b + 1 )
                {
                    /* loop over trial components --> j_c */
                    for (j_c = 0; j_c < dim; j_c = j_c + 1 )
                    {
                        /* set gradU to zero*/
                        for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                        {
                            for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                            {
                                GradU[d1][d2] = 0;
                            }
                        }

                        double aloc = 0;
                        for (q = 0; q < NumQuadPoints; q = q + 1 )
                        {
                            for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                            {
                                GradV[i_c][d2] = gradphi[d2][a][q];
                                GradU[j_c][d2] = gradphi[d2][b][q];
                            }


                            for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                            {
                                for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                                {
                                    F[d1][d2] = Id[d1][d2] + GradU[d1][d2];
                                }
                            }

                            for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                            {
                                for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                                {
                                    EPS[d1][d2] = 0.5 * ( F[d1][d2] + F[d2][d1] ) - Id[d1][d2];
                                }
                            }


                            double trace = Trace(dim, EPS);
                            for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                            {
                                for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                                {
                                    dP[d1][d2] = 2 * mu * EPS[d1][d2] + lambda * trace * Id[d1][d2];
									/*fprintf(fp3, "indici ie=%u q=%u d1=%u d2=%u  \n", ie, q, d1, d2);
									fprintf(fp3, "lambda = %f \n", lambda);
									fprintf(fp4, "dP[d1][d2]=%f \n", dP[d1][d2]);*/
                                }
                            }
                            aloc  = aloc + Mdot( dim, GradV, dP) * w[q];
                        }
                        myArows[ie*nln2*dim*dim+iii] = elements[a+ie*numRowsElements] + i_c * NumNodes;
                        myAcols[ie*nln2*dim*dim+iii] = elements[b+ie*numRowsElements] + j_c * NumNodes;
                        myAcoef[ie*nln2*dim*dim+iii] = aloc*detjac[ie];

                        iii = iii + 1;
                    }
                }

            }
        }
    }
}
/*fclose(fp3);
fclose(fp4);*/

/*************************************************************************/
void LinearElasticMaterial_jacobianFast3D(mxArray* plhs[], const mxArray* prhs[])
{

    double* dim_ptr = mxGetPr(prhs[0]);
    int dim     = (int)(dim_ptr[0]);
    int noe     = mxGetN(prhs[4]);
    double* nln_ptr = mxGetPr(prhs[5]);
    int nln     = (int)(nln_ptr[0]);
    int numRowsElements  = mxGetM(prhs[4]);
    int nln2    = nln*nln;

    plhs[0] = mxCreateDoubleMatrix(nln2*noe*dim*dim,1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(nln2*noe*dim*dim,1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(nln2*noe*dim*dim,1, mxREAL);

    double* myArows    = mxGetPr(plhs[0]);
    double* myAcols    = mxGetPr(plhs[1]);
    double* myAcoef    = mxGetPr(plhs[2]);

    int k,l;
    int q;
    int NumQuadPoints     = mxGetN(prhs[6]);
    int NumNodes          = (int)(mxGetM(prhs[3]) / dim);

    double* U_h   = mxGetPr(prhs[3]);
    double* w   = mxGetPr(prhs[6]);
    double* invjac = mxGetPr(prhs[7]);
    double* detjac = mxGetPr(prhs[8]);
    double* phi = mxGetPr(prhs[9]);
    double* gradrefphi = mxGetPr(prhs[10]);

    double* elements  = mxGetPr(prhs[4]);

    int d1,d2;

    double* material_param = mxGetPr(prhs[2]);
    double mu = material_param[0];
    double lambda = material_param[1];

    /* Assembly: loop over the elements */
    int ie;

#pragma omp parallel for shared(invjac,detjac,elements,myAcols,myArows,myAcoef,U_h) private(ie,k,l,q,d1,d2) firstprivate(phi,gradrefphi,w,numRowsElements,nln2,nln,NumNodes,mu,lambda)

    for (ie = 0; ie < noe; ie = ie + 1 )
    {
        double gradphi[NumQuadPoints][dim][nln];

        for (q = 0; q < NumQuadPoints; q = q + 1 )
        {
            /* Compute Gradient of Basis functions*/
            for (k = 0; k < nln; k = k + 1 )
            {
                for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                {
                    gradphi[q][d1][k] = 0;
                    for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                    {
                        gradphi[q][d1][k] = gradphi[q][d1][k] + INVJAC(ie,d1,d2)*GRADREFPHI(k,q,d2);
                    }
                }
            }
        }

        int iii = 0;
        int a, b, i_c, j_c;

        double aloc[nln][dim][nln][dim];
        /* loop over test functions --> a */
        for (a = 0; a < nln; a = a + 1 )
        {
            /* loop over test components --> i_c */
            for (i_c = 0; i_c < 3; i_c = i_c + 1 )
            {
                /* loop over trial functions --> b */
                for (b = 0; b < nln; b = b + 1 )
                {
                    /* loop over trial components --> j_c */
                    for (j_c = 0; j_c < 3; j_c = j_c + 1 )
                    {
                        aloc[a][i_c][b][j_c]  = 0.0;
                    }
                }
            }
        }
        for (q = 0; q < NumQuadPoints; q = q + 1 )
        {
            /* loop over test functions --> a */
            for (a = 0; a < nln; a = a + 1 )
            {
                /* loop over trial functions --> b */
                for (b = 0; b < nln; b = b + 1 )
                {

                    aloc[a][0][b][0] += ( gradphi[q][0][a]*(lambda*gradphi[q][0][b] + 2.0*mu*gradphi[q][0][b]) + mu*gradphi[q][1][b]*gradphi[q][1][a] + mu*gradphi[q][2][b]*gradphi[q][2][a] ) * w[q];

                    aloc[a][0][b][1] += ( lambda*gradphi[q][1][b]*gradphi[q][0][a] + mu*gradphi[q][0][b]*gradphi[q][1][a] ) * w[q];

                    aloc[a][0][b][2] += ( lambda*gradphi[q][2][b]*gradphi[q][0][a] + mu*gradphi[q][0][b]*gradphi[q][2][a] ) * w[q];

                    aloc[a][1][b][0] += ( lambda*gradphi[q][0][b]*gradphi[q][1][a] + mu*gradphi[q][1][b]*gradphi[q][0][a] ) * w[q];

                    aloc[a][1][b][1] += ( gradphi[q][1][a]*(lambda*gradphi[q][1][b] + 2.0*mu*gradphi[q][1][b]) + mu*gradphi[q][0][b]*gradphi[q][0][a] + mu*gradphi[q][2][b]*gradphi[q][2][a] ) * w[q];

                    aloc[a][1][b][2] += ( lambda*gradphi[q][2][b]*gradphi[q][1][a] + mu*gradphi[q][1][b]*gradphi[q][2][a] ) * w[q];

                    aloc[a][2][b][0] += ( lambda*gradphi[q][0][b]*gradphi[q][2][a] + mu*gradphi[q][2][b]*gradphi[q][0][a] ) * w[q];

                    aloc[a][2][b][1] += ( lambda*gradphi[q][1][b]*gradphi[q][2][a] + mu*gradphi[q][2][b]*gradphi[q][1][a] ) * w[q];

                    aloc[a][2][b][2] += ( gradphi[q][2][a]*(lambda*gradphi[q][2][b] + 2.0*mu*gradphi[q][2][b]) + mu*gradphi[q][0][b]*gradphi[q][0][a] + mu*gradphi[q][1][b]*gradphi[q][1][a] ) * w[q];

                }
            }
        }

        for (a = 0; a < nln; a = a + 1 )
        {
            /* loop over test components --> i_c */
            for (i_c = 0; i_c < 3; i_c = i_c + 1 )
            {
                /* loop over trial functions --> b */
                for (b = 0; b < nln; b = b + 1 )
                {
                    /* loop over trial components --> j_c */
                    for (j_c = 0; j_c < 3; j_c = j_c + 1 )
                    {
                        myArows[ie*nln2*9+iii] = elements[a+ie*numRowsElements] + i_c * NumNodes;
                        myAcols[ie*nln2*9+iii] = elements[b+ie*numRowsElements] + j_c * NumNodes;
                        myAcoef[ie*nln2*9+iii] = aloc[a][i_c][b][j_c]*detjac[ie];
                        iii = iii + 1;
                    }
                }
            }
        }
    }

}
/*************************************************************************/
void LinearElasticMaterial_jacobianFast2D(mxArray* plhs[], const mxArray* prhs[])
{

    double* dim_ptr = mxGetPr(prhs[0]);
    int dim     = (int)(dim_ptr[0]);
    int noe     = mxGetN(prhs[4]);
    double* nln_ptr = mxGetPr(prhs[5]);
    int nln     = (int)(nln_ptr[0]);
    int numRowsElements  = mxGetM(prhs[4]);
    int nln2    = nln*nln;

    plhs[0] = mxCreateDoubleMatrix(nln2*noe*dim*dim,1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(nln2*noe*dim*dim,1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(nln2*noe*dim*dim,1, mxREAL);

    double* myArows    = mxGetPr(plhs[0]);
    double* myAcols    = mxGetPr(plhs[1]);
    double* myAcoef    = mxGetPr(plhs[2]);

    int k,l;
    int q;
    int NumQuadPoints     = mxGetN(prhs[6]);
    int NumNodes          = (int)(mxGetM(prhs[3]) / dim);

    double* U_h   = mxGetPr(prhs[3]);
    double* w   = mxGetPr(prhs[6]);
    double* invjac = mxGetPr(prhs[7]);
    double* detjac = mxGetPr(prhs[8]);
    double* phi = mxGetPr(prhs[9]);
    double* gradrefphi = mxGetPr(prhs[10]);

    double* elements  = mxGetPr(prhs[4]);

    int d1,d2;

    double* material_param = mxGetPr(prhs[2]);
    double mu = material_param[0];
    double lambda = material_param[1];

    /* Assembly: loop over the elements */
    int ie;

#pragma omp parallel for shared(invjac,detjac,elements,myAcols,myArows,myAcoef,U_h) private(ie,k,l,q,d1,d2) firstprivate(phi,gradrefphi,w,numRowsElements,nln2,nln,NumNodes,mu,lambda)

    for (ie = 0; ie < noe; ie = ie + 1 )
    {
        double gradphi[NumQuadPoints][dim][nln];

        for (q = 0; q < NumQuadPoints; q = q + 1 )
        {
            /* Compute Gradient of Basis functions*/
            for (k = 0; k < nln; k = k + 1 )
            {
                for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                {
                    gradphi[q][d1][k] = 0;
                    for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                    {
                        gradphi[q][d1][k] = gradphi[q][d1][k] + INVJAC(ie,d1,d2)*GRADREFPHI(k,q,d2);
                    }
                }
            }
        }

        int iii = 0;
        int a, b, i_c, j_c;

        double aloc[nln][dim][nln][dim];
        /* loop over test functions --> a */
        for (a = 0; a < nln; a = a + 1 )
        {
            /* loop over test components --> i_c */
            for (i_c = 0; i_c < 2; i_c = i_c + 1 )
            {
                /* loop over trial functions --> b */
                for (b = 0; b < nln; b = b + 1 )
                {
                    /* loop over trial components --> j_c */
                    for (j_c = 0; j_c < 2; j_c = j_c + 1 )
                    {
                        aloc[a][i_c][b][j_c]  = 0.0;
                    }
                }
            }
        }
        for (q = 0; q < NumQuadPoints; q = q + 1 )
        {
            /* loop over test functions --> a */
            for (a = 0; a < nln; a = a + 1 )
            {
                /* loop over trial functions --> b */
                for (b = 0; b < nln; b = b + 1 )
                {

                    aloc[a][0][b][0] += ( gradphi[q][0][a]*(lambda*gradphi[q][0][b] + 2.0*mu*gradphi[q][0][b]) + mu*gradphi[q][1][b]*gradphi[q][1][a] ) * w[q];

                    aloc[a][0][b][1] += ( lambda*gradphi[q][1][b]*gradphi[q][0][a] + mu*gradphi[q][0][b]*gradphi[q][1][a] ) * w[q];

                    aloc[a][1][b][0] += ( lambda*gradphi[q][0][b]*gradphi[q][1][a] + mu*gradphi[q][1][b]*gradphi[q][0][a] ) * w[q];

                    aloc[a][1][b][1] += ( gradphi[q][1][a]*(lambda*gradphi[q][1][b] + 2.0*mu*gradphi[q][1][b]) + mu*gradphi[q][0][b]*gradphi[q][0][a] ) * w[q];
                }
            }
        }

        for (a = 0; a < nln; a = a + 1 )
        {
            /* loop over test components --> i_c */
            for (i_c = 0; i_c < 2; i_c = i_c + 1 )
            {
                /* loop over trial functions --> b */
                for (b = 0; b < nln; b = b + 1 )
                {
                    /* loop over trial components --> j_c */
                    for (j_c = 0; j_c < 2; j_c = j_c + 1 )
                    {
                        myArows[ie*nln2*4+iii] = elements[a+ie*numRowsElements] + i_c * NumNodes;
                        myAcols[ie*nln2*4+iii] = elements[b+ie*numRowsElements] + j_c * NumNodes;
                        myAcoef[ie*nln2*4+iii] = aloc[a][i_c][b][j_c]*detjac[ie];
                        iii = iii + 1;
                    }
                }
            }
        }
    }

}
/*************************************************************************/

void LinearElasticMaterial_stress(mxArray* plhs[], const mxArray* prhs[])
{

    double* dim_ptr = mxGetPr(prhs[0]);
    int dim     = (int)(dim_ptr[0]);
    int noe     = mxGetN(prhs[4]);
    double* nln_ptr = mxGetPr(prhs[5]);
    int nln     = (int)(nln_ptr[0]);
    int numRowsElements  = mxGetM(prhs[4]);
    int nln2    = nln*nln;

    plhs[0] = mxCreateDoubleMatrix(noe,dim*dim, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(noe,dim*dim, mxREAL);

    double* P        = mxGetPr(plhs[0]);
    double* Sigma    = mxGetPr(plhs[1]);

    int k,l;
    int q;
    int NumQuadPoints     = mxGetN(prhs[6]);
    int NumNodes          = (int)(mxGetM(prhs[3]) / dim);

    double* U_h   = mxGetPr(prhs[3]);
    double* w   = mxGetPr(prhs[6]);
    double* invjac = mxGetPr(prhs[7]);
    double* detjac = mxGetPr(prhs[8]);
    double* phi = mxGetPr(prhs[9]);
    double* gradrefphi = mxGetPr(prhs[10]);

    double gradphi[dim][nln];
    double* elements  = mxGetPr(prhs[4]);

    double GradUh[dim][dim];

    double Id[dim][dim];
    int d1,d2;
    for (d1 = 0; d1 < dim; d1 = d1 + 1 )
    {
        for (d2 = 0; d2 < dim; d2 = d2 + 1 )
        {
            Id[d1][d2] = 0;
            if (d1==d2)
            {
                Id[d1][d2] = 1;
            }
        }
    }

    double F[dim][dim];
    double EPS[dim][dim];

    double* material_param = mxGetPr(prhs[2]);
    double mu = material_param[0];
    double lambda = material_param[1];

    /* Assembly: loop over the elements */
    int ie;

#pragma omp parallel for shared(invjac,detjac,elements,U_h) private(gradphi,F,EPS,GradUh,ie,k,l,q,d1,d2) firstprivate(phi,gradrefphi,w,numRowsElements,nln2,nln,NumNodes,Id,mu,lambda)

    for (ie = 0; ie < noe; ie = ie + 1 )
    {
        q = 0;
        for (k = 0; k < nln; k = k + 1 )
        {
            for (d1 = 0; d1 < dim; d1 = d1 + 1 )
            {
                gradphi[d1][k] = 0;
                for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                {
                    gradphi[d1][k] = gradphi[d1][k] + INVJAC(ie,d1,d2)*GRADREFPHI(k,q,d2);
                }
            }
        }

        for (d1 = 0; d1 < dim; d1 = d1 + 1 )
        {
            for (d2 = 0; d2 < dim; d2 = d2 + 1 )
            {
                GradUh[d1][d2] = 0;
                for (k = 0; k < nln; k = k + 1 )
                {
                    int e_k;
                    e_k = (int)(elements[ie*numRowsElements + k] + d1*NumNodes - 1);
                    GradUh[d1][d2] = GradUh[d1][d2] + U_h[e_k] * gradphi[d2][k];
                }
                F[d1][d2] = Id[d1][d2] + GradUh[d1][d2];
            }
        }

        for (d1 = 0; d1 < dim; d1 = d1 + 1 )
        {
            for (d2 = 0; d2 < dim; d2 = d2 + 1 )
            {
                EPS[d1][d2] = 0.5 * ( F[d1][d2] + F[d2][d1] ) - Id[d1][d2];
            }
        }

        double trace = Trace(dim, EPS);
        /* For linear elasticity, Cauchy and 1st PK stress tensor coincide */
        for (d1 = 0; d1 < dim; d1 = d1 + 1 )
        {
            for (d2 = 0; d2 < dim; d2 = d2 + 1 )
            {
                Sigma[ie+(d1+d2*dim)*noe] = 2 * mu * EPS[d1][d2] + lambda * trace * Id[d1][d2];
                P[ie+(d1+d2*dim)*noe]     = Sigma[ie+(d1+d2*dim)*noe];
            }
        }

    }
}
/*************************************************************************/

void LinearElasticMaterial_stiffness(mxArray* plhs[], const mxArray* prhs[])
{

    double* dim_ptr = mxGetPr(prhs[0]);
    int dim     = (int)(dim_ptr[0]);
    int noe     = mxGetN(prhs[2]);
    double* nln_ptr = mxGetPr(prhs[3]);
    int nln     = (int)(nln_ptr[0]);
    double* NumNodes_ptr = mxGetPr(prhs[9]);
    int NumNodes= (int)(NumNodes_ptr[0]);
    int numRowsElements  = mxGetM(prhs[2]);
    int nln2     = nln*nln;
    int nln2_dim = (nln*dim)*(nln*dim);

    plhs[0] = mxCreateDoubleMatrix(nln2_dim*noe,1, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(nln2_dim*noe,1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(nln2_dim*noe,1, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(nln2_dim*noe,1, mxREAL);

    double* myKrows       = mxGetPr(plhs[0]);
	double* myKcols       = mxGetPr(plhs[1]);
    double* myKcoef_mu    = mxGetPr(plhs[2]);
    double* myKcoef_lambda= mxGetPr(plhs[3]);

    int k, q;
    int NumQuadPoints     = mxGetN(prhs[4]);

    double* w   = mxGetPr(prhs[4]);
    double* invjac = mxGetPr(prhs[5]);
    double* detjac = mxGetPr(prhs[6]);
    double* gradrefphi = mxGetPr(prhs[8]);

    double gradphi[dim][nln][NumQuadPoints];
    double* elements  = mxGetPr(prhs[2]);

    /* Assembly: loop over the elements */
    int ie, i1;
    int iii = 0;
	int a, b;
	int d1, d2, d3;
	double tmp_mu, tmp_lambda;

	double LocalStiff_mu[nln*dim][nln*dim];
	double LocalStiff_lambda[nln*dim][nln*dim];

	i1 = dim +1;
	double dim_circle[i1];
	for (d1 = 0; d1 < dim; d1 = d1 +1 )
	{
		dim_circle[d1] = d1;
	}
	dim_circle[i1] = 0;

	/*FILE * fp; */


	#pragma omp parallel for shared(invjac,detjac,elements,myKrows,myKcols,myKcoef_mu,myKcoef_lambda,dim_circle) private(gradphi,ie,k,q,d1,d2,d3,tmp_mu,tmp_lambda) firstprivate(gradrefphi,w,numRowsElements,nln2,nln2_dim,nln,NumNodes,LocalStiff_mu,LocalStiff_lambda)

	/* fp = fopen("C:/Users/Rosafalco/Documents/Polimi/3_Dottorato/Codici/Matlab/Quarteroni_Manzoni_Negri/pippoLE.out", "w"); */

	/* Assembly: loop over the elements */
    for (ie = 0; ie < noe; ie = ie + 1 )
    {
		/* Compute the entries of the compatibility matrix gradphi */
        for (k = 0; k < nln; k = k + 1 )
        {
            for (q = 0; q < NumQuadPoints; q = q + 1 )
            {
               	for (d1 = 0; d1 < dim; d1 = d1 + 1 )
               	{
                   gradphi[d1][k][q] = 0;
                   for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                   {
                      gradphi[d1][k][q] = gradphi[d1][k][q] + INVJAC(ie,d1,d2)*GRADREFPHI(k,q,d2);
                   }
                   /*fprintf(fp, "indici k=%u d1=%u; gradphi=%f  \n", k, d1, gradphi[d1][k][q]); */
               	}
            }
        }

		/* Compute the entries of the element matrix contribution LocalStiff_mu and LocalStiff_lambda */
		/* a tes, b trial */
		/* loop over the nodes */
		for (a = 0; a < nln; a = a + 1 )
		{
			/* loop over the nodes */
			for (b = 0; b < nln; b = b + 1 )
			{
				/* loop over the directions - 1st term for LocalStiff_mu
				                            - LocalStiff_lambda*/
				for (d1 = 0; d1 < dim; d1 = d1 + 1 )
				{
			    	for (d2 = 0; d2 < dim; d2 = d2 + 1 )
					{
						tmp_mu = 0;
						tmp_lambda = 0;
						for (q = 0; q < NumQuadPoints; q = q + 1 )
	                	{
							tmp_mu = tmp_mu + gradphi[d1][b][q] * gradphi[d2][a][q]* w[q];
							tmp_lambda = tmp_lambda + gradphi[d1][a][q] * gradphi[d2][b][q]* w[q];
						}
						LocalStiff_mu[dim*a+d1][dim*b+d2] = tmp_mu;
						LocalStiff_lambda[dim*a+d1][dim*b+d2] = tmp_lambda;
						if (d1 != d2)
						{
							LocalStiff_mu[dim*a+d1][dim*b+d2] = LocalStiff_mu[dim*a+d1][dim*b+d2] * 0.5;
						}
					}
				}
				/* loop over the directions - 2nd term for LocalStiff_lambda*/
				for (d1 = 0; d1 < dim; d1 = d1 + 1 )
				{
					for (d2 = 1; d2 < dim; d2 = d2 +1 )
					{
						i1 = d1 + d2;
			    		d3 = dim_circle[i1];
						tmp_mu = 0;
						for (q = 0; q < NumQuadPoints; q = q + 1 )
	                	{
							tmp_mu = tmp_mu + gradphi[d1][a][q] * gradphi[d1][b][q]* w[q];
						}
						LocalStiff_mu[dim*a+d3][dim*b+d3] = LocalStiff_mu[dim*a+d3][dim*b+d3] + 0.5 * tmp_mu;
					}
				}
			}
		}

		iii = 0;

		/* Assembly*/
		/* a tes, b trial */
        for (a = 0; a < nln; a = a + 1 )
        {
            for (b = 0; b < nln; b = b + 1 )
            {
				for (d1 = 0; d1 < dim; d1 = d1 + 1 )
				{
					for (d2 = 0; d2 < dim; d2 = d2 + 1 )
					{
						myKrows[ie*nln2_dim+iii] = elements[a+ie*numRowsElements]+d1*NumNodes;
						myKcols[ie*nln2_dim+iii] = elements[b+ie*numRowsElements]+d2*NumNodes;
		    	   		myKcoef_mu[ie*nln2_dim+iii]     = LocalStiff_mu[dim*a+d1][dim*b+d2]*detjac[ie];
		    	   		myKcoef_lambda[ie*nln2_dim+iii] = LocalStiff_lambda[dim*a+d1][dim*b+d2]*detjac[ie];
	                	iii = iii + 1;
	            	}
				}
        	}
		}
	}
}
/*fclose(fp);*/
/*************************************************************************/