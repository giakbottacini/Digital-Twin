/*   This file is part of redbKIT.
 *   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
 *   Author: Federico Negri <federico.negri@epfl.ch>
 */

#include "mex.h"
#include <stdio.h>
#include <math.h>
#include "blas.h"
#include <string.h>
#define INVJAC(i,j,k) invjac[i+(j+k*dim)*noe]
#define GRADREFPHI(i,j,k) gradrefphi[i+(j+k*NumQuadPoints)*nln]

#include "../../Core/Tools.h"
#include "MaterialModels/LinearElasticMaterial.h"
#include "MaterialModels/LinearElasticMaterialSpatial.h"
#include "MaterialModels/SEMMTMaterial.h"
#include "MaterialModels/NeoHookeanMaterial.h"
#include "MaterialModels/StVenantKirchhoffMaterial.h"
#include "MaterialModels/RaghavanVorpMaterial.h"

#ifdef _OPENMP
#include <omp.h>
#else
#warning "OpenMP not enabled. Compile with mex CSM_assembler_C_omp.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp""
#endif

/*************************************************************************/

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{

    /* Check for proper number of arguments. */

    /*
    if(nrhs!=11) {
        mexErrMsgTxt("11 inputs are required.");
    } else if(nlhs>6) {
        mexErrMsgTxt("Too many output arguments.");
    }
    */

	/*NOTA: in questa fase non � importante che cosa esattamente gli stai passando, ma solo il numero di entrate*/

    char *Material_Model = mxArrayToString(prhs[1]);
    double* dim_ptr = mxGetPr(prhs[0]);
    int dim     = (int)(dim_ptr[0]);  /*dim � la dimensione spaziale della mesh =2 -> piana =3 -> spaziale*/


    if (strcmp(Material_Model, "Linear_stiffness")==0)
    {
            LinearElasticMaterial_stiffness(plhs, prhs);
    }

    if (strcmp(Material_Model, "Linear_forces")==0)
    {
            LinearElasticMaterial_forces(plhs, prhs);
    }

    if (strcmp(Material_Model, "Linear_jacobianSlow")==0)
    {
            LinearElasticMaterial_jacobian(plhs, prhs);
    }

    if (strcmp(Material_Model, "Linear_jacobian")==0)
    {
        if (dim == 2)
        {
            LinearElasticMaterial_jacobianFast2D(plhs, prhs);
        }

        if (dim == 3)
        {
            LinearElasticMaterial_jacobianFast3D(plhs, prhs);
        }
    }

    if (strcmp(Material_Model, "Linear_stress")==0)
    {
            LinearElasticMaterialSpatial_stress(plhs, prhs);
    }

	if (strcmp(Material_Model, "LinearSpatial_stiffness")==0)
    {
            LinearElasticMaterialSpatial_stiffness(plhs, prhs);
    }

    if (strcmp(Material_Model, "LinearSpatial_forces")==0)
    {
            LinearElasticMaterialSpatial_forces(plhs, prhs);
    }

    if (strcmp(Material_Model, "LinearSpatial_jacobianSlow")==0)
    {
            LinearElasticMaterialSpatial_jacobian(plhs, prhs);
    }

    if (strcmp(Material_Model, "LinearSpatial_jacobian")==0)
    {
        if (dim == 2)
        {
            LinearElasticMaterialSpatial_jacobianFast2D(plhs, prhs);
        }

        if (dim == 3)
        {
            LinearElasticMaterialSpatial_jacobianFast3D(plhs, prhs);
        }
    }

    if (strcmp(Material_Model, "LinearSpatial_stress")==0)
    {
            LinearElasticMaterialSpatial_stress(plhs, prhs);
    }


    if (strcmp(Material_Model, "SEMMT_forces")==0)
    {
            SEMMTMaterial_forces(plhs, prhs);
    }

    if (strcmp(Material_Model, "SEMMT_jacobianSlow")==0)
    {
            SEMMTMaterial_jacobian(plhs, prhs);
    }

    if (strcmp(Material_Model, "SEMMT_jacobian")==0)
    {
        if (dim == 2)
        {
            SEMMTMaterial_jacobianFast2D(plhs, prhs);
        }

        if (dim == 3)
        {
            SEMMTMaterial_jacobianFast3D(plhs, prhs);
        }
    }


    if (strcmp(Material_Model, "StVenantKirchhoff_forces")==0)
    {
            StVenantKirchhoffMaterial_forces(plhs, prhs);
    }

    if (strcmp(Material_Model, "StVenantKirchhoff_jacobianSlow")==0)
    {
            StVenantKirchhoffMaterial_jacobian(plhs, prhs);
    }

    if (strcmp(Material_Model, "StVenantKirchhoff_jacobian")==0)
    {
        if (dim == 2)
        {
            StVenantKirchhoffMaterial_jacobianFast2D(plhs, prhs);
        }

        if (dim == 3)
        {
            StVenantKirchhoffMaterial_jacobianFast3D(plhs, prhs);
        }

    }

    if (strcmp(Material_Model, "StVenantKirchhoff_stress")==0)
    {
            StVenantKirchhoffMaterial_stress(plhs, prhs);
    }

    if (strcmp(Material_Model, "NeoHookean_forces")==0)
    {
            NeoHookeanMaterial_forces(plhs, prhs);
    }

    if (strcmp(Material_Model, "NeoHookean_jacobian")==0)
    {
            NeoHookeanMaterial_jacobianFast(plhs, prhs);
    }

    if (strcmp(Material_Model, "NeoHookean_jacobianSlow")==0)
    {
            NeoHookeanMaterial_jacobian(plhs, prhs);
    }

    if (strcmp(Material_Model, "NeoHookean_stress")==0)
    {
            NeoHookeanMaterial_stress(plhs, prhs);
    }

    if (strcmp(Material_Model, "NeoHookean_prestress")==0)
    {
            NeoHookeanMaterial_prestress(plhs, prhs);
    }

    if (strcmp(Material_Model, "RaghavanVorp_forces")==0)
    {
            RaghavanVorpMaterial_forces(plhs, prhs);
    }

    if (strcmp(Material_Model, "RaghavanVorp_jacobian")==0)
    {
            RaghavanVorpMaterial_jacobianFast(plhs, prhs);
    }

    if (strcmp(Material_Model, "RaghavanVorp_jacobianSlow")==0)
    {
            RaghavanVorpMaterial_jacobian(plhs, prhs);
    }

    if (strcmp(Material_Model, "RaghavanVorp_stress")==0)
    {
            RaghavanVorpMaterial_stress(plhs, prhs);
    }

    mxFree(Material_Model);

}
/*************************************************************************/

