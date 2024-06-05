/**
 * g_ordercg calculate the order parameter for lipids in a coarse grained model.
 * 
 * See the description at the end of the program to get more information.
 **/

#include "string.h"

#include <gromacs/copyrite.h>
#include <gromacs/gmx_fatal.h>
#include <gromacs/pbc.h>
#include <gromacs/rmpbc.h>
#include <gromacs/smalloc.h>
#include <gromacs/statutil.h>
#include <gromacs/vec.h>
#include <gromacs/xvgr.h>
#include <gromacs/futil.h>

#include "distances.h"

#ifdef GMX_DOUBLE
real ceilr(double val) {return ceil(val);}
#else
real ceilr(float val) {return ceilf(val);}
#endif

static const char *authors[] = {
    "Written by Jonathan Barnoud (jonathan.barnoud@inserm.fr)",
    "Copyright (c) 2012  Jonathan Barnoud, Luca Monticelli"
};
static const char *gpl[] = {
    "This program is free software; you can redistribute it and/or",
    "modify it under the terms of the GNU General Public License",
    "as published by the Free Software Foundation; either version 2",
    "of the License, or (at your option) any later version."
};
static const char *version[] = {
    ":-) g_ordercg - version 1.0 (-:"
};
void sp_print(FILE *out, const char *txt)
{
    int i, s;
    s = (int) (80 - strlen(txt)) / 2.0;
    for (i=0; i < s; i++) 
        fprintf(out, " ");
    fprintf(out, "%s\n", txt);
}


typedef struct it_status {
    int nmols;
    int nbonds;
    int idx_mol;
    int idx_bond;
    int **def_bonds;
} it_status;

it_status *init_it_status(int nmols, int nbonds, int **def_bonds) {
    it_status *status = NULL;
    snew(status, 1);
    status->nmols = nmols;
    status->nbonds = nbonds;
    status->idx_mol = 0;
    status->idx_bond = -1;
    status->def_bonds = def_bonds;
    return status;
}

/* Iterate on bonds when a topology is provided
 *
 * Parameters :
 *  - x1, x2 : where the coordinates of the bond extremities will be returned
 *  - idx : where the storage index will be returned, this index is the bond
 *      index in this case
 *  - x : the coordinates of all the atoms in the frame
 *  - index : the index of the atoms of interest in the coordinate vector, the
 *      index group should contein the first atom of each molecule of interest
 *  - status : the iterator
 * 
 * Return :
 *  0 until the iterations are done, 1 after that
 * */
int next_bond_mol(rvec **x1, rvec **x2, int *idx, rvec *x, atom_id **index,
  it_status *status) {
    /* Update indices */
    (status->idx_bond)++;
    if (status->idx_bond >= status->nbonds) {
        status->idx_bond = 0;
        (status->idx_mol)++;
    }
    /* Are we at the end of the iterator ? */
    if (status->idx_mol >= status->nmols) {
        return 1;
    }
    /* Return the values */
    *x1 = &(x[index[0][status->idx_mol]+(status->def_bonds)[status->idx_bond][0]]);
    *x2 = &(x[index[0][status->idx_mol]+(status->def_bonds)[status->idx_bond][1]]);
    *idx = status->idx_bond;
    return 0;
}

/* Iterate on bonds in vectors mode.
 *
 * Parameters :
 *  - x1, x2 : where the coordinates of the bond extremities will be returned
 *  - idx : where the storage index will be returned, this index is the bond
 *      index in this case
 *  - x : the coordinates of all the atoms in the frame
 *  - index : the index of the atom of interest in the coordinate vector, the
 *      atoms should be the extremities of each vector of interect
 *  - status : the iterator
 * 
 * Return :
 *  0 until the iterations are done, 1 after that
 * */
int next_bond_vec(rvec **x1, rvec **x2, int *idx, rvec *x, atom_id **index,
  it_status *status) {
    (status->idx_bond)++;
    if (status->idx_bond >= status->nbonds) {
        return 1;
    }
    *x1 = &x[index[0][(status->idx_bond)*2]];
    *x2 = &x[index[0][(status->idx_bond)*2]+1];
    *idx = status->idx_bond;
    return 0;
}

/* Reset the iterator */
void reset_iterator(it_status *status) {
    status->idx_bond = -1;
    status->idx_mol = 0;
}

/*****************************************************************************
 *                            Get the topology                               *
 *****************************************************************************/

/** Integer comparison function for qsort */
int compare (const void * a, const void * b) {
  return ( *(int*)a - *(int*)b );
}

/** Get the index of the first atom of each molecule present in an index group
 *
 * Take a sorted index group and build a list with the index of the first atom
 * of each molecule present in the group. Each molecule will appear only once
 * in the output table even if it appears more in the input group.
 *
 * :Parameters:
 *     - index  : an index table
 *     - nindex : the length of the index
 *     - top    : a pointer to the topology
 *     - molindex : a pointer to the output table
 *     - nmol   : a pointer to the number of found molecules
 */
void first_mol_atom(int *index, int nindex, t_topology *top,
        int **molindex, int *nmol) {
    int imol=0; /* current index in top->mols.index */
    int idx=0;  /* current index in the index table */

    /* molindex maximum size is nindex. We will reduce its size later when we
     * will know how manu molecules we have. */
    snew((*molindex), nindex);
    *nmol=0;

    for (idx=0; idx<nindex; ++idx) {
        /* Go forward to the right molecule if needed. */
        while ((top->mols).index[imol+1] <= index[idx]) {
            imol++;
        }
        /* Is this atom in a new molecule? */
        if (*nmol == 0 || top->mols.index[imol] > (*molindex)[(*nmol)-1]) {
            (*molindex)[(*nmol)] = top->mols.index[imol];
            (*nmol)++;
        }
    }
    
    /* Now we know how many molecules we have so we can adjust the length of
     * the output table. */
    srenew((*molindex), *nmol);
}


/** Generate a non-redundant list of residues of interest from an index group
 * and check if all the residue have the same name
 *
 * :Parameters:
 *     - index   : an index table
 *     - nindex  : the lenght of the index
 *     - top     : a pointer to the topology
 *     - res     : a pointer to the table of the first atom of each residue of
 *                interest. The table will be allocated.
 *     - resname : a pointer to the residue name
 *
 * :Return:
 *     - status : 0 if everything went right and if all the residues have the
 *                same name
 **/
int get_res(int *index, int nindex, t_topology *top, int **res,
  const char **resname);
        

/** Get the bonds for a given residues type relatively to the first atom of the
 * residue
 *
 * Take a residue id and build the topology of the bonds of the given residue
 * relatively to its first atom.
 *
 * The function also build a legend. One entry corresponds to the bond with the
 * same index in the bonds table and take the form "AT1-AT2\0" with AT1 and AT2
 * the atom type of the first and second atom forming the bond respectively.
 * The table contains an extra entry which contains "Total\0" so the length of
 * the legend table is nbonds+1.
 *
 * :Parameters:
 *     - atomid : the index of one atom in the selection used for
 *                order parameter calculation
 *     - top    : a pointer to the topology
 *     - subtop : a pointer to the 2D table which will be created
 *     - legend : a pointer to the legend table which will be created
 *     - nbonds : a pointer to the number of bonds
 *
 * :Return:
 *      - status : 0 if everything went right, 1 else
 **/
int subtopol(int atomid, t_topology *top, int ***subtop, char ***legend,
  int *nbonds) {
    int bond, first, start=0, i=0, j=0, k=0;
    int *single_bond = NULL;
    char *single_leg = NULL;
    int legend_index;
    int imol = 0;

    /* Get the index of the first atom of the residue */
    while (top->mols.index[imol+1] <= atomid) {
        imol++;
    }
    first = top->mols.index[imol];

    /* Get the number of bonds */
    *nbonds = 0;
    start = -1;
    for (bond=0; bond < top->idef.il[0].nr; bond+=3) {
        if (top->idef.il[0].iatoms[bond+1] >= top->mols.index[imol] &&
                top->idef.il[0].iatoms[bond+1] < top->mols.index[imol+1]) {
            (*nbonds)++;
            if (start < 0) start = bond;
        }
    }

    /* Store the bonds and create the legend */
    *subtop = (int**)malloc((*nbonds)*sizeof(int*));
    *legend = (char**)malloc(((*nbonds)+1)*sizeof(char*));
    if (!(subtop || legend)) {
        gmx_fatal(FARGS, "Allocation error.\n");
    }
    snew(single_leg, 8);
    (*legend)[*nbonds] = single_leg;
    strcpy((*legend)[(*nbonds)], "Total\0"); /* Add the "Total" entry */
    i = 0;
    for (bond=0; bond<top->idef.il[0].nr; bond+=3) {
        if (top->idef.il[0].iatoms[bond+1] >= top->mols.index[imol] &&
                top->idef.il[0].iatoms[bond+1] < top->mols.index[imol+1]) {
            snew(single_bond, 2);
            snew(single_leg, 8);
            (*subtop)[i] = single_bond;
            (*legend)[i] = single_leg;
            legend_index = -1;
            for (j=0; j<2; j++) {
                (*subtop)[i][j] = top->idef.il[0].iatoms[bond+j+1] - first;
                for (k=0; k<3; k++) {
                    if ((*(top->atoms.atomname[top->idef.il[0].iatoms[bond+j+1]]))[k] != '\0') {
                        legend_index++;
                        (*legend)[i][legend_index] = \
                            (*(top->atoms.atomname[top->idef.il[0].iatoms[bond+j+1]]))[k];
                    }
                }
                if (j == 0) {
                    legend_index++;
                    (*legend)[i][legend_index] = '-';
                }
                else {
                    legend_index++;
                    (*legend)[i][legend_index] = '\0';
                }
            }
            printf("bond: %3d %3d %s\n",
                (*subtop)[i][0]+1, (*subtop)[i][1]+1, (*legend)[i]);
            i++;
        }
    }

    return 0;
}

/** Read in index file and interpret it as list of bonds relatively to the
 * first atom of a molecule.
 *
 * Atoms in the index file are read by pairs wich represent the extremity of a
 * bond with respect to the first atom of the molecule.
 *
 * A dummy legend is also generated.
 *
 * :Parameter:
 *     - index : index to read
 *     - isize : the lenth of the index table
 *     - subtop : a pointer to the 2D table which will be created
 *     - legend : a pointer to the legend table which will be created
 *     - nbonds : a pointer to the number of bonds
 *
 * :Return:
 *      - status : 0 if everything went right, 1 else
 **/
int read_user_subtopol(atom_id *index, int isize, int ***subtop, char ***legend,
  int *nbonds) {
    int *single_bond = NULL;
    char *single_leg = NULL;
    int i=0;

    /* Check if the index file corresponds to the requirements */
    if (isize % 2 != 0) {
        fprintf(stderr, "User bond list has to be an index file with an even"
                "number of items.\n");
        return(1);
    }
    /* Create the data structures to save the result */
    *nbonds = isize/2;
    *subtop = (int**)malloc((*nbonds)*sizeof(int*));
    *legend = (char**)malloc(((*nbonds)+1)*sizeof(char*));
    if (!(subtop || legend)) {
        gmx_fatal(FARGS, "Allocation error.\n");
    }
    snew(single_leg, 8);
    (*legend)[*nbonds] = single_leg;
    strcpy((*legend)[(*nbonds)], "Total\0"); /* Add the "Total" entry */

    for (i=0; i<(*nbonds); i+=1) {
        snew((*subtop)[i], 2);
    }

    /* Do the actual stuff */
    /* Note that the index read from the index file are minus 1 */
    for (i=0; i<isize; i+=2) {
        (*subtop)[i/2][0] = index[i];
        (*subtop)[i/2][1] = index[i+1];
        snew((*legend)[i/2], 80);
        sprintf((*legend)[i/2],"%i-%i", index[i], index[i+1]);
    }
    return(0);
}

/* TODO sanity check on user bond lists*/

/*****************************************************************************
 *                             Density stuff                                 *
 *****************************************************************************/

/** Get the coordinate of a bond center
 *
 * :Parameters:
 *     - pbc  : the periodic box
 *     - xi   : the coordinates of the first atom
 *     - xj   : the coordinates of the second one
 *     - center : the vector where to write the coordinate of the bound center
 **/
void get_bound_center(const t_pbc *pbc, rvec xi, rvec xj, rvec center) {
    int i;
    for (i = 0; i<DIM; i++) {
        center[i] = (xi[i]+xj[i])/2;
    }
    put_atom_in_box((real (*)[3])pbc->box,center);
}

/*****************************************************************************
 *                              Grid stuff                                   *
 *****************************************************************************/
/** Create a matrix of real numbers
 *
 * The matrix is pre-filled with a value.
 *
 * Parameters:
 *  - d1: number of rows
 *  - d2: number of columns
 *  - defval: the value to put in every cells
 *
 * Return:
 *  The filled matrix.
 */
real **realMatrix(int d1, int d2, real defval) {
    int i, j;
    real **mat;
  
    smalloc(mat, d1 * sizeof(real*));
    for(i = 0; i < d1; i++){
        smalloc(mat[i], d2 * sizeof(real));
        for(j = 0; j < d2; j++) {
            mat[i][j] = defval;
        }
    }
    return mat;
}

/** Destroy a matrix of real numbers
 *
 * Parameters:
 *  - mat: the matrix to destroy
 *  - d1: the number of rows in the matrix
 */
void deleteRealMat(real **mat, int d1) {
    int i;
    for(i = 0; i < d1; i++) {
        sfree(mat[i]);
    }
    sfree(mat);
}

/** Create a matrix of integers
 *
 * The matrix is pre-filled with a value.
 *
 * Parameters:
 *  - d1: number of rows
 *  - d2: number of columns
 *  - defval: the value to put in every cells
 *
 * Return:
 *  The filled matrix.
 */
int **intMatrix(int d1, int d2, int defval){
    int i, j;
    int **mat;
  
    smalloc(mat, d1 * sizeof(int*));
    for(i = 0; i < d1; i++){
        smalloc(mat[i], d2 * sizeof(real));
        for(j = 0; j < d2; j++) {
            mat[i][j] = defval;
        }
    }
    return mat;
}

/** Destroy a matrix of integers
 *
 * Parameters:
 *  - mat: the matrix to destroy
 *  - d1: the number of rows in the matrix
 */
void deleteIntlMat(int **mat, int d1) {
    int i;
    for(i = 0; i < d1; i++) {
        sfree(mat[i]);
    }
    sfree(mat);
}


/*****************************************************************************
 *                             Index stuff                                   *
 *****************************************************************************/
/** Read index file for distance calculation
 *
 * Parameters:
 *  - index_fn: index file name
 *  - top: the gromacs topology
 *  - index: table that will contain indices
 *  - isize: will contain the number of indices in the index table
 *  - grpname: will contain the selected group name
 */
void get_dist_index(const char *index_fn, t_topology *top,
        atom_id ***index, int **isize, char ***grpname) {
    printf("Select the reference for distance calculation:\n");
    snew(*grpname, 1);
    snew(*index,1);
    snew(*isize,1);
    get_index(&top->atoms,index_fn,1,*isize,*index,*grpname);
}

/** Read index groups for order parameter calculation
 */
void get_order_index(gmx_bool vectors, gmx_bool bFirst, const char *index_fn,
        t_topology *top, atom_id ***index, int **isize, char ***grpname) {
    if (!vectors) {
        if (bFirst) {
            printf("Select a group for computation:\n");
        }
        else {
            printf("Select a group with only the first atom for each molecule:\n");
        }
    }
    else {
        printf("Select a group with a list of bonds:\n");
    }
    snew(*grpname,1);
    snew(*index,1);
    snew(*isize,1);
    get_index(&top->atoms,index_fn,1,*isize,*index,*grpname);

    if (vectors && (*isize)[0] % 2 != 0) {
        gmx_fatal(FARGS, "If you do not use automatic bond search, "
                "the index file has to contain an even number of atoms.\n");
    }
}

void set_topology(gmx_bool vectors, gmx_bool bFirst, const char *fnUBL,
        atom_id **index, int *isize, char **grpname, t_topology *top,
        int ***def_bonds, char ***leg, int *n_angles) {
    char    **ugrpname; /* the name of each group */
    int     *uisize;    /* the size of each group */
    atom_id **uindex;   /* the index for the atom numbers */
    atom_id *first_index = NULL;
    int nmol;

    if (!vectors) {
        if (bFirst) {
            /* Fix the index group */
            qsort(index[0], isize[0], sizeof(int), compare);
            first_mol_atom(index[0], isize[0], top, &first_index, &nmol);
            printf("Reading group %s containing %d atoms: "
                    "%d molecules found\n", grpname[0], isize[0], nmol);
            index[0] = first_index;
            isize[0] = nmol;
        }
        /* Get the bonds relatively to the first atom of the residue. */
        if (fnUBL != NULL) {
            /* Get the bonds from user input */
            printf("Select a group with molecule bonds decription:\n");
            snew(ugrpname, 1);
            snew(uindex, 1);
            snew(uisize, 1);
            get_index(&top->atoms,fnUBL,1,uisize,uindex,ugrpname);
            read_user_subtopol(uindex[0], uisize[0], def_bonds, leg, n_angles);
        } else {
            /* Get the bonds from topology */
            subtopol(index[0][0], top, def_bonds, leg, n_angles);
        }
    }
}

/*****************************************************************************
 *                               Main work                                   *
 *****************************************************************************/

/** Compute the order parameter
 *
 * :Parameters:
 *     - sum : the sum of the squared cosines
 *     - pop : the number of records
 *
 * :Return:
 *      - the order parameter
 */
real order_parameter(real sum, real pop) {
    return 0.5*(3*sum/pop-1);
}

/** Do the main job
 *
 * :Parameters:
 *     - fnTRX : the trajectory file name
 *     - fnTPX : the topology file name
 *     - fnNDX : the index file name
 *     - fnXVG : the output file name for grid mode
 *     - fnProf : the output file name for profile mode
 *     - fnGrid : the output file name for grid mode
 *     - fnSampling : the name of the file where to write the sampling in
 *               profile mode if any, NULL else
 *     - fnSamplingGrid : the name of the file where to write the sampling in
 *               grid mode if any, NULL else
 *     - fnUBL : the user bond list file name if any, NULL else
 *     - oenv  : don't really know what it is but it is needed and defined 
 *               during the reading of the arguments
 *     - axis  : the axis normal to the membrane (0 for x, 1 for y and 2 for z)
 *     - nslices : the number of slices. If sl is greater than 0 the order
 *               parameter will be computed as a density on the given axis
 *               else it will be computed classiquely as an average over time on
 *               equivalent bonds
 *     - vectors : if TRUE, the fnNDX index file is read as a list of bonds
 *     - bFirst: if TRUE then the first atom of each molecule is extracted from
 *               the group index and use instead
 *     - bNorm : is the normal mode active?
 *     - bProf : is the profile mode active?
 *     - bGrid : is the grid mode active?
 *
 * This function doesn't return anything but it writes the output files.
 **/
void do_cgorder(const char *fnTRX, const char *fnTPX, const char *fnNDX,
  const char *fnXVG,
  const char *fnProf, const char *fnGrid, const char *fnDistprof, 
  const char *fnSampling, const char *fnSamplingGrid,
  const char *fnSamplingDistprof,
  const char *fnUBL, const output_env_t oenv, int axis, int oaxis,
  int nslices, int nslices2, gmx_bool vectors, gmx_bool bFirst,
  gmx_bool bNorm, gmx_bool bProf, gmx_bool bGrid, gmx_bool bDist, int d2axis) {
    /* Vectors */
    rvec r_u, r_v;
    
    /* Topology */
    t_topology *top=NULL;
    int **def_bonds = NULL;
    char **leg = NULL;
    int n_angles = 0;

    /* Periodic bound conditions */
    int  ePBC;
    t_pbc   *pbc;
    gmx_rmpbc_t  gpbc=NULL;
    matrix box;
    
    /* Output */
    FILE *onorm = NULL;
    FILE *oprof = NULL;
    FILE *ogrid = NULL;
    FILE *odist = NULL;
    FILE *outSampling = NULL;
    FILE *outSamplingGrid = NULL;
    FILE *outSamplingDistprof = NULL;
    char axis_labels[3][2] = {"X\0", "Y\0", "Z\0"};
    real lo, hi, extra;
    real lo_samp, hi_samp;
    t_rgb rlo, rhi, rextra;
    real *ticks_x=NULL, *ticks_y=NULL;
    
    /* Frames */
    rvec *x = NULL;
    int natoms;
    t_trxstatus *status;
    real t;
    int lipid, angle, i, j;
    real *sum_cos2_theta_norm = NULL;
    real *sum_cos2_theta_prof = NULL;
    real *sum_cos2_theta_distprof = NULL;
    int *count_cos2_theta_prof = NULL;
    int *count_cos2_theta_distprof = NULL;
    real **sum_cos2_theta_grid = NULL;
    real **count_cos2_theta_grid = NULL;
    real summation = 0;
    real cosine = 0;
    real sl_width = 0;
    real sl_width_dist = 0;
    int slax1 = 0, slax2 = 0;
    real sl_width_ax1 = 0, sl_width_ax2 = 0;
    real spacing1_sum = 0, spacing2_sum = 0;
    int ax1=0, ax2=0;
    int slice;
    rvec center;
    int nframes = 0;
    gmx_bool profs[3] = {bProf, bDist, bGrid};
    FILE **oprofs[3] = {&oprof, &odist, &ogrid};
    int nprofs = asize(profs);

    /* index reading */
    char    **grpname; /* the name of each group */
    int     *isize;    /* the size of each group */
    atom_id **index;   /* the index for the atom numbers */

    /* reference for distance calculation */
    char    **dgrpname=NULL; /* the name of each group */
    int     *disize=NULL;    /* the size of each group */
    atom_id **dindex=NULL;   /* the index for the atom numbers */

    /* Distance related variables */
    real mass=0;
    rvec *com=NULL;
    
    /* Averaging */
    int angle_per_col = 0;
    int angle_per_row = 0;

    /* Iterator */
    rvec *x1=NULL, *x2=NULL;
    int idx;
    it_status *itstatus=NULL;
    int (*next_bond)(rvec**, rvec**, int*, rvec*, atom_id**, it_status*);

    /* Read the topology */
    top=read_top(fnTPX,&ePBC);

    /* Read the index file */
    if (bDist) {
        get_dist_index(fnNDX, top, &dindex, &disize, &dgrpname);
        mass = get_mass(dindex[0], disize[0], top);
    }
    get_order_index(vectors, bFirst, fnNDX, top, &index,&isize, &grpname);
    set_topology(vectors, bFirst, fnUBL, index, isize, grpname, top,
            &def_bonds, &leg, &n_angles);
    if (!vectors) {
        angle_per_col = isize[0];
        angle_per_row = isize[0]*n_angles;
    } else {
        n_angles = isize[0]/2;
        angle_per_col = 1;
        angle_per_row = n_angles;
    }

    /* Set PBC stuff */
    natoms=read_first_x(oenv,&status,fnTRX,&t,&x,box);
    if (ePBC != epbcNONE)
        snew(pbc,1);
    else
        pbc = NULL;
    gpbc = gmx_rmpbc_init(&top->idef,ePBC,natoms,box);

    /* Prepare output */
    if (bNorm) {
        /* In classic mode */
        onorm = xvgropen(fnXVG, "Order parameter", "Time (ps)",
            "Order parameter", oenv);
        fprintf(onorm, "# You selected the group \"%s\".\n", grpname[0]);
        if (!vectors) {
            fprintf(onorm, "# The group contains %d molecules.\n", isize[0]);
            xvgr_legend(onorm,n_angles+1,(const char**)leg,oenv);
        }
    }
    if (bGrid) {
        ogrid = ffopen(fnGrid, "w");
        if (ogrid == NULL) {
            gmx_fatal(FARGS, "Error opening %s.\n", fnGrid);
        }
        if (fnSamplingGrid != NULL) {
            outSamplingGrid = ffopen(fnSamplingGrid, "w");
            if (outSamplingGrid == NULL) {
                gmx_fatal(FARGS, "Error opening %s.\n", fnSamplingGrid);
            }
        }
    }
    if (bDist) {
        odist = xvgropen(fnDistprof, "Order parameter", "Distance (nm)",
            "Order parameter", oenv);
    }
    if(bProf) {
        /* In density mode */
        oprof = xvgropen(fnProf, "Order parameter", axis_labels[oaxis],
            "Order parameter", oenv);
    }
    for (i=0; i<nprofs; ++i) {
        if (profs[i]) {
            fprintf(*oprofs[i],
                    "# You selected the group \"%s\".\n", grpname[0]);
            if (!vectors) {
                fprintf(*oprofs[i], "# The topology is :\n");
                for (j=0; j<n_angles; ++j) {
                    fprintf(*oprofs[i], "# - %s\n", leg[j]);
                }
            }
        }
    }
    if (fnSampling && bProf) {
        outSampling = xvgropen(fnSampling, "Sampling per slice",
            axis_labels[oaxis], "Number of records", oenv);
    }
    if (fnSamplingDistprof && bDist) {
        outSamplingDistprof = xvgropen(fnSamplingDistprof, "Sampling per slice",
            "Distance (nm)", "Number of records", oenv);
    }
    
    /* Prepare result storage */
    if (bNorm) {
        snew(sum_cos2_theta_norm, n_angles);
    }
    if (bProf) {
        snew(sum_cos2_theta_prof, nslices);
        snew(count_cos2_theta_prof, nslices);
        for (slice=0; slice<nslices; slice++) {
            sum_cos2_theta_prof[slice] = 0;
            count_cos2_theta_prof[slice] = 0;
        }
    }
    if (bDist) {
        snew(sum_cos2_theta_distprof, nslices);
        snew(count_cos2_theta_distprof, nslices);
        for (slice=0; slice<nslices/2; slice++) {
            sum_cos2_theta_distprof[slice] = 0;
            count_cos2_theta_distprof[slice] = 0;
        }
    }
    if (bGrid) {
        sum_cos2_theta_grid = realMatrix(nslices, nslices2, 0.0);
        count_cos2_theta_grid = realMatrix(nslices, nslices2, 0);
        switch (axis) {
            case 0:
                ax1 = 1; ax2 = 2;
                break;
            case 1:
                ax1 = 0; ax2 = 2;
                break;
            case 2:
                ax1 = 0; ax2 = 1;
                break;
            default:
                gmx_fatal(FARGS,"Invalid axes. Terminating. \n");
        }
    }

    /* Define the membrane normal vector */
    for (i=0; i<DIM; i++) {
        r_v[i] = 0;
    }
    r_v[axis] = 1;

    /* Set the iterator */
    if (!vectors) {
        itstatus = init_it_status(isize[0], n_angles, def_bonds);
        next_bond = next_bond_mol;
    }
    else {
        itstatus = init_it_status(0, isize[0]/2, def_bonds);
        next_bond = next_bond_vec;
    }

    /* Read the trajectory */
    do {
        nframes++;
        /* Handle periodic boundary conditions */
        if (pbc) {
            set_pbc(pbc,ePBC,box);
            /* make molecules whole again */
            gmx_rmpbc(gpbc,natoms,box,x);
        }
        
        /* Reset cosine summation in classic mode*/
        if (bNorm) {
            for (angle=0; angle<n_angles; angle++) {
                sum_cos2_theta_norm[angle] = 0;
            }
            summation = 0;
        }
        /* Reset slice width in profile mode */
        if (bProf) {
            sl_width = box[oaxis][oaxis]/nslices;
        }
        /* Reset cell size in grid mode */
        if (bGrid) {
            sl_width_ax1 = box[ax1][ax1]/nslices;
            sl_width_ax2 = box[ax2][ax2]/nslices2;
            spacing1_sum += box[ax1][ax1];
            spacing2_sum += box[ax2][ax2];
        }

        /* Get the center of mass of the distance reference group */
        if (bDist) {
            oaxis=0;
            com = center_of_mass(dindex[0], disize[0], x, top, mass);
            sl_width_dist = box[oaxis][oaxis]/nslices;
        }

        /* Loop on bonds */
        reset_iterator(itstatus);
        while (!next_bond(&x1, &x2, &idx, x, index, itstatus)) {
            /* Get vectors */
            pbc_rvec_sub(pbc, *x1, *x2, r_u);
            /* Update cosine summation */
            cosine = cos_angle(r_u, r_v);
            if (bProf || bGrid || bDist) get_bound_center(pbc,*x1,*x2, center);
            if (bNorm) {
                sum_cos2_theta_norm[idx] += (cosine*cosine);
            }
            if (bProf) {
                slice = center[oaxis]/sl_width;
                sum_cos2_theta_prof[slice] += (cosine*cosine);
                count_cos2_theta_prof[slice]++;
            }
            if (bGrid) {
                slax1 = (int)((double)center[ax1]/(double)sl_width_ax1);
                slax2 = (int)((double)center[ax2]/(double)sl_width_ax2);
                sum_cos2_theta_grid[slax1][slax2] += (cosine*cosine);
                count_cos2_theta_grid[slax1][slax2]++;
            }
            if (bDist) {
                if (d2axis >= 0) {
                    center[d2axis] = 0;
                    (*com)[d2axis] = 0;
                }
                slice = distance(center,*com,pbc)/sl_width_dist;
                sum_cos2_theta_distprof[slice] += (cosine*cosine);
                count_cos2_theta_distprof[slice]++;
            }

        }
        
        /* Output in classic mode */
        if (bNorm) {
            fprintf(onorm, "%12.7f ", t);
            for (angle=0; angle<n_angles; angle++) {
                fprintf(onorm, "%12.7f ",
                    0.5*(3*sum_cos2_theta_norm[angle]/angle_per_col-1));
                summation += sum_cos2_theta_norm[angle];
            }
            fprintf(onorm, "%12.7f\n",
                0.5*(3*summation/angle_per_row-1));
        }
    } while(read_next_x(oenv,status,&t,natoms,x,box));

    /* Output in profile mode */
    if (bProf) {
        for (slice=0; slice<nslices; slice++) {
            if (count_cos2_theta_prof[slice] == 0) {
                fprintf(oprof, "%12.7f %12.7f\n", slice*sl_width, 0.0);
            }
            else {
                fprintf(oprof, "%12.7f %12.7f\n", slice*sl_width,
                    0.5*(3*sum_cos2_theta_prof[slice]/count_cos2_theta_prof[slice]-1));
                    /* count_cos2_theta[slice], 
                    sum_cos2_theta[slice],
                    sum_cos2_theta[slice]/count_cos2_theta[slice]); */
            }
            if (fnSampling) {
                fprintf(outSampling, "%12.7f %d\n", slice*sl_width,
                        count_cos2_theta_prof[slice]);
            }
        }
    }
    /* Output in distance mode */
    if (bDist) {
        for (slice=0; slice<nslices/2; slice++) {
            if (count_cos2_theta_distprof[slice] == 0) {
                fprintf(odist, "%12.7f %12.7f\n", slice*sl_width_dist, 0.0);
            }
            else {
                fprintf(odist, "%12.7f %12.7f\n", slice*sl_width_dist,
                    0.5*(3*sum_cos2_theta_distprof[slice]/count_cos2_theta_distprof[slice]-1));
            }
            if (fnSamplingDistprof) {
                fprintf(outSamplingDistprof, "%12.7f %d\n", slice*sl_width_dist,
                        count_cos2_theta_distprof[slice]);
            }
        }
    }

    /* Output in grid mode */
    if (bGrid) {
        spacing1_sum /= nframes;
        spacing2_sum /= nframes;
        fprintf(ogrid, "@legend order parameter\n");
        fprintf(outSamplingGrid, "@legend sampling\n");
        fprintf(ogrid, "@xlabel %s (nm)\n", axis_labels[ax1]);
        fprintf(ogrid, "@ylabel %s (nm)\n", axis_labels[ax2]);
        fprintf(ogrid, "@xwidth %f\n", spacing1_sum);
        fprintf(ogrid, "@ywidth %f\n", spacing2_sum);
        fprintf(outSamplingGrid, "@xlabel %s (nm)\n", axis_labels[ax1]);
        fprintf(outSamplingGrid, "@ylabel %s (nm)\n", axis_labels[ax2]);
        fprintf(outSamplingGrid, "@xwidth %f\n", spacing1_sum);
        fprintf(outSamplingGrid, "@ywidth %f\n", spacing2_sum);
        for (i=0; i<nslices; i++) {
            for (j=0; j<nslices2; j++) {
                sum_cos2_theta_grid[i][j] = \
                        order_parameter(sum_cos2_theta_grid[i][j],
                                        count_cos2_theta_grid[i][j]);
                if (j > 0) {
                    fprintf(ogrid, "\t");
                    fprintf(outSamplingGrid, "\t");
                }
                fprintf(ogrid, "%7.3f", sum_cos2_theta_grid[i][j]);
                fprintf(outSamplingGrid, "%.0f", count_cos2_theta_grid[i][j]);
            }
            fprintf(ogrid, "\n");
            fprintf(outSamplingGrid, "\n");
        }
        ffclose(ogrid);
        ffclose(outSamplingGrid);
    }

    /* Cleaning */
    gmx_rmpbc_done(gpbc);
    close_trj(status);
    if (bNorm) {
        ffclose(onorm);
        sfree(sum_cos2_theta_norm);
    }
    if (bProf) {
        ffclose(oprof);
        sfree(sum_cos2_theta_prof);
        sfree(count_cos2_theta_prof);
    }
    if (bGrid) {
        deleteRealMat(sum_cos2_theta_grid, nslices);
        deleteRealMat(count_cos2_theta_grid, nslices);
    }
    if (fnSampling) ffclose(outSampling);
}

/** Read the command line arguments and launch the calculation */
int gmx_cgorder(int argc, char **argv) {
    /* Describe the program */
    const char *desc[] = {
        "g_ordercg calculates the second-rank order parameter for",
        "coarse-grained molecules.",
        "[PAR]",
        "The order parameters is defined as",
        "S = 1/2 * (3 * <cos^2(theta)> - 1)",
        "where theta is the angle between a bond and the membrane normal.",
        "[PAR]",
        "The [TT]-d[tt] option allows to choose the axis normal to the",
        "membrane.",
        "[PAR]",
        "Four different outputs are available: (1) the order parameter for",
        "each bond as a function of time using the [TT]-o[tt] option, (2) the",
        "order parameter profile along an axis with the [TT]-op[tt] option,",
        "(3) the order parameter profile as a function of the distance from a",
        "group of atoms with the [TT]-od[tt] option, and (4) the order",
        "parameter landscape on the membrane plane using the [TT]-og[tt]",
        "option.",
        "[PAR]",
        "By default, the program deduces the list of covalent bonds in the",
        "molecule of interest from the topology file. The [TT]-vectors[tt]",
        "option changes this behaviour: when set to TRUE the index file",
        "is read as a list of vectors, and the order parameter is calculated",
        "for those vectors",
        "(which do not need to represent bonds).",
        "[PAR]",
        "It is also possible to give to the programm a description of the",
        "bonds to consider, using the [TT]-t[tt] option. ",
	"The bonds are described by two atom numbers; numbering is ",
        "relative to the first atom of the molecule. For",
        "instance, the bond described by \"1 2\" is the bond between the first",
        "and the second atoms of the molecule.", 
        "[PAR]",
        "The option [TT]-vectors[tt] is not compatible with the option",
        "[TT]-t[tt]. ",
        "[PAR]",
        "See the README for more information on the software usage."
    };

    int i;

    /* Declare the parameters */
    output_env_t oenv;
    static const char *fnUBL = NULL;
    static const char *fnSampling = NULL;
    static const char *fnSamplingGrid = NULL;
    static const char *fnSamplingDistprof = NULL;
    static const char *axtitle[] = { NULL,"z","x","y",NULL };
    static const char *oaxtitle[] = { NULL,"d","z","x","y",NULL };
    static const char *d2axtitle[] = { NULL,"d","n","z","x","y",NULL };
    int axis = 1;
    int d2axis = 1;
    int oaxis = 1;
    gmx_bool bNorm = FALSE;
    gmx_bool bProf = FALSE;
    gmx_bool bGrid = FALSE;
    gmx_bool bDist = FALSE;
    static int nslices = 100;
    static int nslices2 = -1;
    static gmx_bool vectors=FALSE;
    static gmx_bool bFirst=TRUE;
    static t_pargs pa[] = {
        { "-d", FALSE, etENUM, {axtitle}, 
            "The axis normal to the membrane." },
        { "-sl", FALSE, etINT, {&nslices},
            "Number of bins if used with [TT]-op[tt] ot [TT]-od[tt], number"
            "of cell on the first dimension if used with [TT]-og[tt]."},
        { "-sl2", FALSE, etINT, {&nslices2},
            "Number of cells in the second dimension when [TT]-og[tt] is used."
            "If < 0, the value of [TT]-sl[tt] is used."},
        { "-dp", FALSE, etENUM, {oaxtitle}, 
            "Axis for the order parameter profile if used with [TT]-op[tt]."
            "[TT]d[tt] means default, then the value of [TT]-d[tt] is used."},
        { "-vectors", FALSE, etBOOL, {&vectors}, 
            "Read the index file as a list of vectors." },
        { "-first", FALSE, etBOOL, {&bFirst},
            "Extract the first atom of each molecule. Allow to give a regular"
            "index group instead of a group with only the first atom of"
            "each molecule."},
        { "-2d", FALSE, etENUM, {d2axtitle}, 
            "Switch distance to 2D mode, specify here the axis to remove. "
            "If the [TT]d[tt] value is used then the normal axis will be "
            "the one remove. If [TT]n[tt] is used then distance "
            "calculations are done in 3D." },
    };
    #define NPA asize(pa)

	t_filenm fnm[] = {
        { efTRX, "-f", NULL, ffREAD },
        { efTPX, NULL, NULL, ffREAD },
        { efNDX, NULL, NULL, ffOPTRD },
        { efXVG, "-o", "order", ffOPTWR },
        { efXVG, "-op", "order_profile", ffOPTWR },
        { efDAT, "-og", "order_grid", ffOPTWR },
        { efXVG, "-od", "order_distprof", ffOPTWR },
        { efNDX, "-t", "user_bonds", ffOPTRD  },
        { efXVG, "-osp", "profile_sampling", ffOPTWR },
        { efDAT, "-osg", "grid_sampling", ffOPTWR },
        { efXVG, "-osd", "distprof_sampling", ffOPTWR },
    };
    #define NFILE asize(fnm)
    
    /* Display g_ordercg credit */
    for (i=0; i < (int)asize(version); i++)
        sp_print(stderr, version[i]);
    fprintf(stderr, "\n");
    for (i=0; i < (int)asize(authors); i++)
        sp_print(stderr, authors[i]);
    fprintf(stderr, "\n");
    for (i=0; i < (int)asize(gpl); i++)
        sp_print(stderr, gpl[i]);
    fprintf(stderr, "\n");
    /* Display GROMACS credit */
    CopyRight(stderr,argv[0]);

    /* Parse the command line arguments */
    parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_BE_NICE,
	    NFILE,fnm,NPA,pa,asize(desc),desc,0,NULL,&oenv);
    
    /* Defile what mode to activate */
    bNorm = opt2bSet("-o", NFILE, fnm);
    bProf = opt2bSet("-op", NFILE, fnm);
    bGrid = opt2bSet("-og", NFILE, fnm);
    bDist = opt2bSet("-od", NFILE, fnm);

    /* Check is there is something to do */
    if (!bNorm && !bProf && !bGrid && !bDist) {
        gmx_fatal(FARGS, "You need to choose at least one option among "
                         "-o, -op, -od and -og.");
    }

    /* Check if -t and -vectors are not use together */
    fnUBL = opt2fn_null("-t", NFILE, fnm);
    if (fnUBL != NULL && vectors) {
        gmx_fatal(FARGS, "Options -t and -vectors cannot be used at the same time.\n");
    }
    /* Check if the sampling outputs come with the corresponding mode */
    fnSampling = opt2fn_null("-osp",NFILE,fnm);
    if (fnSampling != NULL && !bProf) {
        gmx_fatal(FARGS, "Option -osp cannot be used "
                         "if -op is not set.\n");
    }
    fnSamplingGrid = opt2fn_null("-osg",NFILE,fnm);
    if (fnSamplingGrid != NULL && !bGrid) {
        gmx_fatal(FARGS, "Option -osg cannot be used "
                         "if -og is not set.\n");
    }
    fnSamplingDistprof = opt2fn_null("-osd",NFILE,fnm);
    if (fnSamplingDistprof != NULL && !bDist) {
        gmx_fatal(FARGS, "Option -osd cannot be used "
                         "if -od is not set.\n");
    }
    /* Check if the number of slices is greater than 0 if we use it */
    if ((bProf || bGrid || bDist) && nslices <= 0) {
        gmx_fatal(FARGS, "Number of slices needs to be greater than 0.\n");
    }

    /* Define axis */
    axis = axtitle[0][0] - 'x';

    /* Adjust slices */
    if (nslices2 <= 0) {
        nslices2 = nslices;
    }

    /* Check -do option */
    /* Is it defined without -sl ? Because it do not make sense. */
    if (oaxtitle[0][0] != 'd' && !bProf) {
        gmx_fatal(FARGS, "Option -dp cannot be used "
			 "if -op is not set.\n");
    }
    /* What is the numerical value for the output axis ? */
    if (oaxtitle[0][0] == 'd') {
        oaxis = axis;
    }
    else {
        oaxis = oaxtitle[0][0] - 'x';
    }

    /* What is the numerical value for the 2D mode axis ? */
    if (d2axtitle[0][0] == 'n') {
        d2axis = -1;
    }
    else if (d2axtitle[0][0] == 'd') {
        d2axis = axis;
    }
    else {
        d2axis = d2axtitle[0][0] - 'x';
    }

    /* Do the job */
    do_cgorder(ftp2fn(efTRX,NFILE,fnm), ftp2fn(efTPX,NFILE,fnm),
        ftp2fn(efNDX,NFILE,fnm),
        opt2fn("-o",NFILE,fnm), opt2fn("-op",NFILE,fnm),
        opt2fn("-og",NFILE,fnm), opt2fn("-od",NFILE,fnm),
        opt2fn("-osp",NFILE,fnm), opt2fn("-osg",NFILE,fnm),
        opt2fn("-osd",NFILE,fnm), fnUBL,
        oenv, axis, oaxis, nslices, nslices2,
        vectors, bFirst, bNorm, bProf, bGrid, bDist, d2axis);

    thanx(stderr);
    return 0;
}

/** Entry point of the program */
int main(int argc, char **argv) {
    return gmx_cgorder(argc, argv);
}


