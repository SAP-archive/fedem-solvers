// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include <cstring>
#include <vector>
#include <map>

#include "FFaLib/FFaCmdLineArg/FFaCmdLineArg.H"
#include "FFaLib/FFaOS/FFaFortran.H"

#include "base/base.h"
#include "vdm/vdm.h"

typedef std::map< Vint, std::vector<Vfloat> > StressMap;

static StressMap strmap;


static int loadStress (const char*  setName,
		       vdm_DataFun* datafun,
		       vdm_Library* library,
		       const Vint&  numel,
		       const Vint*  ptr_eid)
{
  strmap.clear();

  // search for stress results datasets
  Vint ndst_s, idst_s[10];
  vdm_LibrarySearchDataset (library,"S.EL:*",10,idst_s,&ndst_s);
  if (ndst_s == 0) return 3;
  if (ndst_s > 10) ndst_s = 10;

  // query first dataset for sizes to allocate memory
  vdm_Dataset* dataset;
  vdm_LibraryGetDataset (library,idst_s[0],&dataset);
  Vchar dsname[65];
  Vint lrec, nrow, ncol, ntyp;
  vdm_DatasetInq (dataset,dsname,&lrec,&nrow,&ncol,&ntyp);
  Vfloat* ptr_s = (Vfloat*)malloc(lrec*sizeof(Vfloat));

  // find linked datasets for pointers and sizes
  Vint iatsize, nat;
  vdm_DatasetSearchAttribute (dataset,"Link.Size",1,&iatsize,&nat);
  vdm_Attribute* attribute;
  vdm_DatasetGetAttribute (dataset,iatsize,&attribute);
  Vchar cvalue[65];
  vdm_AttributeGet (attribute,cvalue);
  Vint idssize, nds;
  vdm_LibrarySearchDataset (library,cvalue,1,&idssize,&nds);
  Vint* ptr_size = (Vint*)malloc(ncol*sizeof(Vint));
  vdm_DataFunReadDataset (datafun,idssize,ptr_size);

  Vint iatpntr, idspntr;
  vdm_DatasetSearchAttribute (dataset,"Link.Pointer",1,&iatpntr,&nat);
  vdm_DatasetGetAttribute (dataset,iatpntr,&attribute);
  vdm_AttributeGet (attribute,cvalue);
  vdm_LibrarySearchDataset (library,cvalue,1,&idspntr,&nds);
  Vint* ptr_pntr = (Vint*)malloc(ncol*sizeof(Vint));
  vdm_DataFunReadDataset (datafun,idspntr,ptr_pntr);

  // loop over stress datasets
  Vint i, j, k, n, elm, nix;
  StressMap::iterator sit;
  for (i = 0; i < ndst_s; i++)
  {
    vdm_LibraryGetDataset (library,idst_s[i],&dataset);
    vdm_DatasetInq (dataset,dsname,&lrec,&nrow,&ncol,&ntyp);
    printf("\nFound stress dataset: %s",dsname);

    // if file contains more than one data set, check that we have the right one
    if (setName && ndst_s > 1)
      if (strcmp(setName,dsname)) continue;

    // read entire dataset and print select elements
    printf("  reading ...");
    vdm_DataFunReadDataset (datafun,idst_s[i],ptr_s);

    // loop over all elements
    for (n = 0; n < numel; n++)
    {
      // compute number of nodes in element from size
      // and assign a stress vector in a map
      elm = ptr_eid[n];
      nix = ptr_size[n]/nrow;
      sit = strmap.insert(make_pair(elm,std::vector<Vfloat>(6*nix,0.0))).first;

      // store nodal stress values in the map
      for (j = 0; j < nix; j++)
	for (k = 0; k < 6; k++)
	  sit->second[k+6*j] = ptr_s[k+nrow*j+ptr_pntr[n]-1];
    }
    break;
  }
  printf("\n\n");

  // free memory
  free(ptr_s);
  free(ptr_size);
  free(ptr_pntr);

  if (strmap.empty())
    if (setName)
      fprintf(stderr,"Dataset %s not found\n",setName);
    else
      fprintf(stderr,"No dataset found\n");

  return strmap.empty() ? 4 : 0;
}


int readStress (const char* fileName, const char* setName)
{
  if (!fileName)
  {
    fprintf(stderr,"readStress: no input file specified!\n");
    return 1;
  }

  // determine file type from file extension
  Vint filetype = 0;
  if (strstr(fileName,".unv"))
    filetype = VDM_SDRC_UNIVERSAL;
  else if (strstr(fileName,".fil"))
    filetype = VDM_ABAQUS_FIL;
  else if (strstr(fileName,".rst") || strstr(fileName,".rth"))
    filetype = VDM_ANSYS_RESULT;
  else if (strstr(fileName,".op2"))
    filetype = VDM_NASTRAN_OUTPUT2;
  else
  {
    fprintf(stderr,"readStress: Bad input file type %s\n",fileName);
    return 2;
  }

  // create data function and file reader objects
  vdm_DataFun* datafun = vdm_DataFunBegin();
  vdm_SDRCLib* sdrclib = 0;
  vdm_ABALib*  abalib  = 0;
  vdm_ANSLib*  anslib  = 0;
  vdm_NASLib*  naslib  = 0;
  switch (filetype)
    {
    case VDM_SDRC_UNIVERSAL:
      sdrclib = vdm_SDRCLibBegin();
      vdm_SDRCLibDataFun (sdrclib,datafun);
      break;
    case VDM_ABAQUS_FIL:
      abalib = vdm_ABALibBegin();
      vdm_ABALibDataFun (abalib,datafun);
      break;
    case VDM_ANSYS_RESULT:
      anslib = vdm_ANSLibBegin();
      vdm_ANSLibDataFun (anslib,datafun);
      break;
    case VDM_NASTRAN_OUTPUT2:
      naslib = vdm_NASLibBegin();
      vdm_NASLibDataFun (naslib,datafun);
      break;
    }

  // open library device
  printf("\nReading external result file: %s\n",fileName);
  vdm_DataFunOpen (datafun,"readStress",(Vchar*)fileName,filetype);

  // get number of elements
  Vint numel;
  vdm_DataFunGetNumEntities (datafun,SYS_ELEM,&numel);

  // get library object
  vdm_Library* library = NULL;
  vdm_DataFunGetLibrary (datafun,&library);

  // search for element number data sets
  Vint idst_eid, ndst_eid;
  vdm_LibrarySearchDataset (library,"EID.E",1,&idst_eid,&ndst_eid);

  // read user element numbers if they exist
  Vint* ptr_eid = (Vint*)malloc(numel*sizeof(Vint));
  if (ndst_eid)
    vdm_DataFunReadDataset (datafun,idst_eid,ptr_eid);
  else
    for (Vint i = 0; i < numel; i++)
      ptr_eid[i] = i+1;

  // load stresses into a global array
  int status = loadStress(setName,datafun,library,numel,ptr_eid);

  // free memory
  free(ptr_eid);

  // close library device
  vdm_DataFunClose (datafun);

  // free objects
  vdm_DataFunEnd (datafun);
  switch (filetype)
    {
    case VDM_SDRC_UNIVERSAL:
      vdm_SDRCLibEnd (sdrclib);
      break;
    case VDM_ABAQUS_FIL:
      vdm_ABALibEnd (abalib);
      break;
    case VDM_ANSYS_RESULT:
      vdm_ANSLibEnd (anslib);
      break;
    case VDM_NASTRAN_OUTPUT2:
      vdm_NASLibEnd (naslib);
      break;
    }

  return status;
}


// Fortran-callable functions

SUBROUTINE(readresstress,READRESSTRESS) (int& ierr)
{
  // Read residual stress file, if specified
  std::string fileName, setName;
  FFaCmdLineArg::instance()->getValue ("resStressFile",fileName);
  FFaCmdLineArg::instance()->getValue ("resStressSet",setName);
  if (fileName.empty())
    ierr = 1;
  else
    ierr = -readStress(fileName.c_str(),setName.c_str());
}


SUBROUTINE(getsolidstress,GETSOLIDSTRESS) (const int& elmNo, const int* mnpc,
					   double* Sigma, int& nstrp, int& ierr)
{
  StressMap::const_iterator sit = strmap.find(elmNo);
  if (sit == strmap.end())
  {
    ierr = -1;
    fprintf(stderr,"getSolidStress: No stresses for element %d\n",elmNo);
    return;
  }

  int nstres = sit->second.size();
  if (nstrp == 0)
    nstrp = nstres/6;
  else if (nstrp*6 != nstres)
  {
    ierr = -2;
    fprintf(stderr,"getSolidStress: Error for solid element %d\n",elmNo);
    fprintf(stderr,"                nstrp*6 != nstres %d %d\n",nstrp*6,nstres);
    return;
  }

  // Assume Nastran ordering for the second order solid elements

  const int nodePermT10[10] = { 1, 3, 5,10, 2, 4, 6, 7, 8, 9 };

  const int nodePermW15[15] = { 1, 3, 5,10,12,14, 2, 4, 6, 7,
				8, 9,11,13,15 };

  const int nodePermH20[20] = { 1, 3, 5, 7,13,15,17,19, 2, 4,
				6, 8, 9,10,11,12,14,16,18,20 };

  ierr = 0;
  int i, j, n, k = 0;
  for (i = 0; i < nstrp; i++)
  {
    if (nstrp == 10)
      n = 6*(nodePermT10[i]-1);
    else if (nstrp == 15)
      n = 6*(nodePermW15[i]-1);
    else if (nstrp == 20)
      n = 6*(nodePermH20[i]-1);
    else
      n = 6*(i-1);

    for (j = 0; j < 6; j++)
      Sigma[n+j] += sit->second[k++];
  }
}


SUBROUTINE(getshellstress,GETSHELLSTRESS) (const int& elmNo, const int* mnpc,
					   double* Sigma, int& nstrp, int& ierr)
{
  ierr = -1;
  fprintf(stderr,"getShellStress: Not implemented yet!\n");
}
