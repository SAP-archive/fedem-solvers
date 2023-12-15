// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFaLib/FFaCmdLineArg/FFaCmdLineArg.H"
#include "FFaLib/FFaOS/FFaFortran.H"
#include "Admin/FedemAdmin.H"
#include <cstring>
#include <cstdio>


#define ADDOPTION FFaCmdLineArg::instance()->addOption

SUBROUTINE(setversion,SETVERSION) (const char* theVer, const int nchar);
SUBROUTINE(setdate,SETDATE) (const char* theDate, const int nchar);
SUBROUTINE(response_pos,RESPONSE_POS) (int& iret);


const size_t bufsiz = 256;


static char* readLine (char* cline, FILE* fp)
{
  return fgets(cline,bufsiz,fp);
}

static size_t writeLine (char* cline, FILE* fp)
{
  return fwrite(cline,sizeof(char),strlen(cline),fp);
}


int main (int argc, char** argv)
{
  // Initialize the command-line parser
  FFaCmdLineArg::init (argc,argv);

  // Initialize only the options that are actually evaluated from this program
  ADDOPTION ("consolemsg",true,"Output error messages to console");
  ADDOPTION ("rdbinc",1,"Increment number for the results database file");
  ADDOPTION ("errorfile","","Name of error message file");
  ADDOPTION ("linkfile","","Name of link input file");
  ADDOPTION ("deformation",false,"Save deformations to results database");
  ADDOPTION ("profile",false,"Print out profiling data at termination",false);

  // Get program version and build date
  const char* fedem_version = FedemAdmin::getVersion();
  const char* build_date = FedemAdmin::getBuildDate();

  // Pass the version tag and date to the fortran modules
  F90_NAME(setversion,SETVERSION) (fedem_version,strlen(fedem_version));
  F90_NAME(setdate,SETDATE) (build_date,strlen(build_date));

  // Lauch the Fortran program
  int iret = 0;
  F90_NAME(response_pos,RESPONSE_POS) (iret);
  FFaCmdLineArg::removeInstance ();
  if (iret) return iret;

  // Merge the two frs-files into a single file response_pos.frs

  FILE* fOut = fopen("response_pos.frs","w");
  FILE* fIn1 = fopen("response_pos_1.frs","r");
  FILE* fIn2 = fopen("gage_pos_2.frs","r");
  if (!fOut || !fIn1 || !fIn2) return -1;

  char cline[bufsiz], dline[bufsiz];

  // First, copy lines from the first file until the first '['-character
  printf("Merging variable definitions ...\n");
  while (readLine(cline,fIn1) && cline[0] != '[')
    writeLine(cline,fOut);

  // Then, insert all lines starting with '<' from second file
  int lcount = 0;
  while (readLine(dline,fIn2) && dline[0] != '[')
    if (dline[0] == '<' && ++lcount > 2)
      writeLine(dline,fOut);

  // Next, copy all lines from first file until 'DATABLOCKS:'
  printf("Merging item group definitions ...\n");
  do
    writeLine(cline,fOut);
  while (readLine(cline,fIn1) && strncmp(cline,"DATABLOCKS:",11));

  // Then, insert all lines from second file until 'DATABLOCKS:'
  do
    writeLine(dline,fOut);
  while (readLine(dline,fIn2) && strncmp(dline,"DATABLOCKS:",11));

  // Next, copy all lines from first file until 'DATA:'
  printf("Merging object group definitions ...\n");
  do
    writeLine(cline,fOut);
  while (readLine(cline,fIn1) && strncmp(cline,"DATA:",5));

  // Then, insert all lines starting with '{' from second file
  while (readLine(dline,fIn2) && dline[0] != 'D')
    if (dline[0] == '{')
      writeLine(dline,fOut);

  // Finally, write the rest of the first file
  printf("Done.\n");
  do
    writeLine(cline,fOut);
  while (readLine(cline,fIn1));
  fprintf(fOut,"\n");

  // Close files
  fclose(fOut);
  fclose(fIn1);
  fclose(fIn2);

  return 0;
}
