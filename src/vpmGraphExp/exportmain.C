// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFpLib/FFpExport/FFpBatchExport.H"
#include "FFaLib/FFaString/FFaTokenizer.H"
#include "FFaLib/FFaCmdLineArg/FFaCmdLineArg.H"
#include "FFaLib/FFaOS/FFaFilePath.H"
#include "Admin/FedemAdmin.H"

extern void cmdLineArgInitStd (int argc, char** argv);
extern void readOptionFilesStd (const char* program);

#define ADDOPTION FFaCmdLineArg::instance()->addOption
#define GETOPTION FFaCmdLineArg::instance()->getValue


int main (int argc, char** argv)
{
  // Print out program name, version and build date
  std::cout <<"\nFedem Curve Export utility "
            << FedemAdmin::getVersion() <<" "<< FedemAdmin::getBuildDate();
  if (!FedemAdmin::is64bit()) std::cout <<"   (32bit)";
  std::cout << std::endl;

  // Lambda function that prints an error message and closes down.
  auto&& finished = [](const char* msg = NULL)
  {
    if (msg) std::cout <<"\n===> "<< msg <<".\n";
    FFaCmdLineArg::removeInstance();
    return msg ? 1 : 0;
  };

  // Initialize the command-line parser
  cmdLineArgInitStd (argc,argv);

  // Define all command-line options for this program
  ADDOPTION ("frsFile","","List of results database files");
  ADDOPTION ("rpcFile","","Get number of repeats, averages, and"
             "\npoints per frame and group, from this RPC-file");
  ADDOPTION ("modelFile","","Name of model file with curve definitions");
  ADDOPTION ("curvePlotFile","response.rsp","Name of curve export output file");
  ADDOPTION ("curvePlotType",3,"Format of curve export output file"
             "\n= 0 : ASCII (separate file for each curve)"
             "\n= 1 : DAC, Windows (separate file for each curve)"
             "\n= 2 : DAC, UNIX (separate file for each curve)"
             "\n= 3 : RPC, Windows (all curves in one file)"
             "\n= 4 : RPC, UNIX (all curves in one file)"
             "\n= 5 : ASCII (all curves in one file)"
             "\n= 9 : Only print file position data for each quantity");
  ADDOPTION ("curvePlotPrec",0,"Output precision for curve data files"
             "\n= 0 : half precision (int*2)"
             "\n= 1 : single precision (real*4)"
             "\n= 2 : double precision (real*8)");

  // Read the option files if any
  readOptionFilesStd ("fedem_graphexp");

  // Evaluate the command-line options
  int format, prec;
  std::string frsFile, rpcFile, modelFile, plotFile;
  GETOPTION ("frsFile",frsFile);
  GETOPTION ("rpcFile",rpcFile);
  GETOPTION ("modelFile",modelFile);
  GETOPTION ("curvePlotFile",plotFile);
  GETOPTION ("curvePlotType",format);
  GETOPTION ("curvePlotPrec",prec);

  if (modelFile.empty())
    return finished("No curve definition file (model file) given");
  else
    FFaFilePath::checkName(modelFile);

  FFaTokenizer frsFiles(frsFile,'<','>',',');
  if (frsFiles.empty() && !FFpBatchExport::readFrsFiles(frsFiles,modelFile))
    return finished("No frs-files given");

  std::cout <<"\n===> Reading results data from the following RDB-files:\n";
  for (std::string& fileName : frsFiles)
  {
    FFaFilePath::checkName(fileName);
    std::cout <<"     "<< fileName << std::endl;
  }

  // Read number of repeats, etc., from the specified input RPC-file
  FFpRPC3Data rpc;
  if (!rpcFile.empty() && format > 2 && format < 5)
    if (!rpc.readDataFromFile(FFaFilePath::checkName(rpcFile)))
      return finished("Exporting Curves failed");

  // Open frs-files and read the curve definitions
  std::cout <<"\n===> Reading curve definitions from "<< modelFile << std::endl;
  FFpBatchExport exporter(frsFiles);
  bool success = exporter.readCurves(modelFile);
  if (!success)
    return finished("Exporting Curves failed");

  if (format == 9)
  {
    // No export, just print out positioning data
    if (plotFile.find(".rsp") != std::string::npos)
      success = exporter.printPosition();
    else
      success = exporter.printPosition(FFaFilePath::checkName(plotFile));
  }
  else
  {
    // Now export the curves
    std::cout <<"\n===> Exporting Curves to "
              << FFaFilePath::checkName(plotFile) << std::endl;
    if (format > 2)
      success = exporter.exportGraph(plotFile,modelFile,prec*10+format%5,rpc);
    else
      success = exporter.exportCurves(plotFile,modelFile,prec*10+format);
    std::cout <<"===> Exporting Curves "<< (success ? "done":"failed") <<".\n";
  }

  return finished(success ? NULL : "Failure");
}
