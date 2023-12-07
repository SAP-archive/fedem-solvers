<!---
  SPDX-FileCopyrightText: 2023 SAP SE

  SPDX-License-Identifier: Apache-2.0

  This file is part of FEDEM - https://openfedem.org
--->

# Step-by-Step guide for sub-model analysis in FEDEM R7.3.2 and later

This file contains a detailed description on how to set up
and manually execute a sub-model analysis in FEDEM, using the
desktop edition on Windows. It's main purpose is to describe
the actual process required, and can be used to develop an
automated process for cloud execution later.

## Restrictions

* Only quasi-static simulation without distributed loads
* No gravity (or neutral buoyancy)

## Pre-requisites

The sub-model analysis in FEDEM assumes that you have a relatively coarse
global model, which can consist of one or more superelements with shell
elements (currently only 3- and 4-noded shell elements are supported,
no solids or beams in the area of the sub-models), and a bunch of
finer sub-models each one covering an area of special interrest.

From each sub-model, extract a list of coordinates of all the nodal points
in which you want to input displacement boundary conditions from the
global model. This list should be stored in a file `SubModelNodes.dat`,
consisting of 5 comma-separated columns of data.

- Column #1 : Interface ID (not used)
- Column #2 : Node ID (not used)
- Column #3 : X-coordinate
- Column #4 : Y-coordinate
- Column #5 : Z-coordinate

### Alternative node file format

A simpler version of the node file is to just list the node numbers as a
single-column file, alternatively preceeded by the interface ID as a two-column file.
The coordinates of the interface nodes are then taken from the sub-model FE data
file instead (specified using the `-subFile` option, see below), which can be on
any of the supported FE data file formats (`.nas`, `.ftl`, `.FEM`, etc.).

## How to do it

1. Create a FEDEM model of the global model by importing the FE part(s)
with the neccessary attributes, loads, etc. Remember to switch off the
gravity vector in the *Model Preferences*. And toggle for "Quasistatic analysis"
and "Complete interval" in the *Dynamics Solver Setup* dialog.

2. (Optional) Enter `-nograv` in the "FE model Reducer" field of the
*Additional Solver Options* dialog. And set "Component modes" to 0 in the
"Reduction Options" tab of the property panel for all FE parts in the model.
This will save computation time in the model reduction since the mass matrix
and gravity force vector calculation is not needed in this simulation.

3. Run the dynamics solver on the global model.

4. Run stress recovery on the FE part(s). In the *Stress Recovery Setup*,
   the "deformation" toggle needs to be on, no others are required, but it is
   maybe OK to check the "von Mises stress" also for verification against the
   subsequent sub-model results.

5. If all goes well, save the model (`GlobalModel.fmm`) and exit the FEDEM GUI.

6. Open a console window, and cd to the folder in which you have the models.

7. Run the `fedem_solmap` utility on each sub-model (`SubModel.nas`):

        C:\<path-to-fedem-installation>\fedem_solmap -fco submodel.fco
   where `submodel.fco` is a file containing all the neccessary input options.
   It typically contains the following:

        -ftlFile GlobalModel_RDB\link_DB\GlobalPart.ftl
        -frsFile GlobalModel_RDB\response_0001\timehist_rcy\1_GlobalPart_0001\GlobalPart_1.frs
        -nodeFile SubModelNodes.dat
        -subFile SubModel.nas
        -outFile SubModel.fnd
        -translate 123.0 0.0 0.0
        -rotate 0.0 0.0 60.0
        -group 1
   The `-translate` and `-rotate` options only need to be specified if
   the sub-model does not have its origin in the origin of the global model,
   and/or the axis directions are not aligned.
   The `-group` option is normally not needed.
   It can be used to specify the ID of an element group to limit the search to.
   The program creates two output files:
   - `SubModel.ftl` will consist of the imported `SubModel.nas` file,
      but with each nodal point that was identified in the `SubModelNodes.dat`
      marked as external node.
   - `SubModel.fnd` is a binary file containing the time history of nodal
      translations of all the nodes listed in `SubModelNodes.dat`

8. Start the FEDEM GUI again with a blank model.

9. Import the sub-model, i.e., the generated `SubModel.ftl` -- not `SubModel.nas`.
   All nodes on the interfaces are then already marked as Triads.

10. In the Objects browser, multi-select all triads, and change the "Constraint Type" to
    "Prescribed" for the "Tx", "Ty" and "Tz" dofs.

11. In the Property panel of the FE part SubModel, add `#Displaced` in the description field.

12. Set number of "Component modes" to 0.

13. If the sub-model is not modelled in an axis system that is alligend with the
    global model axes, you need to rotate it such that the axes do align.
    Otherwise, the displacements from the global model will not be applied correctly.
    There is no need to translate the sub-model.

14. In the *Model Preferences*, set the "Gravitation" vector to 0.0 0.0 0.0

15. In the *Additional Solver Options* dialog
    - Enter `-nograv` in the "FE model Reducer" field
    - Enter `-displacementfile ..\..\SubModel.fnd` in the "Dynamics Solver" field

16. In the *Dynamics Solver Setup* dialog
    - Toggle "Quasistatic analysis" and "Complete interval"
    - Toggle "Ignore integration tolerances" and set "Number of iterations" to -1
    - Toggle off "Geometric stiffness contribution" and "Centripetal force correction"

17. Run the dynamics solver on the sub-model.

18. Set up and run stress or strain gage recovery on the sub-model to extract
    the quantity of interest. This can also be done as part of the dynamics solver
    step, if wanted, by using the `#recover-stress` or `#recover-gage` options.

### Running global model recovery during the system simulation

The steps 3 and 4 above may be combined into one operation, by doing the recovery
during the time step analysis itself, using the `#recover-stress` option.
Since the purpose of this recovery only is to obtain the displacement state at
the sub-model interfaces, the most efficient way to do this is to define a set of
element groups in the global model, containing only those elements that cover the
interface area. It is then possible to perform displacement recovery on those
elements only, which will save both computation time and disk space usage.

1. Let's say you have defined the element groups with ID 11, 12 and 13 for this purpose.
Then in the description field of the FE part containing those element groups
   - Enter the string `#recover-stress <11,12,13>`.
     If only one group, the brackets "<>" can be skipped.

2. In the *Additional Solver Options* dialog
   - Enter `-partDeformation=2 -partVMStress=0` in the "Dynamics Solver" field.

3. Run the dynamics solver.
   The `-partVMStress=0` option will switch off all stress calculations
   during the simulation, such that only the total displacements are calculated.
   This will speed up the overall simulation time considerably, as the stress
   calculation is the task that takes the longest time.

### Running the solution mapping in two stages

If your simulation consists of a series of micro-batches on the same global model,
the solution mapping process can be splitted into two separate stages:

1. Search for matching nodes or elements for the interface nodes.
2. Perform the displacement interpolation from global to sub-model.

Where you typically need to perform the first stage (which often is the most
time consuming) only once, whereas the second stage is performed for each batch.
This can be achieved by storing the search results in an intermediate file,
using the command-line option `-mapFile`. Thus, you can have two setup files
for the `fedem_solmap` utility, one for each of the two stages, as follows:

1. Searching stage, `submodel1.fco` with contents:

        -ftlFile GlobalModel_RDB\link_DB\GlobalPart.ftl
        -nodeFile SubModelNodes.dat
        -subFile SubModel.nas
        -mapFile SubModel.map
        -translate 123.0 0.0 0.0
        -rotate 0.0 0.0 60.0
        -group 1

2. Interpolation stage, `submodel2.fco` with contents:

        -frsFile GlobalModel_RDB\response_0001\timehist_rcy\1_GlobalPart_0001\GlobalPart_1.frs
        -mapFile SubModel.map
        -outFile SubModel.fnd

This should be all. Good luck!
