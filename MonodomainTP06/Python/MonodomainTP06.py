#!/usr/bin/env python

#> \file
#> \author Sander Arens
#> \brief This is an example script to solve a Monodomain problem with the TP06 model using openCMISS calls in python.
#>
#> \section LICENSE
#>
#> Version: MPL 1.1/GPL 2.0/LGPL 2.1
#>
#> The contents of this file are subject to the Mozilla Public License
#> Version 1.1 (the "License"); you may not use this file except in
#> compliance with the License. You may obtain a copy of the License at
#> http://www.mozilla.org/MPL/
#>
#> Software distributed under the License is distributed on an "AS IS"
#> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
#> License for the specific language governing rights and limitations
#> under the License.
#>
#> The Original Code is openCMISS
#>
#> The Initial Developer of the Original Code is University of Auckland,
#> Auckland, New Zealand and University of Oxford, Oxford, United
#> Kingdom. Portions created by the University of Auckland and University
#> of Oxford are Copyright (C) 2007 by the University of Auckland and
#> the University of Oxford. All Rights Reserved.
#>
#>
#> Alternatively, the contents of this file may be used under the terms of
#> either the GNU General Public License Version 2 or later (the "GPL"), or
#> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
#> in which case the provisions of the GPL or the LGPL are applicable instead
#> of those above. if you wish to allow use of your version of this file only
#> under the terms of either the GPL or the LGPL, and not to allow others to
#> use your version of this file under the terms of the MPL, indicate your
#> decision by deleting the provisions above and replace them with the notice
#> and other provisions required by the GPL or the LGPL. if you do not delete
#> the provisions above, a recipient may use your version of this file under
#> the terms of any one of the MPL, the GPL or the LGPL.
#>

#> \example ClassicalField/Diffusion/Diffusion/Python/DiffusionExample.py
## Example script to solve a Diffusion problem using openCMISS calls in python.
#<


# Add Python bindings directory to PATH
import sys, os
sys.path.append(os.path.join((os.environ['OPENCMISS_ROOT'], 'cm', 'bindings', 'python')))

# Intialise OpenCMISS
from opencmiss import CMISS

# Set problem parameters
height = 1.0
width = 1.0
length = 0.05

numberGlobalXElements = 20
numberGlobalYElements = 20
numberGlobalZElements = 1 

numberOfXi = 3

# Set time integration parameters
startTime = 0.0
stopTime = 300.01
timeIncrement = 0.01

# Set point stimulus
activ_r = 3.0e-02 

# Set conduction parameters
FiberD = 0.095298372513562 
TransverseD = 0.0125758411473

#Set user numbers
coordinateSystemUserNumber = 1
regionUserNumber = 1
basisUserNumber = 1
generatedMeshUserNumber = 1
meshUserNumber = 1
decompositionUserNumber = 1
geometricFieldUserNumber = 1
equationsSetFieldUserNumber = 2
equationsSetUserNumber = 1
dependentFieldUserNumber = 3
materialsFieldUserNumber = 4
independentFieldUserNumber = 5
problemUserNumber = 1

#Set output frequency for nodes
outputFrequency = int(1/timeIncrement) 
#outputFrequency = 1 

CMISS.DiagnosticsSetOn(CMISS.DiagnosticTypes.IN, [1,2,3,4,5], "Diagnostics", ["DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE"])
CMISS.ErrorHandlingModeSet(CMISS.ErrorHandlingModes.TRAP_ERROR)

# Get the computational nodes information
numberOfComputationalNodes = CMISS.ComputationalNumberOfNodesGet()
computationalNodeNumber = CMISS.ComputationalNodeNumberGet()

# Creation a RC coordinate system
coordinateSystem = CMISS.CoordinateSystem()
coordinateSystem.CreateStart(coordinateSystemUserNumber)
if numberGlobalZElements == 0:
	coordinateSystem.DimensionSet(2)
else: 
	coordinateSystem.DimensionSet(3)
coordinateSystem.CreateFinish()

# Create a region
region = CMISS.Region()
region.CreateStart(regionUserNumber, CMISS.WorldRegion)
region.LabelSet("MonodomainRegion")
region.CoordinateSystemSet(coordinateSystem)
region.CreateFinish()

# Create a tri-linear lagrange basis
basis = CMISS.Basis()
basis.CreateStart(basisUserNumber)
basis.TypeSet(CMISS.BasisTypes.LAGRANGE_HERMITE_TP) #TP stands for Tensor Product
if numberGlobalZElements == 0:
	basis.NumberOfXiSet(2)
	basis.InterpolationXiSet([CMISS.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*2)
	basis.QuadratureNumberOfGaussXiSet([2]*2)
else: 
	basis.NumberOfXiSet(3)
	basis.InterpolationXiSet([CMISS.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*3)
	basis.QuadratureNumberOfGaussXiSet([2]*3)
basis.CreateFinish()

# Create a generated mesh
generatedMesh = CMISS.GeneratedMesh()
generatedMesh.CreateStart(generatedMeshUserNumber, region)
generatedMesh.TypeSet(CMISS.GeneratedMeshTypes.REGULAR)
generatedMesh.BasisSet([basis])
if numberGlobalZElements == 0:
	generatedMesh.ExtentSet([width,height])
	generatedMesh.NumberOfElementsSet([numberGlobalXElements, numberGlobalYElements])
else: 
	generatedMesh.ExtentSet([width,height,length])
	generatedMesh.NumberOfElementsSet([numberGlobalXElements, numberGlobalYElements, numberGlobalZElements])
mesh = CMISS.Mesh()
generatedMesh.CreateFinish(meshUserNumber, mesh)

# Create a decomposition for the mesh
decomposition = CMISS.Decomposition()
decomposition.CreateStart(decompositionUserNumber, mesh)
decomposition.TypeSet(CMISS.DecompositionTypes.CALCULATED)
decomposition.NumberOfDomainsSet(numberOfComputationalNodes)
decomposition.CreateFinish()

# Create a field for the geometry
geometricField = CMISS.Field()
geometricField.CreateStart(geometricFieldUserNumber, region)
geometricField.MeshDecompositionSet(decomposition)
geometricField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U, 1, 1)
geometricField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U, 2, 1)
if numberGlobalZElements > 0:
	geometricField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U, 3, 1)
geometricField.CreateFinish()

# Set geometry from the generated mesh
generatedMesh.GeometricParametersCalculate(geometricField)

# Create monodomain equations set
equationsSetField = CMISS.Field()
equationsSet = CMISS.EquationsSet()
equationsSet.CreateStart(equationsSetUserNumber, region, geometricField,
        CMISS.EquationsSetClasses.BIOELECTRICS,
        CMISS.EquationsSetTypes.MONODOMAIN_STRANG_SPLITTING_EQUATION,
        CMISS.EquationsSetSubtypes.MONODOMAIN_TENTUSSCHER06,
        equationsSetFieldUserNumber, equationsSetField)
equationsSet.CreateFinish()

# Create dependent field
dependentField = CMISS.Field()
equationsSet.DependentCreateStart(dependentFieldUserNumber, dependentField)
dependentField.VariableLabelSet(CMISS.FieldVariableTypes.U, "Vm")
dependentField.VariableLabelSet(CMISS.FieldVariableTypes.DELUDELN, "del Vm/del n")
equationsSet.DependentCreateFinish()

# Create the equations set material field variables
materialsField = CMISS.Field()
equationsSet.MaterialsCreateStart(materialsFieldUserNumber, materialsField)
materialsField.VariableLabelSet(CMISS.FieldVariableTypes.U, "Material constants")
equationsSet.MaterialsCreateFinish()
nodes = CMISS.Nodes()
region.NodesGet(nodes)
for node in xrange(1, nodes.numberOfNodes + 1):
    nodeDomain = decomposition.NodeDomainGet(node, 1)
    if nodeDomain == computationalNodeNumber:
        position = [
        geometricField.ParameterSetGetNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, node, component + 1)
        for component in range(numberOfXi)]
    else:
        position.append(0.0)
    if abs(position[0]-height*0.5) <= activ_r and abs(position[1]-width*0.5) < activ_r and abs(position[2]-length*0.5) < activ_r:
        activ = 100.0
        print node
    else:
        activ = 0.0
    materialsField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, node, 1, activ) 
materialsField.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 2, FiberD) 
materialsField.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 3, TransverseD)
materialsField.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 4, 1.0) 
materialsField.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 5, 1.0)
if numberGlobalZElements > 0:
	materialsField.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 6, 0.0)
materialsField.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 7, 0.0)

# Create/init independent field
independentField = CMISS.Field()
equationsSet.IndependentCreateStart(independentFieldUserNumber, independentField)
independentField.VariableLabelSet(CMISS.FieldVariableTypes.U, "State variables")
equationsSet.IndependentCreateFinish()

# Create equations
equations = CMISS.Equations()
equationsSet.EquationsCreateStart(equations)
equations.SparsityTypeSet(CMISS.EquationsSparsityTypes.SPARSE)
equations.OutputTypeSet(CMISS.EquationsOutputTypes.NONE)
equationsSet.EquationsCreateFinish()

# Create monodomain problem
problem = CMISS.Problem()
problem.CreateStart(problemUserNumber)
problem.SpecificationSet(CMISS.ProblemClasses.BIOELECTRICS,
        CMISS.ProblemTypes.MONODOMAIN_STRANG_SPLITTING_EQUATION,
        CMISS.ProblemSubTypes.MONODOMAIN_TENTUSSCHER06)
problem.CreateFinish()

# Create control loops
controlLoop = CMISS.ControlLoop()
problem.ControlLoopCreateStart()
problem.ControlLoopGet([CMISS.ControlLoopIdentifiers.NODE], controlLoop)
controlLoop.OutputTypeSet(CMISS.ControlLoopOutputTypes.PROGRESS)
controlLoop.TimesSet(startTime, stopTime, timeIncrement)
controlLoop.TimeOutputSet(outputFrequency)
problem.ControlLoopCreateFinish()

# Create problem solver
solver = CMISS.Solver()
linearSolver = CMISS.Solver()
problem.SolversCreateStart()
problem.SolverGet([CMISS.ControlLoopIdentifiers.NODE], 1, solver)
solver.OutputTypeSet(CMISS.SolverOutputTypes.NONE)
solver.DynamicLinearSolverGet(linearSolver)
solver.DynamicSchemeSet(CMISS.DynamicSchemeTypes.CRANK_NICOLSON)
linearSolver.LinearIterativeMaximumIterationsSet(100)
problem.SolversCreateFinish()

# Create solver equations and add equations set to solver equations
solver = CMISS.Solver()
solverEquations = CMISS.SolverEquations()
problem.SolverEquationsCreateStart()
problem.SolverGet([CMISS.ControlLoopIdentifiers.NODE], 1, solver)
solver.SolverEquationsGet(solverEquations)
solverEquations.sparsityType = CMISS.SolverEquationsSparsityTypes.SPARSE
equationsSetIndex = solverEquations.EquationsSetAdd(equationsSet)
problem.SolverEquationsCreateFinish()

# Create boundary conditions
boundaryConditions = CMISS.BoundaryConditions()
solverEquations.BoundaryConditionsCreateStart(boundaryConditions)
boundaryConditions.NeumannSparsityTypeSet(CMISS.BoundaryConditionSparsityTypes.SPARSE) 
nodes = CMISS.Nodes()
region.NodesGet(nodes)
tol = 1.0e-12
for node in xrange(1, nodes.numberOfNodes + 1):
    nodeDomain = decomposition.NodeDomainGet(node, 1)
    if nodeDomain != computationalNodeNumber:
        continue
    # get x, y and z positions at node
    position = [
        geometricField.ParameterSetGetNodeDP(
            CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,
            1, 1, node, component + 1)
        for component in range(numberOfXi)]
    if (abs(position[0] - width) < tol or abs(position[0]) < tol):
        boundaryConditions.SetNode(dependentField,
                CMISS.FieldVariableTypes.DELUDELN, 1, 1, node, 1,
               CMISS.BoundaryConditionsTypes.NEUMANN_INTEGRATED_ONLY, 0.0)

    elif (abs(position[1] - height) < tol or abs(position[1]) < tol):
        # Set Neumann condition of 0 at boundaries
        boundaryConditions.SetNode(dependentField,
                CMISS.FieldVariableTypes.DELUDELN, 1, 1, node, 1,
               CMISS.BoundaryConditionsTypes.NEUMANN_INTEGRATED_ONLY, 0.0)

    elif (numberGlobalZElements > 0 and abs(position[2] - length) < tol or abs(position[2]) < tol):
        # Set Neumann condition of 0 at boundaries
        boundaryConditions.SetNode(dependentField,
                CMISS.FieldVariableTypes.DELUDELN, 1, 1, node, 1,
               CMISS.BoundaryConditionsTypes.NEUMANN_INTEGRATED_ONLY, 0.0)

solverEquations.BoundaryConditionsCreateFinish()
# Solve the problem
problem.Solve()

# Export results
#baseName = "diffusion"
#dataFormat = "PLAIN_TEXT"
#fml = CMISS.FieldMLIO()
#fml.OutputCreate(mesh, "", baseName, dataFormat)
#fml.OutputAddFieldNoType(baseName+".geometric", dataFormat, geometricField,
#    CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES)
#fml.OutputAddFieldNoType(baseName+".u", dataFormat, dependentField,
#    CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES)
#fml.OutputWrite("DiffusionExample.xml")
#fml.Finalise()

fields = CMISS.Fields()
fields.CreateRegion(region)
fields.NodesExport("MonodomainTP","FORTRAN")
fields.ElementsExport("MonodomainTP","FORTRAN")
fields.Finalise()
CMISS.Finalise()

############################################################################################

# Show the result in cmgui
#os.system("cmgui-wx visualise.com")










