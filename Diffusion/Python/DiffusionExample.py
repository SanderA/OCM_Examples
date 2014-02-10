#!/usr/bin/env python

#> \file
#> \author Sander Arens
#> \brief This is an example script to solve a Laplace problem using openCMISS calls in python.
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
width = 2.0
length = 3.0

numberGlobalXElements = 10
numberGlobalYElements = 10
numberGlobalZElements = 0

CMISS.DiagnosticsSetOn(CMISS.DiagnosticTypes.IN, [1,2,3,4,5], "Diagnostics", ["DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE"])

# Get the computational nodes information
numberOfComputationalNodes = CMISS.ComputationalNumberOfNodesGet()
computationalNodeNumber = CMISS.ComputationalNodeNumberGet()

# Creation a RC coordinate system
coordinateSystem = CMISS.CoordinateSystem()
coordinateSystemUserNumber = 1
coordinateSystem.CreateStart(coordinateSystemUserNumber)
if numberGlobalZElements == 0:
	coordinateSystem.DimensionSet(2)
else: 
	coordinateSystem.DimensionSet(3)
coordinateSystem.CreateFinish()

# Create a region
region = CMISS.Region()
regionUserNumber = 1
region.CreateStart(regionUserNumber, CMISS.WorldRegion)
region.LabelSet("DiffusionRegion")
region.CoordinateSystemSet(coordinateSystem)
region.CreateFinish()

# Create a tri-linear lagrange basis
basis = CMISS.Basis()
basisUserNumber = 1
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
generatedMeshUserNumber = 1
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
meshUserNumber = 1
generatedMesh.CreateFinish(meshUserNumber, mesh)

# Create a decomposition for the mesh
decomposition = CMISS.Decomposition()
decompositionUserNumber = 1
decomposition.CreateStart(decompositionUserNumber, mesh)
decomposition.TypeSet(CMISS.DecompositionTypes.CALCULATED)
decomposition.NumberOfDomainsSet(numberOfComputationalNodes)
decomposition.CreateFinish()

# Create a field for the geometry
geometricField = CMISS.Field()
geometricFieldUserNumber = 1
geometricField.CreateStart(geometricFieldUserNumber, region)
geometricField.MeshDecompositionSet(decomposition)
geometricField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U, 1, 1)
geometricField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U, 2, 1)
if numberGlobalZElements > 0:
	geometricField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U, 3, 1)
geometricField.CreateFinish()

# Set geometry from the generated mesh
generatedMesh.GeometricParametersCalculate(geometricField)

# Create Diffusion equations set
equationsSetField = CMISS.Field()
equationsSetFieldUserNumber = 2
equationsSet = CMISS.EquationsSet()
equationsSetUserNumber = 1
equationsSet.CreateStart(equationsSetUserNumber, region, geometricField,
        CMISS.EquationsSetClasses.CLASSICAL_FIELD,
        CMISS.EquationsSetTypes.DIFFUSION_EQUATION,
        CMISS.EquationsSetSubtypes.NO_SOURCE_DIFFUSION,
        equationsSetFieldUserNumber, equationsSetField)
equationsSet.CreateFinish()

# Create dependent field
dependentField = CMISS.Field()
dependentFieldUserNumber = 3
equationsSet.DependentCreateStart(dependentFieldUserNumber, dependentField)
equationsSet.DependentCreateFinish()

# Create the equations set material field variables
materialsField = CMISS.Field()
materialsFieldUserNumber = 4
equationsSet.MaterialsCreateStart(materialsFieldUserNumber, materialsField)
equationsSet.MaterialsCreateFinish()

# Create the equations set analytic field variables
analyticField = CMISS.Field()
analyticFieldUserNumber = 5
if numberGlobalZElements == 0:
	equationsSet.AnalyticCreateStart(CMISS.DiffusionAnalyticFunctionTypes.DIFFUSION_EQUATION_TWO_DIM_1, analyticFieldUserNumber, analyticField)
else: 
	print "Three dimensional analytical field not implemented"
equationsSet.AnalyticCreateFinish()


# Create equations
equations = CMISS.Equations()
equationsSet.EquationsCreateStart(equations)
equations.SparsityTypeSet(CMISS.EquationsSparsityTypes.SPARSE)
equations.OutputTypeSet(CMISS.EquationsOutputTypes.NONE)
equationsSet.EquationsCreateFinish()

# Create Diffusion problem
problem = CMISS.Problem()
problemUserNumber = 1
problem.CreateStart(problemUserNumber)
problem.SpecificationSet(CMISS.ProblemClasses.CLASSICAL_FIELD,
        CMISS.ProblemTypes.DIFFUSION_EQUATION,
        CMISS.ProblemSubTypes.NO_SOURCE_DIFFUSION)
problem.CreateFinish()

# Create control loops
controlLoop = CMISS.ControlLoop()
problem.ControlLoopCreateStart()
problem.ControlLoopGet([CMISS.ControlLoopIdentifiers.NODE], controlLoop)
controlLoop.OutputTypeSet(CMISS.ControlLoopOutputTypes.PROGRESS)
controlLoop.TimesSet(0.0, 1.001, 0.001)
controlLoop.TimeOutputSet(10)
problem.ControlLoopCreateFinish()

# Create problem solver
solver = CMISS.Solver()
linearSolver = CMISS.Solver()
problem.SolversCreateStart()
problem.SolverGet([CMISS.ControlLoopIdentifiers.NODE], 1, solver)
solver.OutputTypeSet(CMISS.SolverOutputTypes.NONE)
solver.DynamicLinearSolverGet(linearSolver)
linearSolver.LinearIterativeMaximumIterationsSet(1000)
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
solverEquations.BoundaryConditionsAnalytic()
solverEquations.BoundaryConditionsCreateFinish()

# Solve the problem
problem.Solve()

# Output Analytic analysis
CMISS.AnalyticAnalysisOutput(dependentField, "DiffusionAnalytics_x4_y4_q_T1")

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
#fields.NodesExport("Diffusion","FORTRAN")
fields.ElementsExport("Diffusion","FORTRAN")
fields.Finalise()

CMISS.Finalise()

############################################################################################

# Show the result in cmgui
#os.system("cmgui-wx visualise.com")










