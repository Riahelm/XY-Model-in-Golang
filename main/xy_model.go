package main

import (
	"bufio"
	"errors"
	"fmt"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"math"
	"math/rand"
	"os"
	"slices"
	"sync"
)

type Model struct {
	size        int
	coupling    float64
	muB         float64
	temperature float64
	magField    float64
	lattice     [][]float64
}

func newParameterModel(size int, coupling float64, muB float64, temperature float64, magField float64) *Model {
	var m = new(Model)
	m.size = size
	m.coupling = coupling
	m.muB = muB
	m.temperature = temperature
	m.magField = magField
	m.lattice = CreateMatrix[float64](size)
	return m
}

func newModel() *Model {
	var m = new(Model)
	m.setValues()
	return m
}

func (m *Model) setValues() {
	var valCount = 0
	var err error = nil

	for err != nil || valCount != 5 {
		valCount, err = fmt.Scanf("%d %e %e %e %e",
			&m.size,
			&m.coupling,
			&m.muB,
			&m.temperature,
			&m.magField)
		r := bufio.NewReader(os.Stdin)
		_, _ = r.ReadString('\n')
		if err != nil || valCount != 5 {
			panic("Error handling the input")
		}
	}
	m.lattice = CreateMatrix[float64](m.size)
}

func (m *Model) printValues() {
	println("Size: ", m.size)
	println("Coupling: ", m.coupling)
	println("muB: ", m.muB)
	println("Temperature: ", m.temperature)
	println("MagField: ", m.magField)
	m.printLattice()
}

func (m *Model) printLattice() {
	OperateOnCellsWithIndex(&m.lattice, func(cell *float64, i int, j int) {
		print(*cell)
		if j == len(m.lattice[i])-1 {
			println()
		}
	})
}

func (m *Model) assign(value float64) error {
	var err error
	if value < -math.Pi*2 || math.Pi*2 < value {
		err = errors.New("given value is out of bounds")
	}
	OperateOnEachCell[float64](&m.lattice, func(f *float64) {
		*f = value
	})

	return err
}
func (m *Model) randomize() {
	OperateOnEachCell[float64](&m.lattice, func(f *float64) {
		*f = rand.Float64() * 2 * math.Pi
	})
}

func (m *Model) energy() float64 {
	var res = CreateMatrix[float64](m.size)
	OperateOnCellsWithIndex(&res, func(cell *float64, i int, j int) {
		neighbour := math.Cos(m.lattice[i][j]-m.lattice[i][(j+1)%m.size]) + math.Cos(m.lattice[i][j]-m.lattice[(i+1)%m.size][j])
		res[i][j] = -m.coupling * neighbour
	})
	SubtractPMatrixByScalar(&res, m.muB*m.magField*SumMatrix(OperateOnEachCellWithCopy(m.lattice, func(cell *float64) {
		*cell = math.Cos(*cell)
	})))
	return SumMatrix(res)
}

/* This cannot be correct ngl */
func (m *Model) deltaEnergy(newAngles [][]float64) [][]float64 {
	res := CreateMatrix[float64](m.size)

	diff1 := SubtractMatrices(
		OperateOnEachCellWithCopy(SubtractMatrices(newAngles, Roll(m.lattice, 1, 0)), func(cell *float64) {
			*cell = math.Cos(*cell)
		}), OperateOnEachCellWithCopy(SubtractMatrices(m.lattice, Roll(m.lattice, 1, 0)), func(cell *float64) {
			*cell = math.Cos(*cell)
		}))
	diff2 := SubtractMatrices(
		OperateOnEachCellWithCopy(SubtractMatrices(newAngles, Roll(m.lattice, 1, 1)), func(cell *float64) {
			*cell = math.Cos(*cell)
		}), OperateOnEachCellWithCopy(SubtractMatrices(m.lattice, Roll(m.lattice, 1, 1)), func(cell *float64) {
			*cell = math.Cos(*cell)
		}))
	diff3 := SubtractMatrices(
		OperateOnEachCellWithCopy(SubtractMatrices(newAngles, Roll(m.lattice, -1, 0)), func(cell *float64) {
			*cell = math.Cos(*cell)
		}), OperateOnEachCellWithCopy(SubtractMatrices(m.lattice, Roll(m.lattice, -1, 0)), func(cell *float64) {
			*cell = math.Cos(*cell)
		}))
	diff4 := SubtractMatrices(
		OperateOnEachCellWithCopy(SubtractMatrices(newAngles, Roll(m.lattice, -1, 1)), func(cell *float64) {
			*cell = math.Cos(*cell)
		}), OperateOnEachCellWithCopy(SubtractMatrices(m.lattice, Roll(m.lattice, -1, 1)), func(cell *float64) {
			*cell = math.Cos(*cell)
		}))

	MultiplyPMatrixByScalar(&diff1, -m.coupling)
	MultiplyPMatrixByScalar(&diff2, -m.coupling)
	MultiplyPMatrixByScalar(&diff3, -m.coupling)
	MultiplyPMatrixByScalar(&diff4, -m.coupling)

	res = CopyMatrix(diff1)
	AddPMatrices(&res, diff2)
	AddPMatrices(&res, diff3)
	AddPMatrices(&res, diff4)
	OperateOnCellsWithIndex(&res, func(cell *float64, i int, j int) {
		*cell += -m.coupling*m.magField*math.Cos(newAngles[i][j]) - math.Cos(m.lattice[i][j])
	})
	return res
}

func (m *Model) magnetization() float64 {
	var magX float64 = 0
	var magY float64 = 0
	OperateOnEachCell(&m.lattice, func(f *float64) {
		magX += math.Cos(*f)
		magY += math.Sin(*f)
	})
	return math.Sqrt(math.Pow(magX, 2) + math.Pow(magY, 2))
}

func (m *Model) transitionProbability(newAngles [][]float64) [][]float64 {
	res := CreateMatrix[float64](m.size)

	deltaMatrix := m.deltaEnergy(newAngles)
	OperateOnCellsWithIndex(&res, func(cell *float64, i int, j int) {
		*cell = math.Exp(-1.0 / m.temperature * deltaMatrix[i][j])
	})
	return res
}

func (m *Model) evolve(therm int, fraction float64, meas int, drop int) (float64, float64, float64) {
	var age = 0
	var iMeas = 0
	var newAngles = CreateMatrix[float64](m.size)
	var prob = CreateMatrix[float64](m.size)
	var en = 0.
	var en2 = 0.
	var magPerSite float64 = 0
	var mag float64 = 0
	var mag2 float64 = 0

	for i := 0; i < therm+meas; i++ {
		OperateOnCellsWithIndex(&newAngles, func(cell *float64, i int, j int) {
			*cell = rand.NormFloat64() + m.lattice[i][j]
			*cell = math.Mod(newAngles[i][j], 2*math.Pi)
		})

		prob = m.transitionProbability(newAngles)

		OperateOnCellsWithIndex(&m.lattice, func(cell *float64, i int, j int) {
			if rand.Float64() < prob[i][j] && rand.Float64() < fraction {
				*cell = newAngles[i][j]
			}
		})

		age++

		if i > therm && i%drop == 0 {
			var tmpE = 0.
			var tmpM float64
			tmpE = m.energy()
			tmpM = m.magnetization()

			en += tmpE
			en2 += tmpE * tmpE
			magPerSite += tmpM / math.Pow(float64(m.size), 2)
			mag += tmpM
			mag2 += math.Pow(tmpM, 2)
			iMeas++
		}
	}

	en /= float64(iMeas)
	en2 /= float64(iMeas)
	magPerSite /= float64(iMeas)
	mag /= float64(iMeas)
	mag2 /= float64(iMeas)

	energeticVariance := en2 - en*en
	magneticVariance := mag2 - mag*mag

	specificHeat := (1 / (math.Pow(float64(m.size), 2) * math.Pow(m.temperature, 2))) * energeticVariance
	magSusPerSite := (1 / (math.Pow(float64(m.size), 2) * m.temperature)) * magneticVariance

	return specificHeat, magPerSite, magSusPerSite
}

func simulate(m *Model, nTherm int, fraction float64, nMeas int, nDrop int, data *[]float64, group *sync.WaitGroup) {
	defer group.Done()
	m.randomize()
	(*data)[0], (*data)[1], (*data)[2] = m.evolve(nTherm, fraction, nMeas, nDrop)
}

func xyModel() {
	nRealizations := 100
	nTherm := 2000
	nMeasure := 2000
	nDrop := 5
	tMin := 0.01
	tMax := 2.
	var wg sync.WaitGroup
	results := make([][]float64, nRealizations)
	for i := range results {
		results[i] = make([]float64, 3)
	}
	tVals := make([]float64, nRealizations)
	for i := range tVals {
		tVals[i] = rand.Float64() * ((tMax - tMin) + tMin)
	}
	slices.Sort(tVals)
	for i := 0; i < nRealizations; i++ {
		wg.Add(1)
		go simulate(newParameterModel(30, 1, 0.67, tVals[i], 0), nTherm, 0.1, nMeasure, nDrop, &results[i], &wg)
	}

	wg.Wait()

	plotData("Heat Plot", tVals, extractColumn(results, 0))
	plotData("Mag Plot", tVals, extractColumn(results, 1))
	plotData("Mag Sus Plot", tVals, extractColumn(results, 2))
}

func plotData(plotName string, temps []float64, results []float64) {
	newPlot := plot.New()
	newPlot.Title.Text = plotName
	newPlot.X.Label.Text = "Temperature"
	newPlot.Y.Label.Text = plotName

	data := make(plotter.XYs, len(temps))
	for i := range data {
		data[i].X = temps[i]
		data[i].Y = results[i]
	}

	s, _ := plotter.NewScatter(data)

	newPlot.Add(s)
	newPlot.Save(1920, 1080, plotName+".png")
}

func extractColumn[F any](matrix [][]F, index int) []F {
	res := make([]F, len(matrix))
	for i, fs := range matrix {
		res[i] = fs[index]
	}
	return res
}
