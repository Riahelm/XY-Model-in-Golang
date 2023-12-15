package main

import (
	"errors"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"math"
	"math/rand"
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
	m.lattice = CreateMat[float64](size)
	return m
}

func (m *Model) assign(value float64) error {
	var err error
	if value < -math.Pi*2 || math.Pi*2 < value {
		err = errors.New("given value is out of bounds")
	} else {
		OpOnElems[float64](&m.lattice, func(f *float64) {
			*f = value
		})
	}
	return err
}
func (m *Model) randomize() {
	OpOnElems[float64](&m.lattice, func(f *float64) {
		*f = rand.NormFloat64()
		*f = math.Mod(*f, 2*math.Pi)
	})
}

func (m *Model) energy() float64 {

	res := -m.coupling

	res *= SumMat(AddMat(getMCos(SubMat(m.lattice, RollMatrix(m.lattice, 1, 1))), getMCos(SubMat(m.lattice, RollMatrix(m.lattice, 1, 0)))))

	res -= m.muB * m.magField * SumMat(getMCos(m.lattice))

	return res
}

func (m *Model) deltaEnergy(newAngles [][]float64) [][]float64 {
	res := CreateMat[float64](m.size)

	diff1 := SubMat(getMCos(SubMat(newAngles, RollMatrix(m.lattice, 1, 0))), getMCos(SubMat(m.lattice, RollMatrix(m.lattice, 1, 0))))

	diff2 := SubMat(getMCos(SubMat(newAngles, RollMatrix(m.lattice, 1, 1))), getMCos(SubMat(m.lattice, RollMatrix(m.lattice, 1, 1))))

	diff3 := SubMat(getMCos(SubMat(newAngles, RollMatrix(m.lattice, -1, 0))), getMCos(SubMat(m.lattice, RollMatrix(m.lattice, -1, 0))))

	diff4 := SubMat(getMCos(SubMat(newAngles, RollMatrix(m.lattice, -1, 1))), getMCos(SubMat(m.lattice, RollMatrix(m.lattice, -1, 1))))

	MulPMatByScalar(&diff1, -m.coupling)
	MulPMatByScalar(&diff2, -m.coupling)
	MulPMatByScalar(&diff3, -m.coupling)
	MulPMatByScalar(&diff4, -m.coupling)

	res = CopyMat(diff1)
	AddPMat(&res, diff2)
	AddPMat(&res, diff3)
	AddPMat(&res, diff4)
	AddPMat(&res, MulMatByScalar(SubMat(getMCos(newAngles), getMCos(m.lattice)), -m.muB*m.magField))

	return res
}

func (m *Model) magnetization() float64 {
	var magX float64 = 0
	var magY float64 = 0
	OpOnElems(&m.lattice, func(f *float64) {
		magX += math.Cos(*f)
		magY += math.Sin(*f)
	})
	return math.Sqrt(math.Pow(magX, 2) + math.Pow(magY, 2))
}

func (m *Model) transitionProbability(newAngles [][]float64) [][]float64 {
	res := CreateMat[float64](m.size)

	deltaMatrix := m.deltaEnergy(newAngles)
	OpOnElemsWithIndex(&res, func(cell *float64, i int, j int) {
		exp := -1.0 / m.temperature * deltaMatrix[i][j]
		*cell = math.Exp(exp)
	})
	return res
}

func (m *Model) evolve(therm int, fraction float64, meas int, drop int) (float64, float64, float64) {
	var age = 0
	var iMeas = 0
	var newAngles = CreateMat[float64](m.size)
	var prob = CreateMat[float64](m.size)
	var en = 0.
	var en2 = 0.
	var magPerSite float64 = 0
	var mag float64 = 0
	var mag2 float64 = 0

	for i := 0; i < therm+meas; i++ {
		OpOnElemsWithIndex(&newAngles, func(cell *float64, i int, j int) {
			*cell = rand.NormFloat64() + m.lattice[i][j]
			*cell = math.Mod(*cell, 2*math.Pi)
		})

		prob = m.transitionProbability(newAngles)

		OpOnElemsWithIndex(&m.lattice, func(cell *float64, i int, j int) {
			if rand.Float64() < prob[i][j] && rand.Float64() < fraction {
				*cell = newAngles[i][j]
			}
		})

		age++

		if i > therm && i%drop == 0 {

			tmpE := m.energy()
			tmpM := m.magnetization()

			en += tmpE
			en2 += tmpE * tmpE
			magPerSite += tmpM / math.Pow(float64(m.size), 2)
			mag += tmpM
			mag2 += tmpM * tmpM
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

func getMCos[F AcceptedVals](matrix [][]F) [][]float64 {
	res := CreateMat[float64](len(matrix))
	OpOnElemsWithIndex(&res, func(cell *float64, i int, j int) {
		*cell = math.Cos(float64(matrix[i][j]))
	})

	return res
}

func extractColumn[F any](matrix [][]F, index int) []F {
	res := make([]F, len(matrix))
	for i, fs := range matrix {
		res[i] = fs[index]
	}
	return res
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

func xyModel() {
	nRealizations := 500
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
		m := newParameterModel(50, 1, 0.67, tVals[i], 0)
		go simulate(m, nTherm, 0.1, nMeasure, nDrop, &results[i], &wg)
	}

	wg.Wait()

	plotData("Heat Plot", tVals, extractColumn(results, 0))
	plotData("Mag Plot", tVals, extractColumn(results, 1))
	plotData("Mag Sus Plot", tVals, extractColumn(results, 2))
}
