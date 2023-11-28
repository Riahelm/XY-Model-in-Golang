package main

import (
	"bufio"
	"errors"
	"fmt"
	"math"
	"math/rand"
	"os"
)

type Model struct {
	size        int
	coupling    float64
	muB         float64
	temperature float64
	magField    float64
	lattice     [][]float64
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

	m.lattice = make([][]float64, m.size)
	for i := range m.lattice {
		m.lattice[i] = make([]float64, m.size)
	}
	println(m.size)
}

func (m *Model) printValues() {
	println("Size: ", m.size)
	println("Coupling: ", m.coupling)
	println("muB: ", m.muB)
	println("Temperature: ", m.temperature)
	println("MagField: ", m.magField)
	OperateOnEachCell(&m.lattice, func(f *float64) {
		print(f)
	})
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
	if 0 < value && value < math.Pi*2 {
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

func (m *Model) energy() [][]float64 {
	var res = CreateMatrix[float64](m.size)
	OperateOnCellsWithIndex[float64](&res, func(cell *float64, i int, j int) {
		*cell = -m.coupling*(math.Cos(m.lattice[i][j]-m.lattice[i][(j+1)%m.size])) + math.Cos(m.lattice[i][j]-m.lattice[(i+1)%m.size][j]) - m.muB*m.magField*math.Cos(m.lattice[i][j])
	})
	return res
}

/* This cannot be correct ngl */
func (m *Model) deltaEnergy(newAngles [][]float64) [][]float64 {
	res := CreateMatrix[float64](m.size)
	rows, cols := len(m.lattice), len(m.lattice[0])

	OperateOnCellsWithIndex(&res, func(cell *float64, i int, j int) {
		*cell += -m.coupling * (math.Cos(newAngles[i][j]-m.lattice[(i+1)%rows][j]) - math.Cos(m.lattice[i][j]-m.lattice[(i+1)%rows][j]))
		*cell += -m.coupling * (math.Cos(newAngles[i][j]-m.lattice[i][(j+1)%cols]) - math.Cos(m.lattice[i][j]-m.lattice[i][(j+1)%cols]))
		*cell += -m.coupling * (math.Cos(newAngles[i][j]-m.lattice[(i-1+rows)%rows][j]) - math.Cos(m.lattice[i][j]-m.lattice[(i-1+rows)%rows][j]))
		*cell += -m.coupling * (math.Cos(newAngles[i][j]-m.lattice[i][(j-1+cols)%cols]) - math.Cos(m.lattice[i][j]-m.lattice[i][(j-1+cols)%cols]))
		*cell += -m.muB * m.magField * (math.Cos(newAngles[i][j]) - math.Cos(m.lattice[i][j]))
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

func (m *Model) evolve(therm int, fraction float64, meas int, drop int) {
	var age = 0
	var iMeas = 0
	var newAngles = CreateMatrix[float64](m.size)
	var prob = CreateMatrix[float64](m.size)
	var trans = CreateMatrix[float64](m.size)
	var en = CreateMatrix[float64](m.size)
	var en2 = CreateMatrix[float64](m.size)
	var magPerSite float64 = 0
	var mag float64 = 0
	var mag2 float64 = 0
	for i := 0; i < therm+meas; i++ {
		OperateOnEachCell(&newAngles, func(f *float64) {
			*f = rand.NormFloat64()
			*f *= 2 * math.Pi
		})

		prob = m.transitionProbability(newAngles)

		OperateOnCellsWithIndex(&trans, func(cell *float64, i int, j int) {
			transition := rand.Float64() < prob[i][j] && rand.Float64() < fraction
			if transition {
				m.lattice[i][j] = newAngles[i][j]
			}
		})

		age++

		if i > therm && i%drop == 0 {
			tmpE := CreateMatrix[float64](m.size)
			var tmpM float64
			tmpE = m.energy()
			tmpM = m.magnetization()

			AddPMatrices(&en, tmpE)
			MultiplyPMatrices(&tmpE, tmpE)
			AddPMatrices(&en2, tmpE)
			magPerSite += tmpM / math.Pow(float64(m.size), 2)
			mag += tmpM
			mag2 += math.Pow(tmpM, 2)
			iMeas++
		}
	}

	DividePMatrixByScalar(&en, float64(iMeas))
	DividePMatrixByScalar(&en2, float64(iMeas))
	magPerSite /= float64(iMeas)
	mag /= float64(iMeas)
	mag2 /= float64(iMeas)

	energeticVariance := CreateMatrix[float64](m.size)
	energeticVariance = *SubtractMatrices(en2, *MultiplyMatrices(en, en))

	magneticVariance := mag2 - mag*mag

	specificHeat := 1 / (math.Pow(float64(m.size), 2) * math.Pow(m.temperature, 2)) * SumWholeArray(energeticVariance)
	magSusPerSite := 1 / (math.Pow(float64(m.size), 2) * m.temperature) * magneticVariance

	println(specificHeat, magSusPerSite)
}
func xyModel() {
	h := newModel()
	h.randomize()
	h.printLattice()
}
