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
	OperateOnEachCell(&m.lattice, func(f float64) {
		print(f)
	})
}

func (m *Model) assign(value float64) error {
	var err error
	if 0 < value && value < math.Pi*2 {
		err = errors.New("given value is out of bounds")
	}
	OperateOnEachCell[float64](&m.lattice, func(f float64) {
		f = value
	})

	return err
}
func (m *Model) randomize() {
	OperateOnEachCell[float64](&m.lattice, func(f float64) {
		f = rand.Float64() * 2 * math.Pi
	})
}

func (m *Model) energy() [][]float64 {
	var res = CreateMatrix[float64](m.size)

	for i := 0; i < m.size; i++ {
		for j := 0; j < m.size; j++ {
			res[i][j] = -m.coupling*(math.Cos(m.lattice[i][j]-m.lattice[i][(j+1)%m.size])) + math.Cos(m.lattice[i][j]-m.lattice[(i+1)%m.size][j]) - m.muB*m.magField*math.Cos(m.lattice[i][j])
		}
	}
	return res
}

/* This cannot be correct ngl */
func (m *Model) deltaEnergy(newAngles [][]float64) [][]float64 {
	res := CreateMatrix[float64](m.size)
	rows, cols := len(m.lattice), len(m.lattice[0])
	dE := 0.0

	for i := 0; i < rows; i++ {
		for j := 0; j < cols; j++ {
			dE += -m.coupling * (math.Cos(newAngles[i][j]-m.lattice[(i+1)%rows][j]) - math.Cos(m.lattice[i][j]-m.lattice[(i+1)%rows][j]))
			dE += -m.coupling * (math.Cos(newAngles[i][j]-m.lattice[i][(j+1)%cols]) - math.Cos(m.lattice[i][j]-m.lattice[i][(j+1)%cols]))
			dE += -m.coupling * (math.Cos(newAngles[i][j]-m.lattice[(i-1+rows)%rows][j]) - math.Cos(m.lattice[i][j]-m.lattice[(i-1+rows)%rows][j]))
			dE += -m.coupling * (math.Cos(newAngles[i][j]-m.lattice[i][(j-1+cols)%cols]) - math.Cos(m.lattice[i][j]-m.lattice[i][(j-1+cols)%cols]))
			dE += -m.muB * m.magField * (math.Cos(newAngles[i][j]) - math.Cos(m.lattice[i][j]))
		}
	}
	return res
}

func (m *Model) magnetization() float64 {
	var magX float64 = 0
	var magY float64 = 0
	OperateOnEachCell[float64](&m.lattice, func(f float64) {
		magX += math.Cos(f)
		magY += math.Sin(f)
	})
	return math.Sqrt(math.Pow(magX, 2) + math.Pow(magY, 2))
}

func xyModel() {
	h := newModel()
	h.printValues()

}
