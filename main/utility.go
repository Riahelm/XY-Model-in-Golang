package main

import (
	"fmt"
	"golang.org/x/exp/constraints"
)

type AcceptedVals interface {
	constraints.Integer | constraints.Float
}

func CreateMat[F any](size int) [][]F {
	var res = make([][]F, size)
	for side := range res {
		res[side] = make([]F, size)
	}
	return res
}

func CopyMat[F any](m [][]F) [][]F {
	res := CreateMat[F](len(m))
	OpOnElemsWithIndex(&res, func(cell *F, i int, j int) {
		*cell = m[i][j]
	})
	return res
}

func RollMatrix[F any](matrix [][]F, shift int, axis int) [][]F {
	rows, cols := len(matrix), len(matrix[0])
	result := make([][]F, rows)
	for i := range result {
		result[i] = make([]F, cols)
	}

	for i := 0; i < rows; i++ {
		for j := 0; j < cols; j++ {
			var newRow, newCol int
			switch axis {
			case 0:
				newRow = (i + shift + rows) % rows
				newCol = j
			case 1:
				newRow = i
				newCol = (j + shift + cols) % cols
			default:
				newRow = i
				newCol = j
			}
			result[newRow][newCol] = matrix[i][j]
		}
	}

	return result
}

func PrintMat[F any](matrix [][]F) {
	for i := range matrix {
		for j := range matrix[i] {
			fmt.Print("|", matrix[i][j], "\t")
		}
		fmt.Println("|")
	}
}

func OpOnElems[F any](matrix *[][]F, doOperation func(cell *F)) {
	for _, fs := range *matrix {
		for i := range fs {
			doOperation(&fs[i])
		}
	}
}

func OpOnElemsAndCopy[F any](matrix [][]F, doOperation func(cell *F)) [][]F {
	res := CopyMat(matrix)
	for i := range res {
		for j := range res[i] {
			doOperation(&res[i][j])
		}
	}
	return res
}

func OpOnElemsWithIndex[F any](matrix *[][]F, doOperation func(cell *F, i int, j int)) {
	for i, fs := range *matrix {
		for j := range fs {
			doOperation(&fs[j], i, j)
		}
	}
}

func AddMatByScalar[F AcceptedVals](matrix [][]F, scalar F) [][]F {
	res := CopyMat(matrix)
	OpOnElemsWithIndex[F](&res, func(cell *F, i int, j int) {
		*cell += scalar
	})
	return res
}

func AddPMatByScalar[F AcceptedVals](matrix *[][]F, scalar F) {
	OpOnElemsWithIndex[F](matrix, func(cell *F, i int, j int) {
		*cell += scalar
	})
}

func AddMat[F AcceptedVals](matrix [][]F, matrix2 [][]F) [][]F {
	res := CopyMat(matrix)
	OpOnElemsWithIndex[F](&res, func(cell *F, i int, j int) {
		*cell += matrix2[i][j]
	})
	return res
}

func AddPMat[F AcceptedVals](matrix *[][]F, matrix2 [][]F) {
	OpOnElemsWithIndex[F](matrix, func(cell *F, i int, j int) {
		*cell += matrix2[i][j]
	})
}

func SubMatByScalar[F AcceptedVals](matrix [][]F, scalar F) [][]F {
	res := CopyMat(matrix)
	OpOnElemsWithIndex[F](&res, func(cell *F, i int, j int) {
		*cell -= scalar
	})
	return res
}

func SubPMatByScalar[F AcceptedVals](matrix *[][]F, scalar F) {
	OpOnElemsWithIndex[F](matrix, func(cell *F, i int, j int) {
		*cell -= scalar
	})
}

func SubMat[F AcceptedVals](matrix [][]F, matrix2 [][]F) [][]F {
	res := CopyMat(matrix)
	OpOnElemsWithIndex[F](&res, func(cell *F, i int, j int) {
		*cell -= matrix2[i][j]
	})
	return res
}

func SubPMat[F AcceptedVals](matrix *[][]F, matrix2 [][]F) {
	OpOnElemsWithIndex[F](matrix, func(cell *F, i int, j int) {
		*cell -= matrix2[i][j]
	})
}

func MulMatByScalar[F AcceptedVals](matrix [][]F, scalar F) [][]F {
	res := CopyMat(matrix)
	OpOnElemsWithIndex[F](&res, func(cell *F, i int, j int) {
		*cell *= scalar
	})
	return res
}

func MulPMatByScalar[F AcceptedVals](matrix *[][]F, scalar F) {
	OpOnElemsWithIndex[F](matrix, func(cell *F, i int, j int) {
		*cell *= scalar
	})
}

func MulMatElems[F AcceptedVals](matrix [][]F, matrix2 [][]F) [][]F {
	res := CopyMat(matrix)
	OpOnElemsWithIndex[F](&res, func(cell *F, i int, j int) {
		*cell *= matrix2[i][j]
	})
	return res
}

func MulPMatElems[F AcceptedVals](matrix *[][]F, matrix2 [][]F) {
	OpOnElemsWithIndex[F](matrix, func(cell *F, i int, j int) {
		*cell *= matrix2[i][j]
	})
}

func DivMatElemsByScalar[F AcceptedVals](matrix [][]F, scalar F) [][]F {
	res := CopyMat(matrix)
	if scalar == 0 {
		panic("Division by zero error")
	}

	OpOnElemsWithIndex[F](&res, func(cell *F, i int, j int) {
		*cell /= scalar
	})
	return res
}

func DivPMatElemsByScalar[F AcceptedVals](matrix *[][]F, scalar F) {
	if scalar == 0 {
		panic("Division by zero error")
	}

	OpOnElemsWithIndex[F](matrix, func(cell *F, i int, j int) {
		*cell /= scalar
	})
}

func DivMatElems[F AcceptedVals](matrix [][]F, matrix2 [][]F) [][]F {
	res := CopyMat(matrix)
	OpOnElemsWithIndex[F](&res, func(cell *F, i int, j int) {
		if matrix2[i][j] == 0 {
			panic("Division by zero error")
		}
		*cell /= matrix2[i][j]
	})
	return res
}

func DivPMatElems[F AcceptedVals](matrix *[][]F, matrix2 [][]F) {
	OpOnElemsWithIndex[F](matrix, func(cell *F, i int, j int) {
		if matrix2[i][j] == 0 {
			panic("Division by zero error")
		}
		*cell /= matrix2[i][j]
	})
}

func SumMat[F AcceptedVals](matrix [][]F) F {
	var res F
	OpOnElems(&matrix, func(f *F) {
		res += *f
	})
	return res
}
