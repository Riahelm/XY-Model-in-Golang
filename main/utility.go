package main

import "golang.org/x/exp/constraints"

type AcceptedVals interface {
	constraints.Integer | constraints.Float | constraints.Complex
}

func CreateMatrix[F any](size int) [][]F {
	var res = make([][]F, size)
	for side := range res {
		res[side] = make([]F, size)
	}
	return res
}

func CopyMatrix[F any](m [][]F) [][]F {
	res := CreateMatrix[F](len(m))
	OperateOnCellsWithIndex(&res, func(cell *F, i int, j int) {
		*cell = m[i][j]
	})
	return res
}

func Roll[F any](slice [][]F, shiftX int, shiftY int) [][]F {
	length := len(slice)
	x, y := (shiftX%length)+length, (shiftY%length)+length
	res := CreateMatrix[F](length)
	OperateOnCellsWithIndex(&res, func(cell *F, i int, j int) {
		res[i][j] = slice[(i+x)%length][(j+y)%length]
	})
	return res
}

func OperateOnEachCellWithCopy[F any](matrix [][]F, doOperation func(cell *F)) [][]F {
	res := CopyMatrix(matrix)
	OperateOnEachCell(&res, func(cell *F) {
		doOperation(cell)
	})
	return res
}
func OperateOnEachCell[F any](matrix *[][]F, doOperation func(cell *F)) {
	for _, fs := range *matrix {
		for i := range fs {
			doOperation(&fs[i])
		}
	}
}

func OperateOnCellsWithIndex[F any](matrix *[][]F, doOperation func(cell *F, i int, j int)) {
	for i, fs := range *matrix {
		for j := range fs {
			doOperation(&fs[j], i, j)
		}
	}
}

func AddMatrixByScalar[F AcceptedVals](matrix [][]F, scalar F) [][]F {
	res := CopyMatrix(matrix)
	OperateOnCellsWithIndex[F](&res, func(cell *F, i int, j int) {
		*cell += scalar
	})
	return res
}

func AddPMatrixByScalar[F AcceptedVals](matrix *[][]F, scalar F) {
	OperateOnCellsWithIndex[F](matrix, func(cell *F, i int, j int) {
		*cell += scalar
	})
}

func AddMatrices[F AcceptedVals](matrix [][]F, matrix2 [][]F) [][]F {
	res := CopyMatrix(matrix)
	OperateOnCellsWithIndex[F](&res, func(cell *F, i int, j int) {
		*cell += matrix2[i][j]
	})
	return res
}

func AddPMatrices[F AcceptedVals](matrix *[][]F, matrix2 [][]F) {
	OperateOnCellsWithIndex[F](matrix, func(cell *F, i int, j int) {
		*cell += matrix2[i][j]
	})
}

func SubtractMatrixByScalar[F AcceptedVals](matrix [][]F, scalar F) [][]F {
	res := CopyMatrix(matrix)
	OperateOnCellsWithIndex[F](&res, func(cell *F, i int, j int) {
		*cell -= scalar
	})
	return res
}

func SubtractPMatrixByScalar[F AcceptedVals](matrix *[][]F, scalar F) {
	OperateOnCellsWithIndex[F](matrix, func(cell *F, i int, j int) {
		*cell -= scalar
	})
}

func SubtractMatrices[F AcceptedVals](matrix [][]F, matrix2 [][]F) [][]F {
	res := CopyMatrix(matrix)
	OperateOnCellsWithIndex[F](&res, func(cell *F, i int, j int) {
		*cell -= matrix2[i][j]
	})
	return res
}

func SubtractPMatrices[F AcceptedVals](matrix *[][]F, matrix2 [][]F) {
	OperateOnCellsWithIndex[F](matrix, func(cell *F, i int, j int) {
		*cell -= matrix2[i][j]
	})
}

func MultiplyMatrixByScalar[F AcceptedVals](matrix [][]F, scalar F) [][]F {
	res := CopyMatrix(matrix)
	OperateOnCellsWithIndex[F](&res, func(cell *F, i int, j int) {
		*cell *= scalar
	})
	return res
}

func MultiplyPMatrixByScalar[F AcceptedVals](matrix *[][]F, scalar F) {
	OperateOnCellsWithIndex[F](matrix, func(cell *F, i int, j int) {
		*cell *= scalar
	})
}

func MultiplyMatrices[F AcceptedVals](matrix [][]F, matrix2 [][]F) [][]F {
	res := CopyMatrix(matrix)
	OperateOnCellsWithIndex[F](&res, func(cell *F, i int, j int) {
		*cell *= matrix2[i][j]
	})
	return res
}

func MultiplyPMatrices[F AcceptedVals](matrix *[][]F, matrix2 [][]F) {
	OperateOnCellsWithIndex[F](matrix, func(cell *F, i int, j int) {
		*cell *= matrix2[i][j]
	})
}

func DivideMatrixByScalar[F AcceptedVals](matrix [][]F, scalar F) [][]F {
	res := CopyMatrix(matrix)
	OperateOnCellsWithIndex[F](&res, func(cell *F, i int, j int) {
		*cell /= scalar
	})
	return res
}

func DividePMatrixByScalar[F AcceptedVals](matrix *[][]F, scalar F) {
	OperateOnCellsWithIndex[F](matrix, func(cell *F, i int, j int) {
		*cell /= scalar
	})
}

func DivideMatrices[F AcceptedVals](matrix [][]F, matrix2 [][]F) [][]F {
	res := CopyMatrix(matrix)
	OperateOnCellsWithIndex[F](&res, func(cell *F, i int, j int) {
		*cell /= matrix2[i][j]
	})
	return res
}

func DividePMatrices[F AcceptedVals](matrix *[][]F, matrix2 [][]F) {
	OperateOnCellsWithIndex[F](matrix, func(cell *F, i int, j int) {
		*cell /= matrix2[i][j]
	})
}

func SumMatrix[F AcceptedVals](matrix [][]F) F {
	var res F
	OperateOnEachCell(&matrix, func(f *F) {
		res += *f
	})
	return res
}
