package main

func CreateMatrix[F float64 | float32](size int) [][]F {
	var res = make([][]F, size)
	for side := range res {
		res[side] = make([]F, size)
	}
	return res
}

func OperateOnEachCell[F float64 | float32](matrix *[][]F, doOperation func(cell *F)) {
	for _, fs := range *matrix {
		for i := range fs {
			doOperation(&fs[i])
		}
	}
}

func OperateOnCellsWithIndex[F float64 | float32](matrix *[][]F, doOperation func(cell *F, i int, j int)) {
	for i, fs := range *matrix {
		for j := range fs {
			doOperation(&fs[j], i, j)
		}
	}
}

func AddMatrixByScalar[F float64 | float32](matrix [][]F, scalar F) *[][]F {
	res := CreateMatrix[F](len(matrix))
	OperateOnCellsWithIndex[F](&res, func(cell *F, i int, j int) {
		*cell += scalar
	})
	return &res
}

func AddPMatrixByScalar[F float64 | float32](matrix *[][]F, scalar F) {
	OperateOnCellsWithIndex[F](matrix, func(cell *F, i int, j int) {
		*cell += scalar
	})
}

func AddMatrices[F float64 | float32](matrix [][]F, matrix2 [][]F) *[][]F {
	res := CreateMatrix[F](len(matrix))
	OperateOnCellsWithIndex[F](&res, func(cell *F, i int, j int) {
		*cell += matrix2[i][j]
	})
	return &res
}

func AddPMatrices[F float64 | float32](matrix *[][]F, matrix2 [][]F) {
	OperateOnCellsWithIndex[F](matrix, func(cell *F, i int, j int) {
		*cell += matrix2[i][j]
	})
}

func SubtractMatrices[F float64 | float32](matrix [][]F, matrix2 [][]F) *[][]F {
	res := CreateMatrix[F](len(matrix))
	OperateOnCellsWithIndex[F](&res, func(cell *F, i int, j int) {
		*cell -= matrix2[i][j]
	})
	return &res
}

func SubtractPMatrices[F float64 | float32](matrix *[][]F, matrix2 [][]F) {
	OperateOnCellsWithIndex[F](matrix, func(cell *F, i int, j int) {
		*cell -= matrix2[i][j]
	})
}

func MultiplyMatrixByScalar[F float64 | float32](matrix [][]F, scalar F) *[][]F {
	res := CreateMatrix[F](len(matrix))
	OperateOnCellsWithIndex[F](&res, func(cell *F, i int, j int) {
		*cell *= scalar
	})
	return &res
}

func MultiplyPMatrixByScalar[F float64 | float32](matrix *[][]F, scalar F) {
	OperateOnCellsWithIndex[F](matrix, func(cell *F, i int, j int) {
		*cell *= scalar
	})
}

func MultiplyMatrices[F float64 | float32](matrix [][]F, matrix2 [][]F) *[][]F {
	res := CreateMatrix[F](len(matrix))
	OperateOnCellsWithIndex[F](&res, func(cell *F, i int, j int) {
		*cell *= matrix2[i][j]
	})
	return &res
}

func MultiplyPMatrices[F float64 | float32](matrix *[][]F, matrix2 [][]F) {
	OperateOnCellsWithIndex[F](matrix, func(cell *F, i int, j int) {
		*cell *= matrix2[i][j]
	})
}

func DivideMatrixByScalar[F float64 | float32](matrix [][]F, scalar F) *[][]F {
	res := CreateMatrix[F](len(matrix))
	OperateOnCellsWithIndex[F](&res, func(cell *F, i int, j int) {
		*cell /= scalar
	})
	return &res
}

func DividePMatrixByScalar[F float64 | float32](matrix *[][]F, scalar F) {
	OperateOnCellsWithIndex[F](matrix, func(cell *F, i int, j int) {
		*cell /= scalar
	})
}

func DivideMatrices[F float64 | float32](matrix [][]F, matrix2 [][]F) *[][]F {
	res := CreateMatrix[F](len(matrix))
	OperateOnCellsWithIndex[F](&res, func(cell *F, i int, j int) {
		*cell /= matrix2[i][j]
	})
	return &res
}

func DividePMatrices[F float64 | float32](matrix *[][]F, matrix2 [][]F) {
	OperateOnCellsWithIndex[F](matrix, func(cell *F, i int, j int) {
		*cell /= matrix2[i][j]
	})
}

func SumWholeArray[F float64 | float32](matrix [][]F) F {
	var res F
	OperateOnEachCell(&matrix, func(f *F) {
		res += *f
	})
	return res
}
