package main

func CreateMatrix[F float64 | float32](size int) [][]F {
	var res = make([][]F, size)
	for side := range res {
		res[side] = make([]F, size)
	}
	return res
}

func OperateOnEachCell[F float64 | float32](matrix *[][]F, doOperation func(F)) {
	for _, fs := range *matrix {
		for i := range fs {
			doOperation(fs[i])
		}
	}
}
